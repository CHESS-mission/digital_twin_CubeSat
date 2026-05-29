"""HTTP API wrapper for controlled Digital Twin simulation stepping.

This file is intentionally separate from ``main.py``. The original entrypoint
still runs the upstream batch simulation. This wrapper exposes a singleton
real-time controlled run for bridge integration experiments.
"""

from __future__ import annotations

import argparse
import json
import threading
import time
from http import HTTPStatus
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from typing import Any
from urllib.parse import parse_qs, urlparse

from digital_twin.realtime_runner import (
    RealtimeSimulationRunner,
    SimulationFiles,
    build_simulation,
)


PROJECT_ROOT = Path(__file__).resolve().parents[1]


class SimulationService:
    """Thread-safe singleton service around a controlled DT runner."""

    def __init__(self, project_root: Path, history_limit: int = 500) -> None:
        self._project_root = project_root
        self._history_limit = history_limit
        self._lock = threading.RLock()
        self._runner: RealtimeSimulationRunner | None = None
        self._history: list[dict[str, Any]] = []
        self._state = "idle"
        self._error: str | None = None
        self._speed_multiplier = 1.0
        self._max_steps: int | None = None
        self._stop_requested = False
        self._paused = False
        self._paused_since: float | None = None
        self._paused_total_s = 0.0
        self._started_wall_s: float | None = None
        self._thread: threading.Thread | None = None

    def start(self, payload: dict[str, Any]) -> dict[str, Any]:
        with self._lock:
            if self._state in {"running", "paused"}:
                raise ValueError("simulation is already active")

            files = SimulationFiles(
                simulation=payload.get("simulation_file", SimulationFiles.simulation),
                orbit=payload.get("orbit_file", SimulationFiles.orbit),
                spacecraft=payload.get("spacecraft_file", SimulationFiles.spacecraft),
                ground_station=payload.get(
                    "ground_station_file", SimulationFiles.ground_station
                ),
                mission_design=payload.get(
                    "mission_design_file", SimulationFiles.mission_design
                ),
            )
            simulation = build_simulation(
                self._project_root,
                files,
                simulation_overrides=_dict_or_empty(
                    payload.get("simulation_overrides")
                ),
                disable_influxdb=bool(payload.get("disable_influxdb", True)),
                quiet=bool(payload.get("quiet", True)),
            )

            self._runner = RealtimeSimulationRunner(simulation)
            self._history = [self._runner.latest_snapshot.to_dict()]
            self._speed_multiplier = float(payload.get("speed_multiplier", 1.0))
            if self._speed_multiplier <= 0:
                raise ValueError("speed_multiplier must be > 0")
            self._max_steps = _optional_int(payload.get("max_steps"))
            self._stop_requested = False
            self._paused = not bool(payload.get("autostart", True))
            self._paused_since = time.monotonic() if self._paused else None
            self._paused_total_s = 0.0
            self._started_wall_s = time.monotonic()
            self._error = None
            self._state = "paused" if self._paused else "running"

            self._thread = threading.Thread(target=self._run_loop, daemon=True)
            self._thread.start()
            return self.status()

    def pause(self) -> dict[str, Any]:
        with self._lock:
            if self._state != "running":
                raise ValueError("simulation is not running")
            self._paused = True
            self._paused_since = time.monotonic()
            self._state = "paused"
            return self.status()

    def resume(self) -> dict[str, Any]:
        with self._lock:
            if self._state != "paused":
                raise ValueError("simulation is not paused")
            if self._paused_since is not None:
                self._paused_total_s += time.monotonic() - self._paused_since
            self._paused_since = None
            self._paused = False
            self._state = "running"
            return self.status()

    def stop(self) -> dict[str, Any]:
        thread: threading.Thread | None
        with self._lock:
            self._stop_requested = True
            self._paused = False
            self._state = "stopped" if self._runner else "idle"
            thread = self._thread

        if thread and thread.is_alive() and thread is not threading.current_thread():
            thread.join(timeout=2)
        return self.status()

    def step_once(self) -> dict[str, Any]:
        with self._lock:
            if self._runner is None:
                raise ValueError("simulation has not been started")
            if self._state not in {"paused", "running"}:
                raise ValueError(f"cannot step while state is {self._state}")
            snapshot = self._step_locked()
            if self._paused and self._paused_since is not None:
                self._paused_total_s += time.monotonic() - self._paused_since
                self._paused_since = time.monotonic()
            if self._state == "running":
                self._state = "paused"
                self._paused = True
                self._paused_since = time.monotonic()
            return snapshot

    def command(self, payload: dict[str, Any]) -> dict[str, Any]:
        command_name = payload.get("command")
        if not command_name:
            raise ValueError("command is required")
        params = payload.get("params", {})
        if not isinstance(params, dict):
            raise ValueError("params must be an object")

        with self._lock:
            if self._runner is None:
                raise ValueError("simulation has not been started")
            if str(command_name) in {"set_speed", "set_speed_multiplier"}:
                speed_multiplier = float(params.get("speed_multiplier", params.get("speed", 1.0)))
                self._set_speed_locked(speed_multiplier)
                return {
                    "accepted": True,
                    "command": command_name,
                    "params": {"speed_multiplier": self._speed_multiplier},
                }
            self._runner.send_command(str(command_name), params)
            return {"accepted": True, "command": command_name, "params": params}

    def status(self) -> dict[str, Any]:
        with self._lock:
            latest = self._history[-1] if self._history else None
            return {
                "state": self._state,
                "error": self._error,
                "speed_multiplier": self._speed_multiplier,
                "max_steps": self._max_steps,
                "step": latest["step"] if latest else None,
                "simulation_time_s": latest["simulation_time_s"] if latest else None,
                "delta_t_s": self._runner.delta_t_s if self._runner else None,
                "total_steps": self._runner.total_steps if self._runner else None,
                "latest_snapshot": latest,
            }

    def latest_snapshot(self) -> dict[str, Any]:
        with self._lock:
            if not self._history:
                raise ValueError("no snapshot available")
            return self._history[-1]

    def snapshots_since(self, since_step: int | None = None) -> dict[str, Any]:
        with self._lock:
            if since_step is None:
                snapshots = list(self._history)
            else:
                snapshots = [
                    snapshot
                    for snapshot in self._history
                    if int(snapshot["step"]) > since_step
                ]
            return {"snapshots": snapshots, "count": len(snapshots)}

    def _run_loop(self) -> None:
        while True:
            with self._lock:
                if self._stop_requested or self._runner is None:
                    return
                if self._state in {"completed", "failed", "stopped"}:
                    return
                if self._paused:
                    sleep_for = 0.1
                else:
                    sleep_for = self._seconds_until_next_step_locked()

            if sleep_for > 0:
                time.sleep(min(sleep_for, 0.25))
                continue

            with self._lock:
                if self._stop_requested or self._runner is None or self._paused:
                    continue
                self._step_locked()

    def _seconds_until_next_step_locked(self) -> float:
        assert self._runner is not None
        assert self._started_wall_s is not None

        next_step = self._runner.current_step + 1
        target = (
            self._started_wall_s
            + self._paused_total_s
            + (next_step * self._runner.delta_t_s / self._speed_multiplier)
        )
        return target - time.monotonic()

    def _set_speed_locked(self, speed_multiplier: float) -> None:
        if speed_multiplier <= 0:
            raise ValueError("speed_multiplier must be > 0")

        assert self._runner is not None
        now = time.monotonic()
        paused_total_s = self._paused_total_s
        self._speed_multiplier = speed_multiplier
        self._started_wall_s = (
            now
            - paused_total_s
            - (self._runner.current_step * self._runner.delta_t_s / self._speed_multiplier)
        )

    def _step_locked(self) -> dict[str, Any]:
        assert self._runner is not None

        snapshot = self._runner.step().to_dict()
        self._history.append(snapshot)
        if len(self._history) > self._history_limit:
            self._history = self._history[-self._history_limit :]

        if self._runner.status == "failed":
            self._state = "failed"
            self._error = self._runner.error
        elif self._max_steps is not None and self._runner.current_step >= self._max_steps:
            self._state = "completed"
        elif self._runner.status == "completed":
            self._state = "completed"
        elif not self._paused:
            self._state = "running"

        return snapshot


def _optional_int(value: Any) -> int | None:
    if value is None:
        return None
    parsed = int(value)
    if parsed < 0:
        raise ValueError("max_steps must be >= 0")
    return parsed


def _dict_or_empty(value: Any) -> dict[str, Any]:
    if value is None:
        return {}
    if not isinstance(value, dict):
        raise ValueError("simulation_overrides must be an object")
    return value


def make_handler(service: SimulationService) -> type[BaseHTTPRequestHandler]:
    class SimulationRequestHandler(BaseHTTPRequestHandler):
        server_version = "DigitalTwinRealtime/0.1"

        def do_GET(self) -> None:  # noqa: N802 - http.server API
            parsed = urlparse(self.path)
            try:
                if parsed.path == "/health":
                    self._send_json({"ok": True})
                elif parsed.path == "/status":
                    self._send_json(service.status())
                elif parsed.path == "/snapshot":
                    self._send_json(service.latest_snapshot())
                elif parsed.path == "/snapshots":
                    query = parse_qs(parsed.query)
                    since_values = query.get("since_step")
                    since_step = int(since_values[0]) if since_values else None
                    self._send_json(service.snapshots_since(since_step))
                else:
                    self._send_error(HTTPStatus.NOT_FOUND, "unknown endpoint")
            except ValueError as exc:
                self._send_error(HTTPStatus.BAD_REQUEST, str(exc))
            except Exception as exc:
                self._send_error(
                    HTTPStatus.INTERNAL_SERVER_ERROR,
                    f"{type(exc).__name__}: {exc}",
                )

        def do_POST(self) -> None:  # noqa: N802 - http.server API
            parsed = urlparse(self.path)
            try:
                payload = self._read_json()
                if parsed.path == "/start":
                    self._send_json(service.start(payload), HTTPStatus.CREATED)
                elif parsed.path == "/pause":
                    self._send_json(service.pause())
                elif parsed.path == "/resume":
                    self._send_json(service.resume())
                elif parsed.path == "/stop":
                    self._send_json(service.stop())
                elif parsed.path == "/step":
                    self._send_json(service.step_once())
                elif parsed.path == "/command":
                    self._send_json(service.command(payload), HTTPStatus.ACCEPTED)
                else:
                    self._send_error(HTTPStatus.NOT_FOUND, "unknown endpoint")
            except ValueError as exc:
                self._send_error(HTTPStatus.BAD_REQUEST, str(exc))
            except Exception as exc:
                self._send_error(
                    HTTPStatus.INTERNAL_SERVER_ERROR,
                    f"{type(exc).__name__}: {exc}",
                )

        def log_message(self, fmt: str, *args: Any) -> None:
            return

        def _read_json(self) -> dict[str, Any]:
            length = int(self.headers.get("Content-Length", "0"))
            if length == 0:
                return {}
            raw = self.rfile.read(length)
            decoded = json.loads(raw.decode("utf-8"))
            if not isinstance(decoded, dict):
                raise ValueError("request body must be a JSON object")
            return decoded

        def _send_json(
            self,
            payload: dict[str, Any],
            status: HTTPStatus = HTTPStatus.OK,
        ) -> None:
            body = json.dumps(payload).encode("utf-8")
            self.send_response(status)
            self.send_header("Content-Type", "application/json")
            self.send_header("Content-Length", str(len(body)))
            self.end_headers()
            self.wfile.write(body)

        def _send_error(self, status: HTTPStatus, message: str) -> None:
            self._send_json({"error": message}, status)

    return SimulationRequestHandler


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Digital Twin real-time API wrapper")
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", type=int, default=8765)
    parser.add_argument("--history-limit", type=int, default=500)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    service = SimulationService(PROJECT_ROOT, history_limit=args.history_limit)
    handler = make_handler(service)
    server = ThreadingHTTPServer((args.host, args.port), handler)
    print(f"Digital Twin real-time API listening on http://{args.host}:{args.port}")
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        service.stop()
        server.server_close()


if __name__ == "__main__":
    main()
