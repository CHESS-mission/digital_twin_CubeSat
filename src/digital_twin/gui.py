"""Simple GUI for controlling the simulation in real-time."""

import tkinter as tk
from tkinter import ttk
from typing import Optional
import threading


class SimulationGUI:
    """Basic GUI to send commands to the running simulation."""
    
    def __init__(self, simulation=None):
        """Initialize the GUI.
        
        Args:
            simulation: The Simulation instance to control
        """
        self.simulation = simulation
        self.root = None
        
    def _setup_ui(self):
        """Set up the user interface elements."""
        # Title
        title_label = tk.Label(
            self.root, 
            text="Command GUI", 
            font=("Arial", 16, "bold")
        )
        title_label.pack(pady=20)
        
        # Status frame
        status_frame = tk.LabelFrame(self.root, text="Current Status", padx=10, pady=10)
        status_frame.pack(padx=20, pady=10, fill="x")
        
        self.status_label = tk.Label(
            status_frame, 
            text="Simulation not connected", 
            font=("Arial", 10)
        )
        self.status_label.pack()
        
        # Start periodic status updates
        self._update_status_periodically()
        
        # Mode control frame
        mode_frame = tk.LabelFrame(self.root, text="Mode Control", padx=10, pady=10)
        mode_frame.pack(padx=20, pady=10, fill="both", expand=True)
        
        # Mode selection
        tk.Label(mode_frame, text="Select Operating Mode:", font=("Arial", 11)).pack(pady=5)
        
        # Mode buttons
        self.modes = {
            0: "IDLE/CHARGING",
            1: "MEASUREMENT",
            2: "SAFE",
            3: "UHF-COM",
            4: "X-BAND-COM"
        }
        
        button_frame = tk.Frame(mode_frame)
        button_frame.pack(pady=10)
        
        for mode_num, mode_name in self.modes.items():
            btn = tk.Button(
                button_frame,
                text=f"{mode_num}: {mode_name}",
                command=lambda m=mode_num: self.send_mode_command(m),
                width=20,
                height=2,
                font=("Arial", 10)
            )
            btn.pack(pady=5)
        
    def _update_status_periodically(self):
        """Update status display periodically."""
        try:
            current_mode = self.simulation.switch_algo.operating_mode
            mode_name = self.modes.get(current_mode, "Unknown")
            self.status_label.config(
                text=f"Current Mode: {current_mode} ({mode_name})",
                fg="green"
            )
        except Exception:
            self.status_label.config(
                text="Simulation not connected",
                fg="red"
            )
        # Schedule next update
        self.root.after(200, self._update_status_periodically)
    
    def send_mode_command(self, mode: int):
        """Send a set_mode command to the simulation.
        
        Args:
            mode: The mode number to set
        """        
        try:
            self.simulation.send_command("set_mode", {"mode": mode})
        except Exception as e:
            print(f"ERROR: {str(e)}")
    
    def run(self):
        """Internal method to create and run the GUI (called in the thread)."""
        self.root = tk.Tk()
        self.root.title("GUI")
        self.root.geometry("400x550")
        self._setup_ui()
        self.root.mainloop()
    
    def start_non_blocking(self):
        """Start the GUI in a separate thread (non-blocking)."""
        gui_thread = threading.Thread(target=self.run, daemon=True)
        gui_thread.start()
        return gui_thread
