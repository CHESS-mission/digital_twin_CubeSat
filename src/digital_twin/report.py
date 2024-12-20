from typing import Dict
import json

import astropy.units as u
import numpy as np


from digital_twin.orbit_propagator.constants import attractor_string
from digital_twin.plotting import (
    plot_1d,
    plot_orbit_trajectory_3d,
    plot_orbit_2d,
    plot_groundtrack,
    plot_operating_modes,
    find_x_scale,
    plot_boolean_bars,
    plot_dashboard,
    plot_orbital_elem_evolution,
    plot_1d_multiple,
)
from digital_twin.utils import (
    check_and_empty_folder,
)


def produce_report(data: Dict, report_params: Dict) -> None:
    """Produce plots and other report data."""
    folder = report_params["folder"]
    check_and_empty_folder(folder)
    figures_folder = folder + "figures/"
    check_and_empty_folder(figures_folder)
    data_folder = folder + "data/"
    check_and_empty_folder(data_folder)
    generate_figures(data, report_params["figures"], figures_folder)
    save_data(data, report_params["data"], data_folder)


def generate_figures(data: Dict, figure_params: Dict, folder: str) -> None:
    if figure_params["orbital_elem_evolution"] == "yes":
        plot_orbital_elem_evolution(
            data["tofs"],
            data["RAANs"],
            data["AOPs"],
            data["ECCs"],
            data["INCs"],
            data["altitudes"],
            folder,
            data["duration_sim"],
        )
    if figure_params["trajectory_2d"] == "yes":
        plot_orbit_2d(
            figure_params["title_figures"],
            attractor_string,
            folder,
            "initial orbit",
            orbit=data["initial_orbit"],
        )
    if figure_params["trajectory_3d"] == "yes":
        plot_orbit_trajectory_3d(
            figure_params["title_figures"],
            attractor_string,
            folder,
            orbit=data["initial_orbit"],
            label_orbit="initial orbit",
            label_traj="trajectory",
            traj=data["rr"],
        )

    stations_coords = []
    stations_names = []
    stations_colors = []
    for station in data["ground_stations"]:
        name, pos, color = station.get_name_pos_color()
        stations_coords.append(pos)
        stations_names.append(name)
        stations_colors.append(color)
    if figure_params["groundtrack"] == "yes":
        plot_groundtrack(
            figure_params["title_figures"],
            data["rr"],
            data["epochs_array"],
            "CHESS_1 cubeSat",
            folder,
            stations_coords=np.array(stations_coords),
            stations_names=np.array(stations_names),
            stations_colors=np.array(stations_colors),
        )

    if figure_params["modes"] == "yes":
        save_filename = folder + "modes.pdf"
        plot_operating_modes(
            data["modes"],
            data["tofs"].to_value("second"),
            data["duration_sim"],
            save_filename=save_filename,
            show=False,
        )

    if figure_params["dashboard"] == "yes":
        save_filename = folder + "dashboard.pdf"
        if (data["vis"].shape)[1] == 1:
            plot_dashboard(
                data["modes"],
                data["eclipse"],
                data["vis"],
                data["tofs"].to_value("second"),
                data["duration_sim"],
                title="Operating Modes Over Time",
                save_filename=save_filename,
                show=False,
            )
        else:
            print(
                "WARNING: dashboard function not implemented for multiple ground stations"
            )

    x_label, x_label_f = find_x_scale(data["duration_sim"])
    step = int(len(data["tofs"]) / 100)

    if figure_params["battery_energy"] == "yes":
        plot_1d(
            data["tofs"].to_value("second"),
            (data["battery"] * (u.W * u.s)).to(u.W * u.h),
            "Battery Energy Over Time",
            x_label,
            r"Battery Energy ($Wh$)",
            step=1,
            fill_under=False,
            remove_box=True,
            scatter=False,
            x_label_f=x_label_f,
            show=False,
            save_filename=folder + "battery_energy.pdf",
            markersize_plot=0,
        )

    if figure_params["power_consumption"] == "yes":
        plot_1d(
            data["tofs"].to_value("second")[1:],
            data["consumption"][1:],
            "Power Consumption Over Time",
            x_label,
            r"Power Consumption ($W$)",
            step=1,  # no step because we want to capture small intervals of time
            fill_under=False,
            remove_box=True,
            scatter=False,
            x_label_f=x_label_f,
            show=False,
            save_filename=folder + "power_consumption.pdf",
            markersize_plot=0,
        )

    if figure_params["power_generation"] == "yes":
        plot_1d(
            data["tofs"].to_value("second")[1:],
            data["generation"][1:],
            "Power Generation Over Time",
            x_label,
            r"Power Generation ($W$)",
            step=1,
            fill_under=False,
            remove_box=True,
            scatter=False,
            x_label_f=x_label_f,
            show=False,
            save_filename=folder + "power_generation.pdf",
            markersize_plot=0,
        )
    if figure_params["power_balance"] == "yes":
        balance = data["generation"][1:] - data["consumption"][1:]
        # xs = [
        #     data["tofs"].to_value("second")[1:],
        #     data["tofs"].to_value("second")[1:],
        #     data["tofs"].to_value("second")[1:],
        # ]
        # ys = [-data["consumption"][1:], data["generation"][1:], balance]
        # colors = ["blue", "purple", "green"]
        # labels = ["Power Consumption", "Power Generation", "Power Balance"]
        # steps = [1, 1, 1]
        # plot_1d_multiple(
        #     xs,
        #     ys,
        #     "Power Balance Over Time",
        #     x_label,
        #     r"Power ($W$)",
        #     colors=colors,
        #     labels=labels,
        #     step=steps,
        #     fill_under=False,
        #     remove_box=True,
        #     scatter=True,
        #     show=False,
        #     x_label_f=x_label_f,
        #     save_filename=folder + "power_balance.png",
        # )
        plot_1d(
            data["tofs"].to_value("second")[1:],
            balance,
            "Power Balance Over Time",
            x_label,
            r"Power Balance ($W$)",
            step=1,
            fill_under=False,
            remove_box=True,
            scatter=False,
            x_label_f=x_label_f,
            show=False,
            save_filename=folder + "power_balance.pdf",
            markersize_plot=0,
        )

    if figure_params["data_storage"] == "yes":
        plot_1d(
            data["tofs"].to_value("second"),
            data["storage"],
            "Data Storage Over Time",
            x_label,
            r"Data Storage ($Mbit$)",
            step=step,
            fill_under=False,
            remove_box=True,
            scatter=False,
            x_label_f=x_label_f,
            show=False,
            save_filename=folder + "data_storage.pdf",
        )
        plot_1d(
            data["tofs"].to_value("second"),
            data["storage_GNSS_TOF"],
            "Data Storage Over Time (GNSS and TOF)",
            x_label,
            r"GNSS/TOF Data Storage ($Mbit$)",
            step=step,
            fill_under=False,
            remove_box=True,
            scatter=False,
            x_label_f=x_label_f,
            show=False,
            save_filename=folder + "data_storage_GNSS_TOF.pdf",
        )
        plot_1d(
            data["tofs"].to_value("second"),
            data["storage_HK"],
            "Data Storage Over Time (Housekeeping Data)",
            x_label,
            r"HK Data Storage ($Mbit$)",
            step=step,
            fill_under=False,
            remove_box=True,
            scatter=False,
            x_label_f=x_label_f,
            show=False,
            save_filename=folder + "data_storage_HK.pdf",
        )

    if figure_params["visibility_windows"] == "yes":
        save_filename = folder + "visibility_windows.pdf"
        if (data["vis"].shape)[1] == 1:
            plot_boolean_bars(
                data["vis"],
                data["tofs"].to_value("second"),
                data["duration_sim"],
                "Visibility Windows",
                save_filename=save_filename,
                show=False,
            )
        else:
            print(
                "WARNING: dashboard function not implemented for multiple ground stations"
            )

    if figure_params["eclipse_windows"] == "yes":
        save_filename = folder + "eclipse_windows.pdf"
        plot_boolean_bars(
            data["eclipse"],
            data["tofs"].to_value("second"),
            data["duration_sim"],
            "Eclipse Windows",
            save_filename=save_filename,
            show=False,
        )


def save_data(data: Dict, data_params: Dict, folder: str) -> None:
    if data_params["telecom_data"] == "yes":
        with open(folder + "times.npy", "wb") as f:
            np.save(f, data["tofs"].to_value("second"))
        with open(folder + "visibility.npy", "wb") as f:
            np.save(f, data["vis"])
        with open(folder + "modes.npy", "wb") as f:
            np.save(f, data["modes"])
        with open(folder + "data.npy", "wb") as f:
            np.save(f, data["storage"])
        with open(folder + "data_GNSS_TOF.npy", "wb") as f:
            np.save(f, data["storage_GNSS_TOF"])
        with open(folder + "data_HK.npy", "wb") as f:
            np.save(f, data["storage_HK"])

    if data_params["altitude_data"] == "yes":
        with open(folder + "times.npy", "wb") as f:
            np.save(f, data["tofs"].to_value("second"))
        with open(folder + "altitude.npy", "wb") as f:
            np.save(f, data["altitudes"])

    if data_params["orbital_element_data"] == "yes":
        with open(folder + "times.npy", "wb") as f:
            np.save(f, data["tofs"].to_value("second"))
        with open(folder + "altitude.npy", "wb") as f:
            np.save(f, data["altitudes"])
        with open(folder + "RAAN.npy", "wb") as f:
            np.save(f, data["RAANs"])
        with open(folder + "AOP.npy", "wb") as f:
            np.save(f, data["AOPs"])
        with open(folder + "ECC.npy", "wb") as f:
            np.save(f, data["ECCs"])
        with open(folder + "INC.npy", "wb") as f:
            np.save(f, data["INCs"])

    if data_params["eclipse_data"] == "yes":
        with open(folder + "times.npy", "wb") as f:
            np.save(f, data["tofs"].to_value("second"))
        with open(folder + "eclipse.npy", "wb") as f:
            np.save(f, data["eclipse"])

    if data_params["orbit_state"] == "yes":
        with open(folder + "orbit_state.json", "w") as f:
            json.dump(
                data["orbit_state"], f, indent=4
            )  # Save with indentation for readability

    if data_params["spacecraft_state"] == "yes":
        with open(folder + "spacecraft_state.json", "w") as f:
            json.dump(
                data["spacecraft_state"], f, indent=4
            )  # Save with indentation for readability
