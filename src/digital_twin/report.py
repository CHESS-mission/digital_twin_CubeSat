"""Definition of functions to produce the report after the simulation was run."""

import json

import astropy.units as u
import numpy as np

import pandas as pd
import os

from astropy.coordinates import (
    GCRS,
    ITRS,
    CartesianRepresentation,
    SphericalRepresentation,
)

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
)
from digital_twin.utils import (
    check_and_empty_folder,
)
import influxdb_client, os
from influxdb_client import InfluxDBClient, Point, WritePrecision
from influxdb_client.client.write_api import SYNCHRONOUS
import pandas as pd
from datetime import datetime, timedelta
from dotenv import load_dotenv

def produce_report(
    data: dict, report_params: dict, results_folder, env_file, verbose=False
) -> None:
    """Generates various plots based on the provided parameters and saves the report data in the specified folders.

    Args:
        data (dict): A dictionary containing simulation results such as times, battery, visibility, etc.
        report_params (dict): A dictionary specifying the report settings, including folder paths and figure preferences.
    """
    print("Saving the results...") if verbose else None
    folder = results_folder
    check_and_empty_folder(folder)
    figures_folder = folder + "figures/"
    check_and_empty_folder(figures_folder)
    data_folder = folder + "data/"
    data_folder_csv = folder + "csv/"
    check_and_empty_folder(data_folder)
    check_and_empty_folder(data_folder_csv)
    generate_figures(data, report_params["figures"], figures_folder, data_folder_csv)
    save_data(data, report_params["data"], data_folder)
    save_csv(data, report_params["data"], data_folder_csv)
    upload_to_influxdb(env_file, data_folder_csv)



def generate_figures(data: dict, figure_params: dict, folder: str, csv_folder:str) -> None:
    """Generate and save figures based on the simulation results."""

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
        save_filename = folder + "modes.png"
        plot_operating_modes(
            data["modes"],
            data["tofs"].to_value("second"),
            data["duration_sim"],
            save_filename=save_filename,
            show=False,
        )

    if figure_params["dashboard"] == "yes":
        save_filename = folder + "dashboard.png"
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
    step = np.max([int(len(data["tofs"]) / 100), 1])

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
            save_filename=folder + "battery_energy.png",
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
            save_filename=folder + "power_consumption.png",
            markersize_plot=0,
        )
    
    
    if figure_params["solar_cells_efficiency"] == "yes":
        plot_1d(
            data["tofs"].to_value("second")[1:],
            data["solar_cells_efficiency"][1:],
            "Solar cells efficiency Over Time",
            x_label,
            r"Solar cells efficiency",
            step=1,
            fill_under=False,
            remove_box=True,
            scatter=False,
            x_label_f=x_label_f,
            show=False,
            save_filename=folder + "solar_cells_efficiency.png",
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
            save_filename=folder + "power_generation.png",
            markersize_plot=0,
        )
    if figure_params["power_balance"] == "yes":
        balance = data["generation"][1:] - data["consumption"][1:]
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
            save_filename=folder + "power_balance.png",
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
            save_filename=folder + "data_storage.png",
        )
        plot_1d(
            data["tofs"].to_value("second"),
            data["storage_payload"],
            "Data Storage Over Time (GNSS and TOF/Camera)",
            x_label,
            r"GNSS and TOF/Camera Data Storage ($Mbit$)",
            step=step,
            fill_under=False,
            remove_box=True,
            scatter=False,
            x_label_f=x_label_f,
            show=False,
            save_filename=folder + "data_storage_payload.png",
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
            save_filename=folder + "data_storage_HK.png",
        )

    if figure_params["visibility_windows"] == "yes":
        save_filename = folder + "visibility_windows.png"
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
        save_filename = folder + "eclipse_windows.png"
        plot_boolean_bars(
            data["eclipse"],
            data["tofs"].to_value("second"),
            data["duration_sim"],
            "Eclipse Windows",
            save_filename=save_filename,
            show=False,
        )


def save_data(data: dict, data_params: dict, folder: str) -> None:
    """Save the simulation data to specified files."""

    if data_params["telecom_data"] == "yes":
        with open(folder + "times.npy", "wb") as f:
            np.save(f, data["tofs"].to_value("second"))
        with open(folder + "visibility.npy", "wb") as f:
            np.save(f, data["vis"])
        with open(folder + "data.npy", "wb") as f:
            np.save(f, data["storage"])
        with open(folder + "data_payload.npy", "wb") as f:
            np.save(f, data["storage_payload"])
        with open(folder + "data_HK.npy", "wb") as f:
            np.save(f, data["storage_HK"])

    if data_params["eps_data"] == "yes":
        with open(folder + "times.npy", "wb") as f:
            np.save(f, data["tofs"].to_value("second"))
        with open(folder + "battery.npy", "wb") as f:
            np.save(f, data["battery"])
        with open(folder + "consumption.npy", "wb") as f:
            np.save(f, data["consumption"])
        with open(folder + "generation.npy", "wb") as f:
            np.save(f, data["generation"])
        with open(folder + "eclipse.npy", "wb") as f:
            np.save(f, data["eclipse"])
        with open(folder + "solar_cells_efficiency.npy", "wb") as f:
            np.save(f, data["solar_cells_efficiency"])

    if data_params["modes"] == "yes":
        with open(folder + "times.npy", "wb") as f:
            np.save(f, data["tofs"].to_value("second"))
        with open(folder + "modes.npy", "wb") as f:
            np.save(f, data["modes"])

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

    if data_params["density"] == "yes":
        with open(folder + "times.npy", "wb") as f:
            np.save(f, data["tofs"].to_value("second"))
        with open(folder + "density.npy", "wb") as f:
            np.save(f, data["density_array"])



def save_csv(data: dict, data_params: dict, folder: str) -> None:
    """Save all simulation data into a single CSV file with labeled columns."""
    
    df_data = {}

    def to_1d(array):
        """Convert to NumPy array and ensure it's 1D."""
        return np.array(array).flatten()

    # telecom data
    df_data["times_telecom"] = to_1d(data["tofs"].to_value("second"))
    df_data["visibility"] = to_1d(data["vis"])
    df_data["data"] = to_1d(data["storage"])
    df_data["data_payload"] = to_1d(data["storage_payload"])
    df_data["data_HK"] = to_1d(data["storage_HK"])

    # eps data
    # df_data["times_eps"] = to_1d(data["tofs"].to_value("second"))
    df_data["battery"] = to_1d(data["battery"])
    df_data["consumption"] = to_1d(data["consumption"])
    df_data["generation"] = to_1d(data["generation"])
    df_data["eclipse"] = to_1d(data["eclipse"])
    df_data["solar_cells_efficiency"] = to_1d(data["solar_cells_efficiency"])

    # modes data
    # df_data["times_modes"] = to_1d(data["tofs"].to_value("second"))
    df_data["modes"] = to_1d(data["modes"])

    # orbital element data
    # df_data["times_orbital"] = to_1d(data["tofs"].to_value("second"))
    df_data["altitude"] = to_1d(data["altitudes"])
    df_data["RAAN"] = to_1d(data["RAANs"])
    df_data["AOP"] = to_1d(data["AOPs"])
    df_data["ECC"] = to_1d(data["ECCs"])
    df_data["INC"] = to_1d(data["INCs"])

    # spacecraft state data 
    df_data["times_density"] = to_1d(data["tofs"].to_value("second"))
    df_data["density"] = to_1d(data["density_array"])


    # add nan padding to ensure all columns have the same length
    max_length = max(len(v) for v in df_data.values())
    for key in df_data:
        current_length = len(df_data[key])
        if current_length < max_length:
            df_data[key] = np.pad(df_data[key], (0, max_length - current_length), constant_values=np.nan)

    df = pd.DataFrame(df_data)
    csv_filename = os.path.join(folder, "simulation_data.csv")
    df.to_csv(csv_filename, index=False)

    # trajectory coordinates
    raw_xyz = CartesianRepresentation(data["rr"], xyz_axis=-1)
    raw_obstime = data["epochs_array"]
    gcrs_xyz = GCRS(
        raw_xyz, obstime=raw_obstime, representation_type=CartesianRepresentation
    )
    itrs_xyz = gcrs_xyz.transform_to(ITRS(obstime=raw_obstime))  # Converts raw coordinates to ITRS ones.
    itrs_latlon = itrs_xyz.represent_as(SphericalRepresentation)
    
    # Convert to degrees
    latitudes = itrs_latlon.lat.to(u.deg).value
    longitudes = itrs_latlon.lon.to(u.deg).value

    # Create a DataFrame
    df = pd.DataFrame({
        "latitude_deg": latitudes,
        "longitude_deg": longitudes
    })

    # Save to CSV
    csv_filename = os.path.join(folder, "trajectory_coords.csv")
    df.to_csv(csv_filename, index=False)  
    return


def upload_to_influxdb(env_file:str, csv_folder: str) -> None:
    """Upload simulation data to InfluxDB from CSV files."""

    # Securely retrieve credentials from environment variables
    load_dotenv(env_file)
    token = os.environ.get("INFLUXDB_TOKEN")
    org = os.environ.get("INFLUXDB_ORG", "EST")  # Default to "EST" if not set
    url = os.environ.get("INFLUXDB_URL", "http://localhost:8086")  # Default URL

    if not token:
        print("Error: INFLUXDB_TOKEN not found in environment variables. Upload aborted.")
        return

    client = influxdb_client.InfluxDBClient(url=url, token=token, org=org)

    df = pd.read_csv(os.path.join(csv_folder, "simulation_data.csv"), delimiter=',') # all the data
    traj_df = pd.read_csv(os.path.join(csv_folder, "trajectory_coords.csv"), delimiter=',') # the trajectory data

    # we artificially add a timestamp to the data to convert to datetime
    now = datetime(2025, 1, 1, 0, 0, 0)
    df['times_telecom'] = df['times_telecom'].apply(lambda x: now + timedelta(seconds=x))

    # the full df
    df = pd.concat([df.reset_index(drop=True), traj_df.reset_index(drop=True)], axis=1)

    # Write data to InfluxDB
    bucket="NICE"
    write_api = client.write_api(write_options=SYNCHRONOUS)
    delete_api = client.delete_api()

    # Delete all previous data from bucket
    start = "1970-01-01T00:00:00Z"
    stop =  datetime(2070, 1, 1, 0, 0, 0)
    delete_api.delete(start, stop, '', bucket=bucket, org=org)

    # the subsystems' times are not included because for now they are the same for all
    # the .tag are used to index the data, and can be used for filtering
    # Pre-build list of points
    points = [
        Point("satellite_data")
            .tag("mode", int(row["modes"]))
            .tag("visible", int(row["visibility"]))
            .field("visibility", float(row["visibility"]))
            .field("data", float(row["data"]))
            .field("data_payload", float(row["data_payload"]))
            .field("data_HK", float(row["data_HK"]))
            .field("battery", float(row["battery"]))
            .field("consumption", float(row["consumption"]))
            .field("generation", float(row["generation"]))
            .field("eclipse", float(row["eclipse"]))
            .field("modes", float(row["modes"]))
            .field("altitude", float(row["altitude"]))
            .field("RAAN", float(row["RAAN"]))
            .field("AOP", float(row["AOP"]))
            .field("ECC", float(row["ECC"]))
            .field("INC", float(row["INC"]))
            .field("density", float(row["density"]))
            .field("Lat", float(row["latitude_deg"]))
            .field("Lng", float(row["longitude_deg"]))
            .field("solar_cells_efficiency", float(row["solar_cells_efficiency"]))
            .time(row["times_telecom"], write_precision=WritePrecision.NS)
        for _, row in df.iterrows()
    ]

    # Send all points at once
    write_api.write(bucket=bucket, org=org, record=points)
    print("Upload complete.")
    return