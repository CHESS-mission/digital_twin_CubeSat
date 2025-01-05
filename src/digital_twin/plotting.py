"""Gather all functions relating to plotting."""

from typing import Callable, Any

from astropy import units as u
from astropy.coordinates import (
    GCRS,
    ITRS,
    CartesianRepresentation,
    SphericalRepresentation,
)
from astropy.units import Quantity
from matplotlib import pyplot as plt
from matplotlib.ticker import FuncFormatter
import matplotlib.dates as mdates
import numpy as np
import plotly.graph_objects as go
from poliastro.earth.plotting import GroundtrackPlotter
from poliastro.plotting import OrbitPlotter3D
from poliastro.twobody import Orbit

from digital_twin.constants import mode_dict


def plot_1d(
    x: np.ndarray,
    y: np.ndarray,
    title: str,
    xlabel: str,
    ylabel: str,
    step: int = 1,
    fill_under: bool = True,
    remove_box: bool = True,
    x_range: tuple = None,
    y_range: tuple = None,
    scatter: bool = False,
    x_label_f: Callable = None,
    custom_y_ticks: bool = False,
    y_ticks: np.ndarray = None,
    y_tick_labels: np.ndarray = None,
    custom_x_ticks: bool = False,
    x_ticks: np.ndarray = None,
    x_tick_labels: np.ndarray = None,
    save_filename: str = None,
    show: bool = True,
    markersize_plot: int = 4,
    date_x_axis: bool = False,
    date_interval: int = 1,
    date_format: str = "%d/%m/%Y",
) -> None:
    """Function for general plotting in 2d with x and y arrays as input.
    Many parameters are optional to change the visualization design."""
    # Downsample the x and y arrays
    x_downsampled = x[::step]
    y_downsampled = y[::step]

    # Create the plot
    plt.figure(figsize=(6, 4))
    if scatter:
        plt.scatter(x_downsampled, y_downsampled, marker=".", color="blue", s=10)
    else:
        plt.plot(
            x_downsampled,
            y_downsampled,
            marker=".",
            linestyle="-",
            color="blue",
            markersize=markersize_plot,
        )

    # Beautifying the plot
    plt.title(title, fontsize=13, fontweight="medium")
    plt.xlabel(xlabel, fontsize=11)
    plt.ylabel(ylabel, fontsize=11)
    plt.grid(True, which="both", linestyle="--", linewidth=0.5)
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)

    # Setting x and y axis ranges if provided
    if x_range is not None:
        plt.xlim(x_range)  # Set the x-axis range (tuple)
    if y_range is not None:
        plt.ylim(y_range)  # Set the y-axis range (tuple)

    # Adding customizations
    ax = plt.gca()  # Get the current axis
    if x_label_f is not None:
        ax.xaxis.set_major_formatter(FuncFormatter(x_label_f))
    if custom_y_ticks:
        if y_ticks is None or y_tick_labels is None:
            raise ValueError("Must provide y ticks and y tick labels")
        ax.set_yticks(y_ticks)
        ax.set_yticklabels(y_tick_labels)
    if custom_x_ticks:
        if x_ticks is None or x_tick_labels is None:
            raise ValueError("Must provide x ticks and x tick labels")
        ax.set_xticks(x_ticks)
        ax.set_xticklabels(x_tick_labels)
    if fill_under:
        plt.fill_between(x_downsampled, y_downsampled, color="lightblue", alpha=0.3)
    if remove_box:
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_visible(False)
        ax.spines["bottom"].set_visible(False)
        # Remove ticks
        ax.yaxis.set_ticks_position("none")
        ax.xaxis.set_ticks_position("none")
    else:
        ax.spines["top"].set_color((0.8, 0.8, 0.8))
        ax.spines["right"].set_color((0.8, 0.8, 0.8))
        ax.spines["left"].set_color((0.8, 0.8, 0.8))
        ax.spines["bottom"].set_color((0.8, 0.8, 0.8))

    # If want to have dates as the x axis
    if date_x_axis:
        # Set major locator to show every other day (or adjust as needed)
        ax.xaxis.set_major_locator(mdates.DayLocator(interval=date_interval))
        ax.xaxis.set_major_formatter(
            mdates.DateFormatter(date_format)
        )  # Display the year only
        plt.xticks(rotation=45)

    # Show the plot
    plt.tight_layout()
    if save_filename is not None:
        plt.savefig(save_filename, dpi=300, bbox_inches="tight")

    if show:
        plt.show()
    else:
        plt.close()


def plot_1d_multiple(
    x: list,
    y: list,
    title: str,
    xlabel: str,
    ylabel: str,
    colors: list,
    labels: list,
    step: list,
    fill_under: bool = True,
    remove_box: bool = True,
    x_range: tuple = None,
    y_range: tuple = None,
    scatter: bool = False,
    x_label_f: Callable = None,
    custom_y_ticks: bool = False,
    y_ticks: np.ndarray = None,
    y_tick_labels: np.ndarray = None,
    custom_x_ticks: bool = False,
    x_ticks: np.ndarray = None,
    x_tick_labels: np.ndarray = None,
    save_filename: str = None,
    show: bool = True,
    markersize_plot: int = 4,
    date_x_axis: bool = False,
    date_interval: int = 1,
    date_format: str = "%d/%m/%Y",
    log_y_scale: bool = False,
    legend_title: str = "",
) -> None:
    """Function for general plotting in 2d with x and y arrays as input for MULTIPLE PLOTS on the same ax."""
    nb_plots = len(y)  # = len(x)

    # Downsample the x and y arrays
    x_downsampled = []
    y_downsampled = []
    for i in range(nb_plots):
        x_downsampled.append(x[i][:: step[i]])
        y_downsampled.append(y[i][:: step[i]])

    # Create the plot
    plt.figure(figsize=(6, 4))
    if scatter:
        for i in range(nb_plots):
            plt.scatter(
                x_downsampled[i],
                y_downsampled[i],
                marker=".",
                color=colors[i],
                s=10,
                label=labels[i],
            )
    else:
        for i in range(nb_plots):
            plt.plot(
                x_downsampled[i],
                y_downsampled[i],
                marker=".",
                linestyle="-",
                color=colors[i],
                markersize=markersize_plot,
                label=labels[i],
            )

    # Beautifying the plot
    if log_y_scale:
        plt.yscale("log")
    plt.title(title, fontsize=13, fontweight="medium")
    plt.xlabel(xlabel, fontsize=11)
    plt.ylabel(ylabel, fontsize=11)
    plt.grid(True, which="both", linestyle="--", linewidth=0.5)
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)

    # Setting x and y axis ranges if provided
    if x_range is not None:
        plt.xlim(x_range)  # Set the x-axis range (tuple)
    if y_range is not None:
        plt.ylim(y_range)  # Set the y-axis range (tuple)

    # Adding customizations
    ax = plt.gca()  # Get the current axis
    if x_label_f is not None:
        ax.xaxis.set_major_formatter(FuncFormatter(x_label_f))
    if custom_y_ticks:
        if y_ticks is None or y_tick_labels is None:
            raise ValueError("Must provide y ticks and y tick labels")
        ax.set_yticks(y_ticks)
        ax.set_yticklabels(y_tick_labels)
    if custom_x_ticks:
        if x_ticks is None or x_tick_labels is None:
            raise ValueError("Must provide x ticks and x tick labels")
        ax.set_xticks(x_ticks)
        ax.set_xticklabels(x_tick_labels)
    if fill_under:
        for i in range(nb_plots):
            plt.fill_between(
                x_downsampled[i], y_downsampled[i], color="light" + colors[i], alpha=0.3
            )
    if remove_box:
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_visible(False)
        ax.spines["bottom"].set_visible(False)
        # Remove ticks
        ax.yaxis.set_ticks_position("none")
        ax.xaxis.set_ticks_position("none")
    else:
        ax.spines["top"].set_color((0.8, 0.8, 0.8))
        ax.spines["right"].set_color((0.8, 0.8, 0.8))
        ax.spines["left"].set_color((0.8, 0.8, 0.8))
        ax.spines["bottom"].set_color((0.8, 0.8, 0.8))

    # if want to have dates as the x axis
    if date_x_axis:
        # Set major locator to show every other day (or adjust as needed)
        ax.xaxis.set_major_locator(mdates.DayLocator(interval=date_interval))
        ax.xaxis.set_major_formatter(
            mdates.DateFormatter(date_format)
        )  # display the year only
        plt.xticks(rotation=45)

    # plt.legend()
    legend = plt.legend(title=legend_title)
    # legend.get_title().set_fontsize(12)

    # Show the plot
    plt.tight_layout()
    if save_filename is not None:
        plt.savefig(save_filename, dpi=300, bbox_inches="tight")

    if show:
        plt.show()
    else:
        plt.close()


def seconds_to_days(x: Any, pos: Any) -> str:
    """Convert seconds to days and return a string cast."""
    return f"{x / 86400:.0f}"


def seconds_to_hours(x: Any, pos: Any) -> str:
    """Convert seconds to hours and return a string cast."""
    return f"{x / 3600:.0f}"


def seconds_to_minutes(x: Any, pos: Any) -> str:
    """Convert seconds to minutes and return a string cast."""
    return f"{x / 60:.0f}"


def plot_orbit_trajectory_3d(
    title: str,
    attractor_name: str,
    folder: str,
    orbit: Orbit = None,
    label_orbit: str = "",
    label_traj: str = "",
    traj: np.ndarray = None,
    option: int = 1,
    orbit_color: str = "red",
    orbit_width: float = 7,
    orbit_type: str = "solid",
    traj_color: str = "black",
    traj_width: float = 0.1,
    traj_type: str = "solid",
    attractor_color: str = "mediumblue",
    sat_color: str = "red",
    show: bool = False,
) -> None:
    """Plot an orbit in 3d around its attractor.

    Option 1: with a grid and tick values but no tick lines.
    Option 2: with no grid not tick values.
    """

    if (orbit is None) and (traj is None):
        raise AttributeError("Need to provide at least 1 orbit or trajectory to plot!")

    plotter = OrbitPlotter3D()
    if orbit is not None:
        fig = plotter.plot(orbit)
    if traj is not None:
        coords = CartesianRepresentation(
            traj[:, 0] << u.km, traj[:, 1] << u.km, traj[:, 2] << u.km
        )
        fig = plotter.plot_trajectory(coords, label=label_traj)

    if option == 1:
        scene_dic = dict(
            xaxis=dict(
                showgrid=True,
                gridcolor="grey",
                zerolinecolor="grey",
                showbackground=True,
                title="<i>X</i> (km)",
                ticks="",
                showticklabels=True,
            ),
            yaxis=dict(
                showgrid=True,
                gridcolor="grey",
                zerolinecolor="grey",
                showbackground=True,
                title="<i>Y</i> (km)",
                ticks="",
                showticklabels=True,
            ),
            zaxis=dict(
                showgrid=True,
                gridcolor="grey",
                zerolinecolor="grey",
                showbackground=True,
                title="<i>Z</i> (km)",
                ticks="",
                showticklabels=True,
            ),
            camera=dict(eye=dict(x=1.8, y=-1.3, z=0.2)),
        )
    if option == 2:
        scene_dic = dict(
            xaxis=dict(
                backgroundcolor="rgba(105,132,163, 1)",
                showgrid=False,
                zeroline=False,
                showbackground=True,
                title="X",
                ticks="",
                showticklabels=False,
            ),
            yaxis=dict(
                backgroundcolor="rgba(138,164,197, 1)",
                showgrid=False,
                zeroline=False,
                showbackground=True,
                title="Y",
                ticks="",
                showticklabels=False,
            ),
            zaxis=dict(
                backgroundcolor="rgba(171,196,231, 1)",
                showgrid=False,
                zeroline=False,
                showbackground=True,
                title="Z",
                ticks="",
                showticklabels=False,
            ),
            camera=dict(eye=dict(x=1.8, y=-1.3, z=0.2)),
        )
    if option == 3:
        scene_dic = dict(
            xaxis=dict(
                visible=False,  # Hide the x-axis
            ),
            yaxis=dict(
                visible=False,  # Hide the y-axis
            ),
            zaxis=dict(
                visible=False,  # Hide the z-axis
            ),
            camera=dict(eye=dict(x=1.5, y=-1, z=0.1)),
        )

    fig.update_layout(
        title={
            "text": title,
            "y": 0.92,
            "x": 0.5,
            "xanchor": "center",
            "yanchor": "top",
        },
        margin=dict(l=0, r=0, b=0, t=0),  # Set margins to zero
        height=400,
        width=450,
        scene=scene_dic,  # depending on option chosen
        legend=dict(
            title="Legend:",
            orientation="h",
            yanchor="bottom",
            y=-0.1,
            xanchor="center",
            x=0.5,
        ),
    )

    for trace in fig.data:
        if isinstance(trace, go.Scatter3d):
            if trace.name == label_traj:
                trace.line.color = traj_color
                trace.line.width = traj_width
                trace.line.dash = traj_type
                trace.showlegend = True
            else:
                trace.line.color = orbit_color
                trace.line.width = orbit_width
                trace.line.dash = orbit_type
                trace.name = label_orbit
                trace.showlegend = True
        elif isinstance(trace, go.Surface):
            if trace.name == attractor_name:
                trace.colorscale = [[0, attractor_color], [1, attractor_color]]
            else:
                trace.colorscale = [[0, sat_color], [1, sat_color]]

    fig.write_html(folder + "orbit_3d.html")
    if show:
        fig.show()


def plot_orbit_2d(
    title: str,
    folder: str,
    label: str,
    orbit: Orbit,
    orbit_color: str = "black",
    orbit_width: float = 2,
    orbit_type: str = "solid",
    attractor_color: str = "mediumblue",
    sat_color: str = "red",
    show: bool = False,
) -> None:
    """Plot an Orbit object in 2d (orbit plane) around its attractor."""

    fig = orbit.plot(interactive=True)

    fig.update_layout(
        title={
            "text": title,
            "y": 0.95,
            "x": 0.5,
            "xanchor": "center",
            "yanchor": "top",
        },
        margin=dict(l=0, r=0, b=60, t=60),
        height=400,
        width=450,
        legend=dict(
            title="Legend:",
            orientation="h",
            yanchor="bottom",
            y=-0.3,
            xanchor="center",
            x=0.5,
        ),
        xaxis=dict(
            title="<i>X</i> (km)",
            zerolinewidth=1,
            range=[-8000, 8000],
        ),
        yaxis=dict(
            title="<i>Y</i> (km)",
            zerolinewidth=1,
            range=[-8000, 8000],
        ),
    )

    fig.data[0].line.dash = orbit_type
    fig.data[0].line.color = orbit_color
    fig.data[0].line.width = orbit_width
    fig.data[0].showlegend = True
    fig.data[0].name = label
    fig.data[1].marker.color = sat_color

    # Update the central circle properties
    fig.layout.shapes[1].fillcolor = attractor_color
    fig.layout.shapes[1].line.color = attractor_color

    fig.write_image(folder + "orbit_2d.png", scale=3)
    if show:
        fig.show()


def plot_groundtrack(
    title: str,
    rr: np.ndarray,
    epochs: np.ndarray,
    label: str,
    folder: str,
    traj_color: str = "purple",
    traj_width: float = 1,
    stations_coords: np.ndarray = None,
    stations_names: np.ndarray = None,
    stations_colors: str = "red",
    show: bool = False,
) -> None:
    """Plot a satellite groundtrack on Earth."""
    # For building geo traces
    gp = GroundtrackPlotter()

    # Add satellite groundtrack
    raw_xyz = CartesianRepresentation(rr, xyz_axis=-1)
    raw_obstime = epochs
    gcrs_xyz = GCRS(
        raw_xyz, obstime=raw_obstime, representation_type=CartesianRepresentation
    )
    itrs_xyz = gcrs_xyz.transform_to(
        ITRS(obstime=raw_obstime)
    )  # Converts raw coordinates to ITRS ones.
    itrs_latlon = itrs_xyz.represent_as(SphericalRepresentation)
    gp.add_trace(
        go.Scattergeo(
            lat=itrs_latlon.lat.to(u.deg),
            lon=itrs_latlon.lon.to(u.deg),
            mode="lines",
            name=label,
            line={"color": traj_color, "width": traj_width},
        )
    )

    # Add ground station
    if stations_coords is not None:
        for index, station in enumerate(stations_coords):
            coord = station * u.deg  # [LAT LON]
            name = stations_names[index]
            gp.add_trace(
                go.Scattergeo(
                    lat=coord[0],
                    lon=coord[-1],
                    name=name,
                    marker={"color": stations_colors[index]},
                    showlegend=True,
                )
            )

    # Customize plots
    gp.update_layout(
        title={
            "text": "Groundtrack: " + title,
            "y": 0.92,
            "x": 0.5,
            "xanchor": "center",
            "yanchor": "top",
        },
        margin=dict(l=10, r=10, b=0, t=50),  # Set margins to zero
        height=500,
        width=800,
        legend=dict(
            title="Legend:",
            orientation="h",
            yanchor="bottom",
            y=-0,
            xanchor="center",
            x=0.5,
        ),
    )
    gp.fig.write_image(folder + "groundtrack.png", scale=3)
    if show:
        gp.fig.show()


def plot_operating_modes(
    modes: np.ndarray,
    tofs: np.ndarray,
    duration_sim: Quantity["time"],
    save_filename: str = None,
    show: bool = True,
) -> None:
    """Plot operating modes during simulation to visualize them throughout time."""

    modes = [int(mode) for mode in modes]

    # update what the x_axis scale should be
    x_label, x_label_f = find_x_scale(duration_sim)

    # Mode names
    mode_labels = [mode_dict[key] for key in sorted(mode_dict.keys())]
    colors = ["#5cb85c", "#d9534f", "#5bc0de", "#f0ad4e", "#a0522d", "#6f42c1"]

    plt.figure(figsize=(6, 3.5))
    plt.ylim(-0.6, len(mode_labels) - 0.3)  # Adjust limits to avoid space

    # Manually add horizontal grid lines at each mode level to make grid
    for y in range(len(mode_labels)):
        plt.axhline(
            y=y, color="gray", linestyle="--", linewidth=0.8
        )  # Add horizontal grid lines

    # Create the step plot using horizontal lines (hlines)
    for i in range(1, len(tofs)):
        plt.hlines(
            modes[i - 1],
            tofs[i - 1],
            tofs[i],
            colors=colors[modes[i - 1]],
            linewidth=10,
        )
    length_unit = duration_sim.value / len(modes)
    plt.hlines(
        modes[-1],
        tofs[-1],
        tofs[-1] + length_unit,
        colors=colors[modes[-1]],
        linewidth=10,
    )  # Extend last line

    plt.title("Satellite Modes Over Time", fontsize=13, fontweight="medium")
    plt.xlabel(x_label, fontsize=11)

    # Add a legend
    handles = [
        plt.Line2D([0], [0], color=colors[i], linewidth=5)
        for i in range(len(mode_labels))
    ]
    plt.legend(
        handles,
        mode_labels,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.5),
        ncol=3,
        title="Modes",
    )

    # Hide y-ticks and y-axis label
    plt.yticks([])
    plt.ylabel("")

    ax = plt.gca()
    ax.xaxis.set_ticks_position("none")  # Remove x-ticks
    ax.spines["top"].set_color((0.8, 0.8, 0.8))
    ax.spines["right"].set_color((0.8, 0.8, 0.8))
    ax.spines["left"].set_color((0.8, 0.8, 0.8))
    ax.spines["bottom"].set_color((0.8, 0.8, 0.8))
    ax.xaxis.set_major_formatter(FuncFormatter(x_label_f))

    # Show the plot
    plt.tight_layout()
    if save_filename is not None:
        plt.savefig(save_filename, dpi=300, bbox_inches="tight")
    if show:
        plt.show()
    else:
        plt.close()


def find_x_scale(duration_sim: Quantity["time"]) -> tuple[str, Callable]:
    """Determine the appropriate x-axis scale and label based on the simulation duration.

    Args:
        duration_sim (Quantity["time"]): The total duration of the simulation, represented as a time quantity.

    Returns:
        tuple: A tuple containing the x-axis label (str) and the corresponding function to convert time values (function).

    The function updates the x-axis scale and label based on the total simulation time:
    - For durations less than or equal to 8 hours, the x-axis is labeled in minutes.
    - For durations less than or equal to 3 days, the x-axis is labeled in hours.
    - For durations greater than 3 days, the x-axis is labeled in days.
    """
    if duration_sim <= 8 * u.h:
        x_label_f = seconds_to_minutes
        x_label = r"Time $\left[min\right]$"
    elif duration_sim <= 3 * u.day:
        x_label_f = seconds_to_hours
        x_label = r"Time $\left[hour\right]$"
    else:
        x_label_f = seconds_to_days
        x_label = r"Time $\left[day\right]$"
    return (x_label, x_label_f)


def plot_boolean_bars(
    bool_array: np.ndarray,
    tofs: np.ndarray,
    duration_sim: Quantity["time"],
    title,
    save_filename: str = None,
    show: bool = True,
) -> None:
    """
    Plot horizontal bars when the boolean variable is True (or 1),
    and nothing when False (or 0).
    """
    # Convert the boolean array to integers (0, 1)
    bool_array = [int(b) for b in bool_array]

    # Update the x-axis scale
    x_label, x_label_f = find_x_scale(duration_sim)

    plt.figure(figsize=(6, 2))

    # Manually add horizontal grid lines at each mode level to make grid
    plt.axhline(y=0, color="gray", linestyle="--", linewidth=0.8)

    # Create the step plot using horizontal lines (hlines)
    for i in range(1, len(tofs)):
        if bool_array[i - 1] == 1:
            plt.hlines(
                0, tofs[i - 1], tofs[i], colors="blue", linewidth=8
            )  # Plot a bar when True (1)

    # Extend the last line if the last value is True (1)
    length_unit = duration_sim.value / len(bool_array)
    if bool_array[-1] == 1:
        plt.hlines(0, tofs[-1], tofs[-1] + length_unit, colors="blue", linewidth=8)

    plt.title(title, fontsize=13, fontweight="medium")
    plt.xlabel(x_label, fontsize=11)

    # Hide y-ticks and y-axis label, focus on horizontal bars
    plt.yticks([])
    plt.ylabel("")

    # Adjust the axes
    ax = plt.gca()
    ax.spines["top"].set_color((0.8, 0.8, 0.8))
    ax.spines["right"].set_color((0.8, 0.8, 0.8))
    ax.spines["left"].set_color((0.8, 0.8, 0.8))
    ax.spines["bottom"].set_color((0.8, 0.8, 0.8))
    ax.xaxis.set_major_formatter(FuncFormatter(x_label_f))

    plt.tight_layout()

    # Save the plot if a filename is provided
    if save_filename is not None:
        plt.savefig(save_filename, dpi=300, bbox_inches="tight")

    # Show or close the plot
    if show:
        plt.show()
    else:
        plt.close()


def plot_dashboard(
    modes: np.ndarray,
    eclipses: np.ndarray,
    visibility: np.ndarray,
    tofs: np.ndarray,
    duration_sim: Quantity["time"],
    title: str,
    save_filename: str = None,
    show: bool = True,
) -> None:
    """Plot operating modes during simulation to visualize them throughout time with the ecplipse status and visibility windows."""

    mode_remap = {
        0: 1,
        1: 2,
        2: 0,
        3: 4,
        4: 5,
        5: 3,
    }  # to change plot order on the y axis, so that charging mode is next to eclipse bar and xband next to the visibility bar
    new_mode_order = [2, 0, 1, 5, 3, 4]

    # Apply the remapping
    modes = [mode_remap[int(mode)] for mode in modes]

    # Update what the x_axis scale should be
    x_label, x_label_f = find_x_scale(duration_sim)

    # Mode names
    mode_labels = [mode_dict[key] for key in sorted(mode_dict.keys())]
    mode_labels = [mode_dict[key] for key in new_mode_order]
    colors = ["#5cb85c", "#d9534f", "#5bc0de", "#f0ad4e", "#a0522d", "#6f42c1"]
    colors = [colors[i] for i in new_mode_order]

    fig, axes = plt.subplots(
        3,
        1,
        figsize=(6, 3.2),
        gridspec_kw={"height_ratios": [1, 6, 1]},
        sharex=True,
    )

    # AX 0
    for i in range(1, len(tofs)):
        if visibility[i - 1] == 1:
            axes[0].hlines(0, tofs[i - 1], tofs[i], colors="pink", linewidth=5)
        else:
            axes[0].hlines(0, tofs[i - 1], tofs[i], colors="grey", linewidth=5)

    # Extend the last line if the last value is True (1)
    length_unit = duration_sim.value / len(modes)
    if visibility[-1] == 1:
        axes[0].hlines(0, tofs[-1], tofs[-1] + length_unit, colors="pink", linewidth=5)

    # Hide y-ticks and y-axis label, focus on horizontal bars
    axes[0].set_yticks([])
    axes[0].set_ylabel("Visibility", fontsize=8)
    axes[0].yaxis.label.set(rotation="horizontal", ha="right", va="center")

    axes[0].xaxis.set_ticks_position("none")  # Remove x-ticks
    axes[0].spines["top"].set_color((0.8, 0.8, 0.8))
    axes[0].spines["right"].set_color((0.8, 0.8, 0.8))
    axes[0].spines["left"].set_color((0.8, 0.8, 0.8))
    axes[0].spines["bottom"].set_color((0.8, 0.8, 0.8))
    axes[0].xaxis.set_major_formatter(FuncFormatter(x_label_f))
    axes[0].set_title(title, fontsize=11, fontweight="medium")

    visibility_handle = plt.Line2D(
        [0], [0], color="orchid", linewidth=5, label="Satellite is visible"
    )

    # AX 1
    axes[1].set_ylim(-0.6, len(mode_labels) - 0.3)  # Adjust limits to avoid space

    # Manually add horizontal grid lines at each mode level to make grid
    for y in range(len(mode_labels)):
        axes[1].axhline(
            y=y, color="gray", linestyle="--", linewidth=0.8
        )  # Add horizontal grid lines

    # Create the step plot using horizontal lines (hlines)
    for i in range(1, len(tofs)):
        axes[1].hlines(
            modes[i - 1],
            tofs[i - 1],
            tofs[i],
            colors=colors[modes[i - 1]],
            linewidth=10,
        )
    axes[1].hlines(
        modes[-1],
        tofs[-1],
        tofs[-1] + length_unit,
        colors=colors[modes[-1]],
        linewidth=10,
    )  # Extend last line

    # Add a legend
    mode_handles = [
        plt.Line2D([0], [0], color=colors[i], linewidth=5)
        for i in range(len(mode_labels))
    ]

    # Hide y-ticks and y-axis label
    axes[1].set_yticks([])
    axes[1].set_ylabel("Operating\nMode", fontsize=8)
    axes[1].yaxis.label.set(rotation="horizontal", ha="right", va="center")
    axes[1].xaxis.set_ticks_position("none")  # Remove x-ticks
    axes[1].spines["top"].set_color((0.8, 0.8, 0.8))
    axes[1].spines["right"].set_color((0.8, 0.8, 0.8))
    axes[1].spines["left"].set_color((0.8, 0.8, 0.8))
    axes[1].spines["bottom"].set_color((0.8, 0.8, 0.8))
    axes[1].xaxis.set_major_formatter(FuncFormatter(x_label_f))

    # AX 2
    for i in range(1, len(tofs)):
        if eclipses[i - 1] == 1:
            axes[2].hlines(0, tofs[i - 1], tofs[i], colors="grey", linewidth=5)
        else:
            axes[2].hlines(0, tofs[i - 1], tofs[i], colors="gold", linewidth=5)

    # Extend the last line if the last value is True (1)
    length_unit = duration_sim.value / len(modes)
    if eclipses[-1] == 1:
        axes[2].hlines(0, tofs[-1], tofs[-1] + length_unit, colors="gold", linewidth=5)

    # Hide y-ticks and y-axis label, focus on horizontal bars
    axes[2].set_yticks([])
    axes[2].set_ylabel("Sunlight", fontsize=8)
    axes[2].yaxis.label.set(rotation="horizontal", ha="right", va="center")

    axes[2].xaxis.set_ticks_position("none")  # Remove x-ticks
    axes[2].spines["top"].set_color((0.8, 0.8, 0.8))
    axes[2].spines["right"].set_color((0.8, 0.8, 0.8))
    axes[2].spines["left"].set_color((0.8, 0.8, 0.8))
    axes[2].spines["bottom"].set_color((0.8, 0.8, 0.8))
    axes[2].xaxis.set_major_formatter(FuncFormatter(x_label_f))
    axes[2].tick_params(axis="x", labelsize=8)
    axes[2].set_xlabel(x_label, fontsize=9)

    plt.rcParams["legend.title_fontsize"] = 9
    axes[2].legend(
        mode_handles,
        mode_labels,
        loc="upper center",
        bbox_to_anchor=(0.5, -2.5),
        ncol=3,
        title="Modes",
        fontsize=8,
    )

    # Show the plot
    fig.tight_layout()
    plt.subplots_adjust(
        left=None, bottom=None, right=None, top=None, wspace=None, hspace=0
    )
    if save_filename is not None:
        plt.savefig(save_filename, dpi=300, bbox_inches="tight")
    if show:
        plt.show()
    else:
        plt.close()


def plot_orbital_elem_evolution(
    tofs: np.ndarray,
    RAANs: np.ndarray,
    AOPs: np.ndarray,
    ECCs: np.ndarray,
    INCs: np.ndarray,
    altitudes: np.ndarray,
    folder: str,
    duration_sim: Quantity["time"],
):
    """Plot the evolution of RAAN, AOP, ECC, INC, and altitude during the simulation."""
    ticks_angle = np.array([0, 1 / 2 * np.pi, np.pi, 3 / 2 * np.pi, 2 * np.pi])
    tick_labels_angle = np.array([r"$0$", r"$\pi/2$", r"$\pi$", r"$3\pi/2$", r"$2\pi$"])

    x_label, x_label_f = find_x_scale(duration_sim)

    # update how many values to step to obtain clear plotting (approximately 100 points)
    step = int(len(tofs) / 100)

    plot_1d(
        tofs.to_value("second"),
        RAANs,
        "RAAN Over Time",
        x_label,
        r"RAAN ($rad$)",
        step=step,
        fill_under=False,
        remove_box=True,
        y_range=(min(RAANs) - 0.2, max(RAANs) + 0.2),
        scatter=True,
        custom_y_ticks=True,
        y_ticks=ticks_angle,
        y_tick_labels=tick_labels_angle,
        x_label_f=x_label_f,
        show=False,
        save_filename=folder + "RAAN.png",
    )

    plot_1d(
        tofs.to_value("second"),
        AOPs,
        "Argument of Periapsis Over Time",
        x_label,
        r"Argp ($rad$)",
        step=step,
        fill_under=False,
        remove_box=True,
        y_range=(min(AOPs) - 0.2, max(AOPs) + 0.2),
        custom_y_ticks=True,
        y_ticks=ticks_angle,
        y_tick_labels=tick_labels_angle,
        x_label_f=x_label_f,
        show=False,
        save_filename=folder + "AOP.png",
    )

    plot_1d(
        tofs.to_value("second"),
        ECCs,
        "Eccentricity Over Time",
        x_label,
        r"Ecc",
        step=step,
        fill_under=False,
        remove_box=True,
        y_range=(min(ECCs) - 0.1, max(ECCs) + 0.1),
        scatter=True,
        x_label_f=x_label_f,
        show=False,
        save_filename=folder + "ECC.png",
    )

    plot_1d(
        tofs.to_value("second"),
        INCs,
        "Inclination over time",
        x_label,
        r"Inc ($rad$)",
        step=step,
        fill_under=False,
        remove_box=True,
        y_range=(min(INCs) - 0.2, max(INCs) + 0.2),
        custom_y_ticks=True,
        y_ticks=ticks_angle,
        y_tick_labels=tick_labels_angle,
        x_label_f=x_label_f,
        show=False,
        save_filename=folder + "INC.png",
    )

    plot_1d(
        tofs.to_value("second"),
        altitudes,
        "Spacecraft Altitude Over Time",
        x_label,
        r"Altitude ($km$)",
        step=step,
        fill_under=False,
        remove_box=True,
        y_range=(min(altitudes) - 20, max(altitudes) + 20),
        x_label_f=x_label_f,
        show=False,
        save_filename=folder + "altitude.png",
    )
