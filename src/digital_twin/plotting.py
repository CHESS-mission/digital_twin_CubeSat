from typing import List

from matplotlib import pyplot as plt
from matplotlib.ticker import FuncFormatter
import numpy as np
import plotly.graph_objects as go

from astropy.time import Time, TimeDelta
from astropy import units as u
from astropy.coordinates import (
    GCRS,
    ITRS,
    CartesianDifferential,
    CartesianRepresentation,
    SphericalRepresentation,
)
from poliastro.twobody import Orbit
from poliastro.plotting import StaticOrbitPlotter, OrbitPlotter3D
from poliastro.earth.plotting import GroundtrackPlotter


# Function for general plotting in 2d with x and y arrays as input
def plot_1d(
    x: np.array,
    y: np.array,
    title: str,
    xlabel: str,
    ylabel: str,
    step: int = 1,
    fill_under: bool = True,
    remove_box: bool = True,
    x_range: tuple = None,
    y_range: tuple = None,
    scatter: bool = False,
    x_label_f: callable = None,
    custom_y_ticks: bool = False,
    y_ticks: np.array = None,
    y_tick_labels: np.array = None,
    custom_x_ticks: bool = False,
    x_ticks: np.array = None,
    x_tick_labels: np.array = None,
    save_filename: str = None,
    show: bool = True,
):
    # Downsample the x and y arrays
    x_downsampled = x[::step]
    y_downsampled = y[::step]

    # Create the plot
    plt.figure(figsize=(6, 4))
    if scatter:
        plt.scatter(x_downsampled, y_downsampled, marker=".", color="blue")
    else:
        plt.plot(x_downsampled, y_downsampled, marker=".", linestyle="-", color="blue")

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

    # Show the plot
    plt.tight_layout()
    if save_filename is not None:
        plt.savefig(save_filename)
    if show:
        plt.show()
    else:
        plt.close()


# Function to convert seconds to days to change label of x axis
def seconds_to_days(x, pos):
    """Convert seconds to days."""
    return f"{x / 86400:.0f}"


# Function to convert seconds to hours to change label of x axis
def seconds_to_hours(x, pos):
    """Convert seconds to hours."""
    return f"{x / 3600:.0f}"


# Function to convert seconds to minutes to change label of x axis
def seconds_to_minutes(x, pos):
    """Convert seconds to minutes."""
    return f"{x / 60:.0f}"


def plot_orbit_trajectory_3d(
    title: str,
    attractor_name: str,
    folder: str,
    orbit: Orbit = None,
    label_orbit: str = "",
    label_traj: str = "",
    traj: np.array = None,
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
):
    """Plot an orbit in 3d around its attractor.

    Option 1: with a grid and tick values but no tick lines.
    Option 2: with no grid not tick values

    Args:
        orbit (Orbit): The orbit from the class Orbit to plot
        title (str): Title of the figure
        option (int, optional): Which kind of plotting desired, as explained above. Defaults to 1.
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
    attractor_name: str,
    folder: str,
    label: str,
    orbit: Orbit,
    orbit_color: str = "black",
    orbit_width: float = 2,
    orbit_type: str = "solid",
    attractor_color: str = "mediumblue",
    sat_color: str = "red",
    show: bool = False,
):

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
    rr: np.array,
    epochs: np.array,
    label: str,
    folder: str,
    traj_color: str = "purple",
    traj_width: float = 1,
    station_coords: np.array = None,
    station_name: str = "station",
    station_color: str = "red",
    show: bool = False,
):
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
    if station_coords is not None:
        station = station_coords * u.deg  # [LAT LON]
        gp.add_trace(
            go.Scattergeo(
                lat=station[0],
                lon=station[-1],
                name=station_name,
                marker={"color": station_color},
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
