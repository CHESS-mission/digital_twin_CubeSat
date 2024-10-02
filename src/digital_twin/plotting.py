from typing import List

from matplotlib import pyplot as plt
from matplotlib.ticker import FuncFormatter
import numpy as np


# Function for general plotting in 2d
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
    x_tick_label: np.array = None,
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
