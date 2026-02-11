"""
Module: Analytics Utilities

Author: Souren Ishkhanian

This module is mainly used for development and analysis of the underlying dataset
and is not directly related to Trojan analysis.
"""


import time
import numpy as np
from Visualizer import Visualizer
from matplotlib import pyplot as plt
from collections import defaultdict
from scipy.stats import rayleigh
from Alignment_Utilities import align_cells
from Annotation_Helpers import *


def check_efficiency(function: Callable[..., Any], args: tuple[Any, ...]) -> Any:
    """Measure and print the execution time of a function.

    Executes `function` with the provided `args` tuple, prints the elapsed
    time in seconds, and returns the function's result.

    Args:
        function (Callable[..., Any]): The function to be executed and timed.
        args (Tuple[Any, ...]): Positional arguments to pass to the function.

    Returns:
        result (Any): The return value of the executed function.
    """
    start = time.perf_counter()
    ret = function(*args)
    print(f"Time: {time.perf_counter() - start :.4f} seconds.")
    return ret


def alignment_test(
        cell_type: str = "BLS",
        sorted_cells: dict[str, list[Cell]] = {},
        representatives: dict[str, list[Point]] = {}
    ) -> None:
    """Test and visualize alignment for a specific cell type.

    Computes alignment distances for the given cell type, prints basic
    statistics, displays the first 1000 aligned cells, and plots the
    distribution of distances.

    Args:
        cell_type (str, optional): Name of the cell type to test. Defaults to "BLS".
        sorted_cells (dict[str, list[Cell]], optional): Mapping from cell type
            to all available cells. Defaults to an empty dict.
        representatives (dict[str, list[Point]], optional): Representative
            points for each cell type. Defaults to an empty dict.
    """
    aligned, dists = check_efficiency(align_cells, (sorted_cells[cell_type], representatives[cell_type]))
    dists = np.array(dists)
    print("min:", dists.min())
    print("max:", dists.max())
    print("mean:", dists.mean())
    print("median:", np.median(dists))
    print(f"Top 10: {np.sort(dists)[-10:]}")
    v = Visualizer(aligned[:1000], representatives[cell_type])
    v.display_all()
    plt.plot(dists)
    plt.xlabel("Distance")
    plt.ylabel("Count")
    plt.title("Distribution of distances")
    plt.show()


def show_boxsize_distribution(
        boxes: dict[str, tuple[int | float, int | float]]
    ) -> None:
    """Visualize the distribution of box widths and heights.

    Generates two bar charts showing how many cells fall into each width
    and height value.

    Args:
        boxes (dict[str, tuple[int | float, int | float]]): Mapping from cell
            type to its (width, height).
    """
    widths: dict[int | float, int] = defaultdict(int)
    heights: dict[int | float, int] = defaultdict(int)

    for (width, height) in boxes.values():
        widths[width] += 1
        heights[height] += 1

    plt.figure()
    plt.bar(widths.keys(), widths.values()) # type: ignore
    plt.title("Width distribution")
    plt.xlabel("Width")
    plt.ylabel("Count")
    plt.show()

    plt.figure()
    plt.bar(heights.keys(), heights.values()) # type: ignore
    plt.title("Height distribution")
    plt.xlabel("Height")
    plt.ylabel("Count")
    plt.show()


def plot_distance_distributions(
        dists_list: list[list[Optional[float]]],
        threshold: int = 2
    ) -> None:
    """Plot Rayleigh distributions of multiple distance datasets.

    For each list of distances in `dists_list`, fits a Rayleigh distribution
    and plots its probability density function. Distances with fewer than
    `threshold` valid entries are skipped.

    Args:
        dists_list (Sequence[Sequence[float | None]]): Lists of distances to plot.
        threshold (int, optional): Minimum number of valid distances required
            to include a dataset. Defaults to 2.
    """
    for dists in dists_list:
        d = np.array([x for x in dists if x is not None])
        if len(d) < threshold:
            continue
        params = rayleigh.fit(d)
        x = np.linspace(0, d.max(), 200)
        y = rayleigh.pdf(x, *params)
        plt.plot(x, y)

    plt.xlabel("Distance")
    plt.ylabel("Density")
    plt.title("Rayleigh Distributions")
    plt.show()


def distance_distribution_own_label(
        aligned_cells: dict[str, tuple[list[Optional[float]], list[Optional[float]]]],
        cell_types: list[str],
        threshold: int
    ) -> None:
    """Plot Rayleigh distance distributions for specified cell types.

    For each cell type in `cell_types`, filters out `None` distances, fits a
    Rayleigh distribution, and plots its probability density function if the
    number of valid distances exceeds `threshold`.

    Args:
        aligned_cells (dict[str, tuple[list[float | None], list[float | None]]]):
            Mapping from cell type to a tuple of (cells, distances).
        cell_types (Sequence[str]): Cell types to include in the plot.
        threshold (int): Minimum number of valid distances required to plot.
    """
    for cell_type in cell_types:
        _, dists = aligned_cells[cell_type]
        valid_dists = [dist for dist in dists if dist is not None]
        if len(valid_dists) > threshold:
            x = np.linspace(0, max(valid_dists), 100)
            params = rayleigh.fit(valid_dists)
            pdf = rayleigh.pdf(x, *params)
            plt.plot(x, pdf, label=cell_type)
    plt.xlabel("Distance")
    plt.ylabel("Density")
    plt.legend()
    plt.show()