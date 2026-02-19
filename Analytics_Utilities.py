"""
Module: Analytics Utilities

Author: Souren Ishkhanian

This module is mainly used for development and analysis of the underlying dataset
and is not directly related to Trojan analysis.
"""


import time
import numpy as np
import pandas as pd
from scipy import stats
from Visualizer import Visualizer
from matplotlib import pyplot as plt
from collections import Counter, defaultdict
from Alignment_Utilities import align_cells
from Annotation_Helpers import *
pd.set_option('display.max_rows', None)


def data_overview(sorted_cells: dict[str, list[Cell]]) -> pd.DataFrame:
    """Return an overview of sorted cells as a table with the following columns:
    - type: Cell type name
    - amount: Number of cells of this type
    - vias: Most common via count per cell for this type
    - total_vias: amount times vias

    Args:
        sorted_cells (dict[str, list[Cell]]):
            Mapping from cell type to all available cells of that type.

    Returns:
        overview (pd.DataFrame): Summary table with one row per cell type.
    """
    rows = []
    for cell_type, cells in sorted_cells.items():
        cell_vias = [len(cell["vias"]) for cell in cells]
        majority_vias = Counter(cell_vias).most_common(1)[0][0] if cells else 0
        rows.append({
            "type": cell_type,
            "amount": len(cells),
            "vias": majority_vias,
            "total_vias" : len(cells) * majority_vias
        })
    df = pd.DataFrame(rows)
    df = df.sort_values("total_vias", ascending=False)
    return df


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
        params = stats.rayleigh.fit(d)
        x = np.linspace(0, d.max(), 200)
        y = stats.rayleigh.pdf(x, *params)
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
            params = stats.rayleigh.fit(valid_dists)
            pdf = stats.rayleigh.pdf(x, *params)
            plt.plot(x, pdf, label=cell_type)
    plt.xlabel("Distance")
    plt.ylabel("Density")
    plt.legend()
    plt.show()


def fit_comparison(sorted_cells, aligned_cells, row, dists_list, percents=[85,90,95]):
    """
    Plot and compare fitted probability distributions against empirical distance data.

    For the cell type selected by `row` in `sorted_cells`, this function extracts
    alignment distances, fits multiple SciPy distributions, overlays their PDFs,
    and marks selected upper-percentile cutoffs.

    Args:
        sorted_cells:
            Output used by `data_overview`; selects the cell type via `row`.
        aligned_cells:
            Mapping from cell type to (cells, distances), where distances may contain None.
        row:
            Row index into `data_overview(sorted_cells)` identifying the cell type.
        dists_list:
            Iterable of (name, scipy.stats distribution) pairs
            (e.g. [("F", stats.f), ("Log-Laplace", stats.loglaplace)]).
        percents:
            Percentiles (0–100) for which vertical cutoff lines are drawn.
    """
    print("Plotting distributions of distances for")
    print(data_overview(sorted_cells).iloc[row])
    cell_type = data_overview(sorted_cells).iloc[row]['type']
    _, dists = aligned_cells[cell_type]
    x = np.array([d for d in dists if d is not None])
    xs = np.linspace(x.min(), x.max(), 10000)

    plt.hist(x, bins=500, density=True, alpha=0.35, label="Data")

    for dist_name, dist_func in dists_list:
        params = dist_func.fit(x)
        line, = plt.plot(xs, dist_func.pdf(xs, *params), label=f"{dist_name} fit")
        for p in percents:
            q = dist_func.ppf(p/100, *params)
            ls = "-" if p==max(percents) else "--" if p==90 else ":"
            plt.axvline(q, color=line.get_color(), lw=1, ls=ls, label=f"{dist_name} {p}%")

    plt.legend()
    plt.show()