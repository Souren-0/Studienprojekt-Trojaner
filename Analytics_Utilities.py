import time
import numpy as np
from Visualizer import Visualizer
from matplotlib import pyplot as plt
from collections import defaultdict
from scipy.stats import rayleigh
from Alignment_Utilities import align_cells


def check_efficiency(function, args):
    start = time.perf_counter()
    ret = function(*args)
    print(f"Time: {time.perf_counter() - start :.4f} seconds.")
    return ret


def alignment_test(cell_type="BLS", sorted_cells=[], representatives=[]):
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


def show_boxsize_distribution(boxes):
    widths = defaultdict(int)
    heights = defaultdict(int)

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


def plot_distance_distributions(dists_list, threshold=2):
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


def distance_distribution_own_label(aligned_cells, cell_types, threshold):
    for cell_type in cell_types:
        _, dists = aligned_cells[cell_type]
        dists = [dist for dist in dists if dist is not None]
        if len(dists) > threshold:
            x = np.linspace(0, max(dists), 100)
            params = rayleigh.fit(dists)
            pdf = rayleigh.pdf(x, *params)
            plt.plot(x, pdf, label=cell_type)
    plt.xlabel("Distance")
    plt.ylabel("Density")
    plt.legend()
    plt.show()