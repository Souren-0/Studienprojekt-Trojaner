# from collections import defaultdict, Counter
import numpy as np
import pandas as pd
from pprint import pprint
from Data_Manager import *
from collections import Counter
from scipy.stats import rayleigh
from matplotlib import pyplot as plt
from Visualizer import Visualizer
from Cell_Inspector import Cell_Inspector
from Labeling_Utilities import predict_trojans, prune_predicted_trojans, confirm_predicted_trojans
from Alignment_Utilities import align_cells
from Distance_Measures import *

DATASETS = {
"28nm_chip" : {"path" : Path("./Data/Chip_Data_28nm.pickle")},
"65nm_chip" : {"path" : Path("./Data/Chip_Data_65nm.pickle")}
}


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


def data_overview(sorted_cells):
    rows = []
    for cell_type, cells in sorted_cells.items():
        cell_vias = [len(cell["vias"]) for cell in cells]
        majority_vias = Counter(cell_vias).most_common(1)[0][0]
        rows.append({
            "type": cell_type,
            "amount": len(cells),
            "vias": majority_vias,
            "total_vias" : len(cells) * majority_vias
        })
    df = pd.DataFrame(rows)
    df = df.sort_values("total_vias", ascending=False)
    return df


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


def group_boxsizes(boxes, w_gap=5, h_gap=5):
    widths, heights = map(sorted, map(list, map(set, zip(*boxes.values()))))

    width_ranges = []
    w_start = widths[0]
    for i in range(len(widths) - 1):
        if widths[i+1] - widths[i] > w_gap:
            width_ranges.append((w_start, widths[i]))
            w_start = widths[i+1]
    width_ranges.append((w_start, widths[-1]))

    height_ranges = []
    h_start = heights[0]
    for i in range(len(heights) - 1):
        if heights[i+1] - heights[i] > h_gap:
            height_ranges.append((h_start, heights[i]))
            h_start = heights[i+1]
    height_ranges.append((h_start, heights[-1]))

    return width_ranges, height_ranges


def group_by_boxsize(boxes, w_gap=5, h_gap=5):
    widths, heights = group_boxsizes(boxes, w_gap, h_gap)
    grouped = defaultdict(list)
    for cell_type, (width, height) in boxes.items():
        key = (
            next(w_min for w_min, w_max in widths if w_min <= width <= w_max),
            next(h_min for h_min, h_max in heights if h_min <= height <= h_max)
        )
        grouped[key].append(cell_type)
    return grouped


def group_by_width(sorted_cells, boxes):
    grouped = group_by_boxsize(boxes)

    rows = []
    for (w, h), types in grouped.items():
        total_vias = 0
        for t in types:
            cells = sorted_cells[t]
            majority = Counter(len(c["vias"]) for c in cells).most_common(1)[0][0]
            total_vias += len(cells) * majority

        rows.append({
            "Width": w,
            "Height": h,
            "Types_amount": len(types),
            "Types": types,
            "Total_vias": total_vias
        })

    df = pd.DataFrame(rows).sort_values("Total_vias", ascending=False)
    return df


def get_cell_type_info(df, cell_types):
    return df[df['type'].isin(cell_types)]


def get_mapping(dataset="28nm"):
    with open("./Data/Cell_Mapping.pickle", "rb") as f:
        data = pickle.load(f)
    return data[dataset]


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


def start_trojan_scan(sorted_cells, aligned_cells, representatives, boxes, cache: DataCache):
    groups = group_by_boxsize(boxes, w_gap=15, h_gap=15)
    
    distances = {}
    for cell_type, (cells, dists) in aligned_cells.items():
        distances[cell_type] = [d for d in dists if d is not None]

    print("Starting trojan scan...")
    start = time.perf_counter()
    trojans = cache.get_trojans()
    for i, (width, group) in enumerate(groups.items()):
        group_representatives = {k : v for (k, v) in representatives.items() if k in group}
        cells = []
        dists = []
        for cell_type in group:
            c, d = aligned_cells[cell_type]
            cells.extend(c)
            dists.extend(d)
        print(f"Scanning boxsize {width} ({i}/{len(groups)})\nThere are {len(group)} different types in this group.")
        possible_trojans = predict_trojans(cells, dists, group_representatives)
        possible_trojans = prune_predicted_trojans(possible_trojans)
        possible_trojans = confirm_predicted_trojans(possible_trojans, distances, representatives)
        trojans[width] = possible_trojans
        print(f"Found {len(possible_trojans)} potential trojans.")
        cache.save_trojans(trojans)
        print("Trojans cached.")
    print("Trojan scan completed. Retrieve the results from the cache")
    print(f"Scanning took {time.perf_counter() - start :.4f} seconds.")


def read_trojans(path):
    with open(path, "rb") as f:
        trojans = pickle.load(f)
    return trojans


if __name__ == "__main__":
    dataset = "65nm_chip"
    chip_data_file = DATASETS[dataset]["path"]
    cache = DataCache(chip_data_file)

    # cache.update_sorted_cells()
    # cache.group_cells(get_mapping("65nm"))

    sorted_cells = cache.get_sorted_cells()
    cell_types = list(sorted_cells.keys())

    # cache.update_boxes()
    # cache.update_representatives(cell_types, reset=True)
    # cache.update_aligned_cells(cell_types, reset=True)

    boxes = cache.get_boxes()
    aligned_cells = cache.get_aligned_cells()
    representatives = cache.get_representatives()

    print(
f"""Cell types: {len(sorted_cells)}
Total Boxes: {len(boxes)}
Total Representatives: {len(representatives)}
Total aligned cells: {len(aligned_cells)}""")

    # start_trojan_scan(sorted_cells, aligned_cells, representatives, boxes, cache)
    all_trojans = cache.get_trojans()

    # all_trojans = read_trojans("trojans_65nm.pickle")

    # # For each width range show how many cells were predicted as trojans
    # total = 0
    # for w, trojans in all_trojans.items():
    #     n = len(trojans)
    #     if not n: continue
    #     print(f"{str(w):<15}: {len(trojans)}")
    #     for trojan in trojans: print(f"    {trojan['actual']} -> {trojan['predicted']}")
    #     total += len(trojans)
    # print(f"In total: {total} predicted trojans")

    # For each width range start inspecting trojans if there are any
    # After inspection confirmed trojans are updated in the cache
    for w, trojans in all_trojans.items():
        if not trojans: continue
        confirmed: list[Predicted_Cell] = []
        inspector = Cell_Inspector(trojans, confirmed, representatives, chamfer)
        inspector.start_interactive()
        all_trojans[w] = confirmed
    # cache.save_trojans(all_trojans)

# Notes:
# Original: 372
# Pruning: 85
# Confirming: 223 (directed)
# Confirming: 200 (undirected)
# Pruning + Confirming: 36