# from collections import defaultdict, Counter
import numpy as np
import pandas as pd
from pprint import pprint
from Data_Manager import *
from collections import Counter
from Visualizer import Visualizer
from Cell_Inspector import Cell_Inspector
from Labeling_Utilities import predict_trojans, prune_predicted_trojans, confirm_predicted_trojans
from Alignment_Utilities import sample
from Distance_Measures import *


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


def group_by_boxsize(boxes, w_gap=5, h_gap=5):
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
    
    widths, heights = group_boxsizes(boxes, w_gap, h_gap)
    grouped = defaultdict(list)
    for cell_type, (width, height) in boxes.items():
        key = (
            next(w_min for w_min, w_max in widths if w_min <= width <= w_max),
            next(h_min for h_min, h_max in heights if h_min <= height <= h_max)
        )
        grouped[key].append(cell_type)
    return grouped


def show_alignment_example(sorted_cells, aligned_cells, representatives, sample_num=1000, example_num=0):
    if not sorted_cells:
        print("No example can be shown. Did you forget to fill the cache?")
        return
    
    cell_type = data_overview(sorted_cells).iloc[
        min(example_num, len(sorted_cells) - 1)
    ]['type']

    if cell_type not in representatives:
        print(f"Representative for the cell type \"{cell_type}\" is not present.")
        representative = None
    else: representative = representatives[cell_type]

    if cell_type not in aligned_cells:
        print(f"The cell type \"{cell_type}\" was not aligned yet. Considering unaligned cells.")
        cells = sorted_cells[cell_type]
    else: cells = aligned_cells[cell_type][0]

    print(f"There are {len(cells)} cells of type \"{cell_type}\". Showing {sample_num} random cells...")
    v = Visualizer(sample(cells, sample_num), representative)
    v.display_all()


def show_trojan_predictions(trojans, cache=None):
    for w, trojans in all_trojans.items():
        if not trojans: continue
        confirmed = []
        inspector = Cell_Inspector(trojans, confirmed, representatives, chamfer)
        inspector.start_interactive()
        all_trojans[w] = confirmed
    if cache: cache.save_trojans(all_trojans)


def trojan_prediction_summary(trojans):
    print("Trojan summary" if trojans else "No trojans present. Did you run a trojan scan?")
    total = 0
    for w, trojan_set in trojans.items():
        n = len(trojan_set)
        if not n: continue

        print(f"{str(w):<15}: {n}")
        for trojan in trojan_set: print(f"    {trojan['actual']} -> {trojan['predicted']}")
        total += n
    print(f"In total: {total} predicted trojans")


def read_trojans(path):
    with open(path, "rb") as f:
        trojans = pickle.load(f)
    return trojans


def get_mapping(dataset="28nm", path="./Data/Cell_Mapping.pickle"):
    if not Path(path).exists:
        print("Cell types could not be grouped.")
        return {}
    
    with open(path, "rb") as f:
        data = pickle.load(f)
    return data[dataset]


def bool_prompt(msg,
        accepting=["y", "ye", "yea", "yeah", "yes", "yep", "yup", "10-4", "roger", "affirmative", "1"],
        rejecting=["n", "no", "nop", "nope", "nah", "never", "negative", "cancel", "abort", "stop", "0"]):
    input_str = input(msg + " (y/n): ")
    valid_input = accepting + rejecting
    while input_str.lower() not in valid_input:
        input_str = input("Invalid input. Please try again: ")
    return input_str in accepting


def start_trojan_scan(aligned_cells, representatives, boxes, cache: DataCache, bias_strength=0.5, safe_cell_thr=0.9):
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
        possible_trojans = predict_trojans(cells, dists, group_representatives, bias_strength)
        possible_trojans = prune_predicted_trojans(possible_trojans)
        possible_trojans = confirm_predicted_trojans(possible_trojans, distances, representatives)
        trojans[width] = possible_trojans
        print(f"Found {len(possible_trojans)} potential trojans.")
        cache.save_trojans(trojans)
        print("Trojans cached.")
    print("Trojan scan completed. Retrieve the results from the cache")
    print(f"Scanning took {time.perf_counter() - start :.4f} seconds.")


def fill_cache(cache):
    cache.update_sorted_cells()
    cache.group_cells(get_mapping(cache.data_info["name"]))
    sorted_cells = cache.get_sorted_cells()
    cell_types = list(sorted_cells.keys())
    cache.update_boxes()
    cache.update_representatives(cell_types, reset=True)
    cache.update_aligned_cells(cell_types, reset=True)


DATASETS = {
    "28nm_chip" : {"name": "28nm", "path" : Path("./Data/type_bins_new_vias_28nm.pickle")},
    "40nm_chip" : {"name": "40nm", "path" : Path("./Data/type_bins_new_vias_40nm.pickle")},
    "65nm_chip" : {"name": "65nm", "path" : Path("./Data/type_bins_new_vias_65nm.pickle")},
    "90nm_chip" : {"name": "90nm", "path" : Path("./Data/type_bins_new_vias_90nm.pickle")}
}


FULL_PROGRAM = False
if __name__ == "__main__":
    if FULL_PROGRAM: start = time.perf_counter()

    dataset = DATASETS["28nm_chip"]
    cache = DataCache(dataset)

    if FULL_PROGRAM or bool_prompt("Do you want to fill the cache?"):
        fill_cache(cache)

    sorted_cells = cache.get_sorted_cells()
    cell_types = list(sorted_cells.keys())
    boxes = cache.get_boxes()
    aligned_cells = cache.get_aligned_cells()
    representatives = cache.get_representatives()

    print(f"Cache content summary:\n\
    Total Cell types: {len(sorted_cells)}\n\
    Total Boxes: {len(boxes)}\n\
    Total Representatives: {len(representatives)}\n\
    Total aligned cells: {len(aligned_cells)}")
    
    if not FULL_PROGRAM:
        show_alignment_example(sorted_cells, aligned_cells, representatives, example_num=0)
    if FULL_PROGRAM or bool_prompt("Run Trojan scan?"):
        start_trojan_scan(aligned_cells, representatives, boxes, cache)
    print("Trojan scan completed.")
    
    all_trojans = cache.get_trojans()
    # all_trojans = read_trojans("trojans_28nm_new.pickle")
    trojan_prediction_summary(all_trojans)

    if FULL_PROGRAM: print(f"Total time: {time.perf_counter() - start :.4f} seconds.") #type: ignore

    if FULL_PROGRAM or bool_prompt("Inspect predictions?"):
        show_trojan_predictions(all_trojans)
