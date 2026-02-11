"""
Module: Trojan Detection

Author: Souren Ishkhanian

This module provides the main pipeline for detecting Trojan cells in chip datasets.
It integrates cache management, cell alignment, representative selection, Trojan
prediction, pruning, and interactive confirmation. The module also includes utilities
for visualizing distributions, inspecting alignments, and summarizing predictions.
"""

import pandas as pd
from Data_Manager import *
from Distance_Measures import *
from collections import Counter
from Visualizer import Visualizer
from Cell_Inspector import Cell_Inspector
from collections import defaultdict, Counter
from Alignment_Utilities import sample
from Labeling_Utilities import predict_trojans, prune_predicted_trojans, confirm_predicted_trojans


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


def group_by_boxsize(
        boxes: dict[str, tuple[int | float, int | float]],
        w_gap: int = 5, h_gap: int = 5
    ) -> dict[tuple[int | float, int | float], list[str]]:
    """Group cell types by similar box dimensions.

    Cell boxes are grouped into width and height ranges. A new range is started whenever
    the difference between consecutive widths or heights exceeds the given gap. Each cell type
    is then assigned to the corresponding (`width_range_start`, `height_range_start`) group.

    Args:
        boxes (dict[str, tuple[int | float, int | float]]):
            Mapping from cell type to its (width, height).
        w_gap (int, optional):
            Maximum allowed width difference within a group. Defaults to 5.
        h_gap (int, optional):
            Maximum allowed height difference within a group. Defaults to 5.

    Returns:
        grouped (dict[tuple[int | float, int | float], list[str]]):
            Mapping from (`width_range_start`, `height_range_start`) to the list of
            cell types that fall into that size group.
    """
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


def show_alignment_example(
        sorted_cells: dict[str, list[Cell]],
        aligned_cells: dict[str, tuple[list[Cell], list[float | None]]],
        representatives: dict[str, list[tuple[int | float, int | float]]],
        sample_num: int = 1000,
        example_num: int = 0
    ) -> None:
    """Display a visual alignment example for a selected cell type.

    The cell type is chosen based on a summary ordering of `sorted_cells`.
    If aligned cells are available, those are shown; otherwise, unaligned
    cells are used. A random subset of cells is visualized.

    Args:
        sorted_cells (dict[str, list[Cell]]):
            Mapping from cell type to all available cells of that type.
        aligned_cells (dict[str, tuple[list[Cell], list[float | None]]]):
            Mapping from cell type to aligned cells and their alignment scores.
        representatives (dict[str, list[tuple[int | float, int | float]]]):
            Representative vias for each cell type.
        sample_num (int, optional):
            Number of cells to visualize. Defaults to 1000.
        example_num (int, optional):
            Index of the cell type to visualize in the overview ordering.
            Defaults to 0.
    """
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


def show_trojan_predictions(
        all_trojans: dict[tuple[int | float, int | float], list[Predicted_Cell]],
        representatives: dict[str, list[tuple[int | float, int | float]]],
        cache: Optional[DataCache] = None
    ) -> None:
    """Interactively inspect and confirm predicted Trojan cells.

    For each box-size group containing predictions, an interactive inspector
    is launched to confirm or reject detected Trojan candidates. The confirmed
    results replace the original predictions in-place.

    Args:
        all_trojans (dict[tuple[int | float, int | float], list[Predicted_Cell]]):
            Mapping from box size to predicted Trojan cells.
        representatives (dict[str, list[tuple[int | float, int | float]]]):
            Representative vias for each cell type.
        cache (DataCache, optional):
            Cache used to persist confirmed Trojan predictions.
            If provided, results are saved after inspection.
    """
    for w, trojans in all_trojans.items():
        if not trojans: continue
        confirmed: list[Predicted_Cell] = []
        inspector = Cell_Inspector(trojans, confirmed, representatives, chamfer)
        inspector.start_interactive()
        all_trojans[w] = confirmed
    if cache: cache.save_trojans(all_trojans)


def trojan_prediction_summary(
        all_trojans: dict[tuple[int | float, int | float], list[Predicted_Cell]]
    ) -> None:
    """Print a summary of predicted Trojan cells.

    For each box-size group, the number of predicted Trojans is shown along
    with their actual and predicted labels. A total count is printed at the end.

    Args:
        all_trojans (dict[tuple[int | float, int | float], list[Predicted_Cell]]):
            Mapping from box size to predicted Trojan cells.
    """
    print("Trojan summary" if all_trojans else "No trojans present. Did you run a trojan scan?")
    total = 0
    for w, trojan_set in all_trojans.items():
        n = len(trojan_set)
        if not n: continue

        print(f"{str(w):<15}: {n}")
        for trojan in trojan_set: print(f"    {trojan['actual']} -> {trojan['predicted']}")
        total += n
    print(f"In total: {total} predicted trojans")


def read_trojans(path: str) -> dict[tuple[int | float, int | float], list[Predicted_Cell]]:
    """Load predicted Trojan cells from a pickle file.
    If path does not exist an empty dictionary is returned.

    Args:
        path (str): Path to the pickle file containing Trojan predictions.

    Returns:
        trojans (dict[tuple[int | float, int | float], list[Predicted_Cell]]):
            Mapping from box size to predicted Trojan cells.
    """
    if not Path(path).exists():
        return {}
    
    with open(path, "rb") as f:
        trojans = pickle.load(f)
    return trojans


def make_mapping(
        src: str = "./Data/Cell_Groups.pickle",
        dest: str = "./Data/Cell_Mapping.pickle"
    ) -> dict[str, dict[str, str]]:
    """Create and save a mapping from each cell type to a representative cell type.

    Reads grouped cell data from `src`, then assigns each cell type to the first
    cell type in its group. Saves the resulting mapping to `dest`.

    Args:
        src (str, optional): Path to the source pickle file containing cell groups.
            Defaults to "./Data/Cell_Groups.pickle".
        dest (str, optional): Path where the resulting cell mapping will be saved.
            Defaults to "./Data/Cell_Mapping.pickle".

    Returns:
        mapping (dict[str, dict[str, str]]): For each dataset, a mapping from each cell type
        to its assigned representative cell type.
    """
    if not Path(src).exists():
        print("Cannot make a cell mapping. Invalid source Path.")
        return {}
    
    with open(src, "rb") as f:
        groups = pickle.load(f)
    
    type_mapping = {}
    for dataset, groupings in groups.items():
        mappings = {}
        for group in groupings:
            mapping = {cell_type: group[0] for cell_type in group} if group else {}
            mappings.update(mapping)
        type_mapping[dataset] = mappings

    with open(dest, "wb") as f:
        pickle.dump(type_mapping, f)
    
    return type_mapping


def get_mapping(
        dataset: str = "28nm",
        path: str = "./Data/Cell_Mapping.pickle"
    ) -> dict[str, str]:
    """Load the cell type mapping for a specific dataset.

    Reads the mapping pickle from `path` and returns the mapping for the
    requested dataset. If the file does not exist, an empty mapping is returned.

    Args:
        dataset (str, optional): Name of the dataset to retrieve. Defaults to "28nm".
        path (str, optional): Path to the cell mapping pickle. Defaults to
            "./Data/Cell_Mapping.pickle".

    Returns:
        mapping (dict[str, str]): Mapping from each cell type to its representative
            cell type for the specified dataset.
    """
    if not Path(path).exists():
        print("Cell types could not be grouped.")
        return {}
    
    with open(path, "rb") as f:
        data = pickle.load(f)
    return data[dataset]


def bool_prompt(
        msg: str,
        accepting: list[str] = ["y", "ye", "yea", "yeah", "yes", "yep", "yup", "10-4", "roger", "affirmative", "1"],
        rejecting: list[str] = ["n", "no", "nop", "nope", "nah", "never", "negative", "cancel", "abort", "stop", "0"]
    ) -> bool:
    """Prompt the user for a yes/no answer and return a boolean result.

    The function continues to prompt until the user enters a recognized response. Responses are case-insensitive.

    Args:
        msg (str): Message displayed to the user.
        accepting (list[str], optional): List of input strings interpreted as "yes".
            Defaults to ["y", "ye", "yea", "yeah", yes, ...].
        rejecting (list[str], optional): List of input strings interpreted as "no".
            Defaults to ["n", "no", "nop", "nope", "nah", ...].

    Returns:
        answer (bool): True if the user entered an accepting value, False if a rejecting value was entered.
    """
    input_str = input(msg + " (y/n): ")
    valid_input = accepting + rejecting
    while input_str.lower() not in valid_input:
        input_str = input("Invalid input. Please try again: ")
    return input_str.lower() in accepting


def start_trojan_scan(
        aligned_cells: dict[str, tuple[list[Cell], list[float | None]]],
        representatives: dict[str, list[tuple[int | float, int | float]]],
        boxes: dict[str, tuple[int | float, int | float]],
        cache: DataCache,
        bias_strength: float = 0.5,
        safe_cell_thr: float = 0.9):
    """Perform a Trojan cell scan across all box-size groups.

    Groups cells by box size, predicts potential Trojan cells using the
    provided aligned cells and representatives, prunes and confirms predictions,
    and saves results in the provided cache.

    Args:
        aligned_cells (dict[str, tuple[list[Cell], list[float | None]]]):
            Mapping from cell type to its aligned cells and alignment distances.
        representatives (dict[str, list[tuple[int | float, int | float]]]):
            Representative box dimensions for each cell type.
        boxes (dict[str, tuple[int | float, int | float]]):
            Box sizes for all cell types.
        cache (DataCache):
            Cache object used to store and retrieve Trojan predictions.
        bias_strength (float, optional):
            Strength of the bias in Trojan prediction. Defaults to 0.5.
        safe_cell_thr (float, optional):
            Threshold for treating a cell as safe. Defaults to 0.9.
    """
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
        possible_trojans = predict_trojans(cells, dists, group_representatives, bias_strength, safe_cell_thr)
        possible_trojans = prune_predicted_trojans(possible_trojans)
        possible_trojans = confirm_predicted_trojans(possible_trojans, distances, representatives)
        trojans[width] = possible_trojans
        print(f"Found {len(possible_trojans)} potential trojans.")
        cache.save_trojans(trojans)
        print("Trojans cached.")
    print("Trojan scan completed. Retrieve the results from the cache")
    print(f"Scanning took {time.perf_counter() - start :.4f} seconds.")


def fill_cache(cache: DataCache) -> None:
    """Populate a DataCache with all necessary cell information.

    Updates sorted cells, groups them according to the dataset mapping,
    computes box sizes, sets representatives, and aligns the cells.
    This prepares the cache for downstream analysis and visualization.

    Args:
        cache (DataCache): The cache object to populate.
    """
    cache.update_sorted_cells()
    cache.group_cells(get_mapping(cache.data_info["name"]))
    sorted_cells = cache.get_sorted_cells()
    cell_types = list(sorted_cells.keys())
    cache.update_boxes()
    cache.update_representatives(cell_types, reset=True)
    cache.update_aligned_cells(cell_types, reset=True)


# Set correct paths before running
DATASETS = {
    "28nm_chip" : {"name": "28nm", "path" : Path("./Data/type_bins_new_vias_28nm.pickle")},
    "40nm_chip" : {"name": "40nm", "path" : Path("./Data/type_bins_new_vias_40nm.pickle")},
    "65nm_chip" : {"name": "65nm", "path" : Path("./Data/type_bins_new_vias_65nm.pickle")},
    "90nm_chip" : {"name": "90nm", "path" : Path("./Data/type_bins_new_vias_90nm.pickle")}
}

# Change accordingly
CELL_GROUPS_PATH = "./Data/Cell_Groups.pickle"

# Set true if you want to skip intermediate prompts/checks
# Full program includes:
# - 
FULL_PROGRAM = True
"""
Set to `True` if you want to skip intermediate prompts/checks

Full program includes:
- Creating a cell mapping to group identical cell types (valid `CELL_GROUPS_PATH` required)
- Filling the cache (which includes computation of representatives and aligned cells)
- All steps of trojan scan (prediction -> pruning -> confirming)
"""

DATASET = DATASETS["40nm_chip"]

if __name__ == "__main__":
    if FULL_PROGRAM:
        start = time.perf_counter()
        make_mapping(CELL_GROUPS_PATH)

    cache = DataCache(DATASET)

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
        show_trojan_predictions(all_trojans, representatives)
