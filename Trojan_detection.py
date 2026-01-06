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
from Labeling_Utilities import predict_trojans
from Alignment_Utilities import align_cells
from Distance_Measures import *

DATASETS = {
"28nm_chip" : {"path" : Path("./Data/Chip_Data_28nm.pickle")}
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


def group_by_width(sorted_cells=[], boxes=[]):
    groups = defaultdict(list)

    for cell_type, (width, _) in boxes.items():
        key = width
        for w in groups:
            if abs(w - width) <= 10:
                key = w
                break
        groups[key].append(cell_type)
    
    grouped = {}
    for group in groups.values():
        widths = [boxes[cell_type][0] for cell_type in group]
        grouped[(min(widths), max(widths))] = group

    rows = []
    for w, types in grouped.items():
        total_vias = 0
        for t in types:
            cells = sorted_cells[t]
            majority = Counter(len(c["vias"]) for c in cells).most_common(1)[0][0]
            total_vias += len(cells) * majority

        rows.append({
            "Width": w,
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
    width_df = group_by_width(sorted_cells, boxes)
    width_groups = {row['Width']: row['Types'] for _, row in width_df.iterrows()}

    trojans = cache.get_confirmed_trojans()
    for width, group in width_groups.items():
        group_representatives = {k : v for (k, v) in representatives.items() if k in group}
        cells = []
        dists = []
        for cell_type in group:
            c, d = aligned_cells[cell_type]
            cells.extend(c)
            dists.extend(d)
        trojans[width] = predict_trojans(cells, dists, group_representatives)
    cache.save_confirmed_trojans(trojans)


if __name__ == "__main__":
    dataset = "28nm_chip"
    chip_data_file = DATASETS[dataset]["path"]
    cache = DataCache(chip_data_file)

    sorted_cells = cache.get_sorted_cells()
    cell_types = list(sorted_cells.keys())
    boxes = cache.get_boxes()
    aligned_cells = cache.get_aligned_cells()
    representatives = cache.get_representatives()

    print(
f"""Cell types: {len(sorted_cells)}
Total Boxes: {len(boxes)}
Total Representatives: {len(representatives)}
Total aligned cells: {len(aligned_cells)}""")
    
    width_df = group_by_width(sorted_cells, boxes)
    width_df = width_df[width_df["Width"] == (268, 268)]
    print(width_df["Types"].tolist()[0])
    print(get_cell_type_info(data_overview(sorted_cells), width_df["Types"].tolist()[0]))

    all_trojans = cache.get_confirmed_trojans()

    # For each width range show how many cells were predicted as trojans
    for w, trojans in all_trojans.items():
        print(f"{str(w):<15}: {len(trojans)}")

    # For each width range start inspecting trojans if there are any
    # After inspection confirmed trojans are updated in the cache
    for w, trojans in all_trojans.items():
        if not trojans: continue
        confirmed: list[Predicted_Cell] = []
        inspector = Cell_Inspector(trojans, confirmed, representatives, chamfer)
        inspector.start_interactive()
        all_trojans[w] = confirmed
    # cache.save_confirmed_trojans(all_trojans)
