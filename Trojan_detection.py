# from collections import defaultdict, Counter
import pandas as pd
from matplotlib import pyplot as plt
# from Visualizer import Visualizer
from pprint import pprint
# import Distance_Measures
import statistics
from Cell_Inspector import Cell_Inspector
from Data_Manager import *
from Cell_Via_Utilities import *
from collections import Counter
from scipy.stats import rayleigh

# Benutze Sphinx für Dokumentation

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
    grouped = defaultdict(list)

    for cell_type, (width, _) in boxes.items():
        key = width
        for w in grouped:
            if abs(w - width) <= 3:
                key = w
                break
        grouped[key].append(cell_type)

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


def get_cell_type_info(df, cells):
    return df[df['type'].isin(cells)]


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

    overview_df = data_overview(sorted_cells)
    # print(overview_df[:20])

    width_df = group_by_width(sorted_cells, boxes)
    cell_type_group = width_df.iloc[2]["Types"][1:]
    cell_type_group_representatives = {k : v for (k, v) in representatives.items() if k in cell_type_group}
    cell_type_group_overview = overview_df[overview_df["type"].isin(cell_type_group)]
    # cell_type_group_overview["aligned"] = (cell_type_group_overview["type"].isin(aligned_cells))
    # print(cell_type_group_overview)

    predicted_trojans = []

    for cell_type in cell_type_group:
        cells, dists = aligned_cells[cell_type]
        predicted_trojans.extend(predict_trojans(cells, dists, cell_type_group_representatives))
    approved_trojans = []
    inspector = Cell_Inspector(predicted_trojans, approved_trojans, representatives, chamfer)
    inspector.start_interactive()
    
    print(len(approved_trojans))
    print(approved_trojans[0]["cell"])
