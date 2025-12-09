# from collections import defaultdict, Counter
import pandas as pd
from matplotlib import pyplot as plt
# from Visualizer import Visualizer
from pprint import pprint
# import Distance_Measures
import statistics
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


def alignment_test(cell_type="BLS"):
    aligned, dists = check_efficiency(align_cells, (sorted_cells[cell_type], representatives[cell_type]))
    dists = np.array(dists)
    print("min:", dists.min())
    print("max:", dists.max())
    print("mean:", dists.mean())
    print("median:", np.median(dists))
    print(f"Top 10: {np.sort(dists)[-10:]}")
    v = Visualizer(boxes[cell_type], aligned[:1000], representatives[cell_type])
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


def group_by_width():
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
    width_df = group_by_width()
    cell_type_group = width_df.iloc[2]["Types"][1:]
    
    cell_type = "BKI"
    cells_dists = aligned_cells[cell_type]
    # trojans = check_cells_for_trojan(*aligned_cells[cell_type])
    # print(len(trojans))
    reps = {rep : representatives[rep] for rep in cell_type_group}
    cells_to_check = check_cells_for_trojan(*cells_dists, 0.8)[:10]
    assigned = defaultdict(list)
    for cell in tqdm(cells_to_check):
        label = assign_cell_type_biased(cell, reps, 0.5)
        assigned[label].append(cell)
    
    # pprint(assigned[list(assigned.keys())[0]])
    counts = {rep: len(cells) for rep, cells in assigned.items()}
    plt.bar(counts.keys(), counts.values())
    plt.xticks(rotation=45)
    plt.ylabel("Count")
    plt.show()

    v = Visualizer(assigned['BJI'], representatives['BJI'], representatives[cell_type])
    v.display_all()

    # cell_type = "BKK"
    # trojans = check_cells_for_trojan(*aligned_cells[cell_type])
    # print(len(trojans))
    # v = Visualizer(boxes[cell_type], trojans, representatives[cell_type])
    # v.display_all()

    # print(width_df[:10])
    # print(get_cell_type_info(overview_df, cell_type_group))
    # cache.update_aligned_cells(cell_type_group, replace=True)

    # dist_matrix = {}
    # for rep1 in cell_type_group:
    #     rep_vias = representatives[rep1]
    #     other_reps = [{"vias": representatives[rep2]} for rep2 in cell_type_group]
    #     _, dists = align_cells(other_reps, representatives[rep1], boxes[rep1]) # type: ignore
    #     dist_matrix[rep1] = dists
    # pprint(dist_matrix)

    # reps = list(cell_type_group)
    # # convert dict of lists into a square matrix
    # M = np.array([dist_matrix[r] for r in reps], dtype=float)
    # plt.imshow(M, cmap="viridis", interpolation="nearest")
    # plt.colorbar(label="Distance")
    # plt.xticks(range(len(reps)), reps, rotation=90)
    # plt.yticks(range(len(reps)), reps)
    # plt.title("Representative Distance Matrix")
    # plt.tight_layout()
    # plt.show()

    # flatten
    # all_dists = []
    # for d in dist_matrix.values():
    #     all_dists.extend([x for x in d if x is not None and x >= 0])

    # all_dists = np.array(all_dists)
    # if len(all_dists) == 0:
    #     print("No valid distances")
    # else:
    #     params = rayleigh.fit(all_dists)
    #     x = np.linspace(0, all_dists.max(), 300)
    #     pdf = rayleigh.pdf(x, *params)

    #     plt.plot(x, pdf)
    #     plt.show()

    # for cell_type in cell_type_group[3:5]:
    #     print(cell_type)
    #     v = Visualizer(boxes[cell_type], cache.get_sorted_cells()[cell_type][:1000], representatives[cell_type])
    #     v.display_all()
    #     v = Visualizer(boxes[cell_type], cache.get_aligned_cells()[cell_type][0][:1000], representatives[cell_type])
    #     v.display_all()

    # distance_distribution_own_label(aligned_cells, cell_type_group, 2)

    # cell_type = "BKE"
    # temp_BKI = ["BKI", "BKK", "BJI", "AMK", "BJS", "BAS", "ALE"]
    # temp_BKE = ["BKE", "BKI", "BKK", "BHE", "BJS", "BAS", "ALE"]
    # dist_dict = {rep: align_cells(sorted_cells[cell_type], representatives[rep], boxes[rep])[1] for rep in temp_BKE}

    # for rep, dists in dist_dict.items():
    #     d = [x for x in dists if x is not None and x >= 0]
    #     if len(d) < 2: 
    #         continue

    #     params = rayleigh.fit(d)
    #     x = np.linspace(0, max(d), 200)
    #     pdf = rayleigh.pdf(x, *params)

    #     plt.plot(x, pdf, label=rep)

    # plt.legend()
    # plt.xlabel("Distance")
    # plt.ylabel("Density")
    # plt.show()

    # cell_type = cell_type_group[1]
    # print(cell_type)
    # aligned, _ = aligned_cells[cell_type]
    # aligned = sorted_cells[cell_type]
    # aligned = aligned[:1000]
    # v = Visualizer(boxes[cell_type], aligned, representatives[cell_type])
    # v.display_all()

    # grouped_representatives: list[dict[str, list[Point]]] = []
    # for rep_group in grouped_width_celltypes.values():
    #     grouped_representatives.append({rep : representatives[rep] for rep in rep_group})
    # grouped_representatives.sort(key=lambda x: len(x), reverse=True)
    
    # group = grouped_representatives[10]
    # group_cells = []
    # for rep in group:
    #     group_cells.extend(sorted_cells[rep])
    #     print(rep, ":", len(sorted_cells[rep]))
    # print(list(group.keys()))
    # print(len(group_cells))

    # valid = 0
    # possible_trojan = []
    # for cell in tqdm(group_cells[:10]):
    #     label = assign_cell_type(cell, group)
    #     cell_type = cell["data"]["name"]
    #     print(label, cell_type)
    #     if label == cell_type: valid += 1
    #     else: possible_trojan.append(cell)
    # print(valid)
    # print(len(possible_trojan))

    # with open("Data/Cell_Mapping.pickle", "rb") as f:
    #     data = pickle.load(f)

    # grouping = data["28nm"]
    # for group in grouping:
    #     if "BIW" in group: print(group)
    #     # if "ACQ" in group: print(group)
    # print(len(boxes))

    # mapping = data["28nm"]
    # cache.group_cells(mapping)
    # print(len(cache.get_sorted_cells()))

    cell_type = "BEU" # Same as ARE
    # v = Visualizer(boxes[cell_type], sorted_cells[cell_type], representatives[cell_type])
    # v.display_all()
    # find_representative_vias(sorted_cells[cell_type], plot=True, multiprocess=True)
    
    cell_type = "ARE" # Same as BEU
    # v = Visualizer(boxes[cell_type], sorted_cells[cell_type], representatives[cell_type])
    # v.display_all()
    # find_representative_vias(sorted_cells[cell_type], plot=True, multiprocess=True)
    
    # trojan = possible_trojan[0]
    # v = Visualizer(trojan["box"], possible_trojan, representatives[label], representatives[trojan["data"]["name"]])
    # v.display_all()

    cell_type = "BIW" # type with the most cell instances (BIW, BDA, GU)
    # cells = sorted_cells[cell_type]
    # possible_trojans = check_cells_for_trojan(cells, representatives[cell_type], 0.9999999)
    # print(len(possible_trojans))
    
    # v = Visualizer(boxes[cell_type], possible_trojans[:3000], representatives[cell_type])
    # v.display_all()

    # sorted_widths = list(grouped_width_celltypes.keys())
    # sorted_widths.sort()
    # diffs = []
    # for i in range(len(sorted_widths[:-1])):
    #     diffs.append(sorted_widths[i+1] - sorted_widths[i])
    # # pprint(diffs)
    # width_type_count = []
    # for _, types in grouped_width_celltypes.items():
    #     width_type_count.append(len(types))
    # width_type_count = list(filter(lambda x: x > 1, width_type_count))
    # width_type_count.sort()
    # pprint(width_type_count)
    # print(len(width_type_count))

    # Not so important anymore
    # Via-cell count ratios
    # total_vias = []
    # for cell_type, cells in sorted_cells.items():
    #     total_cells = len(cells)
    #     vias_per_cell = int(statistics.median([len(cell["vias"]) for cell in cells]))
    #     if vias_per_cell > 0:
    #         total_vias.append((cell_type, total_cells * vias_per_cell, total_cells, vias_per_cell))
    # total_vias.sort(key = lambda x: x[2])
    # pprint(total_vias)

    # cell_type = "BIW"
    # cells = [reset_transform(cell) for cell in sorted_cells[cell_type]]
    # # box = boxes[cell_type]

    # via_counts = Counter()
    # for cell in cells:
    #     via_counts[len(cell["vias"])] += 1
        
    # labels, counts = zip(*via_counts.items())

    # plt.bar(labels, counts)
    # plt.xlabel("Number of Vias")
    # plt.ylabel("Count of Cells")
    # plt.title("Vias per Cell")
    # plt.show()

    # start = time.perf_counter()
    # representative = find_representative_vias(cells, plot=True, multiprocess=True)
    # print(f"Finding a representative took {time.perf_counter() - start:.4f} seconds.")

    # cache.update_aligned_cells(["BLG", "BKK"])
    # for cell_type in aligned_cells:
    #     trojans = sort_cells(*aligned_cells[cell_type])
    #     v = Visualizer(boxes[cell_type], trojans, representatives[cell_type])
    #     v.display_all()

    # cell_type = "BLS"
    # cache.update_aligned_cells([cell_type])
    # cache.update_aligned_cells([cell_type], replace=True)
    # trojans = check_cells_for_trojan(*aligned_cells[cell_type])
    # print(len(sorted_cells[cell_type]))
    # print(len(trojans))
    # v = Visualizer(boxes[cell_type], trojans, representatives[cell_type])
    # v.display_all()
    # print(len(sorted_cells[cell_type]))
