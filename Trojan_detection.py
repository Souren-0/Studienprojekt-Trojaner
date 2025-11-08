# from collections import defaultdict, Counter
from matplotlib import pyplot as plt
# from Visualizer import Visualizer
from pprint import pprint
# import Distance_Measures
import statistics
from Data_Manager import *
from Cell_Via_Utilities import *


# Benutze Sphinx für Dokumentation
# Rigit registration für alignment

DATASETS = {
"28nm_chip" : {"path" : Path("./Data/Chip_Data_28nm.pickle")}
}


if __name__ == "__main__":
    dataset = "28nm_chip"
    chip_data_file = DATASETS[dataset]["path"]

    cache = DataCache(chip_data_file)

    sorted_cells = cache.get_sorted_cells()
    aligned_cells = cache.get_aligned_cells()
    cell_types = list(sorted_cells.keys())
    boxes = cache.get_boxes()
    representatives = cache.get_representatives()

    print(len(sorted_cells))
    print(len(boxes), len(representatives))
    # print("BIW:", len(sorted_cells["BIW"]))
    # print("BDA:", len(sorted_cells["BDA"]))
    # print("GU:", len(sorted_cells["GU"]))

    # v = Visualizer(boxes["BIW"], sorted_cells["BDA"])
    # v.display_all()

    # !!! Important for next step (hopefully tomorrow) !!!
    # Grouping types by their widths +/- 3 margin
    # grouped_width_celltypes: dict[int | float, list[str]] = defaultdict(list)
    # for cell_type, (width, height) in boxes.items():
    #     key = width
    #     for w in grouped_width_celltypes.keys():
    #         if abs(w - width) <= 3:
    #             key = w
    #             break
    #     grouped_width_celltypes[key].append(cell_type)
    
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
    for cell_type in aligned_cells:
        trojans = sort_cells(*aligned_cells[cell_type])
        v = Visualizer(boxes[cell_type], trojans, representatives[cell_type])
        v.display_all()
