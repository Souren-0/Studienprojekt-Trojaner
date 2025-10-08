# from collections import defaultdict, Counter
# from matplotlib import pyplot as plt
# from Visualizer import Visualizer
# from Cell_Via_Utilities import *
# from pprint import pprint
# import Distance_Measures
# import statistics
from Data_Manager import *


DATASETS = {
"28nm_chip" : {"path" : Path("./Data/Chip_Data_28nm.pickle")}
}


if __name__ == "__main__":
    dataset = "28nm_chip"
    chip_data_file = DATASETS[dataset]["path"]

    cache = DataCache(chip_data_file)

    sorted_cells = cache.get_sorted_cells()
    cell_types = list(sorted_cells.keys())
    boxes = cache.get_boxes()
    representatives = cache.get_representatives()

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
    #     total_vias.append((cell_type, total_cells * vias_per_cell, total_cells, vias_per_cell))
    # # total_vias.sort(key = lambda x: x[3])
    # pprint(total_vias)

    # cell_type = "BGY"
    # cells = [reset_transform(cell) for cell in sorted_cells[cell_type]]
    # box = boxes[cell_type]

    # via_counts = Counter()
    # for cell in cells:
    #     via_counts[len(cell["vias"])] += 1

    # labels, counts = zip(*via_counts)

    # plt.bar(labels, counts)
    # plt.xlabel("Number of Vias")
    # plt.ylabel("Count of Cells")
    # plt.title("Vias per Cell")
    # plt.show()

    # start = time.perf_counter()
    # representative = find_representative_vias(cells, plot=True, multiprocess=True)
    # print(f"Finding a representative took {time.perf_counter() - start:.4f} seconds.")
