# from collections import defaultdict, Counter
# from multiprocessing import Process, Queue
from matplotlib import pyplot as plt
from Cell_Via_Utilities import *
# import Distance_Measures
from Data_Manager import *
from Visualizer import Visualizer
from pprint import pprint


DATASETS = {
"28nm_chip" : {"path" : Path("./Data/Chip_Data_28nm.pickle")}
}

if __name__ == "__main__":
    dataset = "28nm_chip"
    chip_data_file = DATASETS[dataset]["path"]

    cache = DataCache(chip_data_file)

    sorted_cells = cache.get_sorted_cells()
    cell_types = list(sorted_cells.keys())

    # for cell_type in cell_types:
    #     cells = sorted_cells[cell_type]
    #     print(f"{cell_type}: {len(cells):03}")

    # cells = sorted_cells['BLA']
    # cells = [reset_transform(cell) for cell in cells]
    # vias = get_aligned_vias(cells, num_cells=10000)[1]
    # vias = [via for cell in cells for via in cell['vias']]
    # plt.scatter(*zip(*vias), s=5)
    # plt.show()
    # visualizer = Visualizer(cells[0]["box"], cells[:500])
    # visualizer.display_all()


    # filtered_cells = {}
    # for cell_type, cells in sorted_cells.items():
    #     n = len(cells)
    #     if n > 2 and cells[0]["vias"]:
    #         if n <= 5:
    #             print(f"{cell_type} contains only {n} cells. Consider ignoring it.")
    #         filtered_cells[cell_type] = cells

    # all_cells = [cell for cells in filtered_cells.values() for cell in cells]
    # counts = Counter([len(cell['vias']) for cell in all_cells])
            
    # pprint(counts)
    # print(f"Filtering reduced the cell types from {len(sorted_cells)} to {len(filtered_cells)}.")
    # # del counts[0]
    # max_count = max(counts)
    # for cell in all_cells:
    #     if len(cell["vias"]) == max_count:
    #         name = cell["data"]["name"]
    #         print(f"Cell of type {name} contains {max_count} vias")
    #         print(f"There are in total of {len(sorted_cells[name])} cells of that type.")
    #         v = Visualizer(cell["box"], [cell])
    #         v.display_all()
    
    # vias_counts, freq = zip(*counts.items())

    # plt.bar(vias_counts, freq)
    # plt.xlabel("Number of Vias")
    # plt.ylabel("Frequency")
    # plt.title("Distribution of Vias per Cell")
    # plt.show()



    boxes = cache.get_boxes()
    # widths, heights = zip(*boxes.values())
    
    # width_counts = Counter(widths)
    # height_counts = Counter(heights)

    # # Plot widths
    # plt.bar(width_counts.keys(), width_counts.values())
    # plt.title("Width Counts")
    # plt.xlabel("Width")
    # plt.ylabel("Count")
    # plt.show()

    # # Plot heights
    # plt.bar(height_counts.keys(), height_counts.values())
    # plt.title("Height Counts")
    # plt.xlabel("Height")
    # plt.ylabel("Count")
    # plt.show()


    # Grouping types by their widths +/- 3 margin
    grouped_width_celltypes: dict[int | float, list[str]] = defaultdict(list)
    for cell_type, (width, height) in boxes.items():
        key = width
        for w in grouped_width_celltypes.keys():
            if abs(w - width) <= 3:
                key = w
                break
        grouped_width_celltypes[key].append(cell_type)

    sorted_widths = list(grouped_width_celltypes.keys())
    sorted_widths.sort()
    diffs = []
    for i in range(len(sorted_widths[:-1])):
        diffs.append(sorted_widths[i+1] - sorted_widths[i])
    # pprint(diffs)
    width_type_count = []
    for _, types in grouped_width_celltypes.items():
        width_type_count.append(len(types))
    width_type_count = list(filter(lambda x: x > 1, width_type_count))
    width_type_count.sort()
    pprint(width_type_count)
    print(len(width_type_count))

