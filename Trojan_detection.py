from Data_Manager import *
from pathlib import Path

DATASETS = {
"28nm_chip" : {"path" : Path("./Data/Chip_Data_28nm.pickle")}
}

if __name__ == "__main__":
    dataset = "28nm_chip"
    chip_data_file = DATASETS[dataset]["path"]

    cache = DataCache(chip_data_file)

    sorted_cells = cache.get_sorted_cells()
    cell_types = list(sorted_cells.keys())
