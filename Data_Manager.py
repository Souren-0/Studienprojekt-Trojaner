"""
Module: Data Manager
Author: Souren Ishkhanian

As the chip dataset is large, and redunant computations - such as sorting cells or computing representatives - is costly, this module computes
and caches data for each dataset, that remain consistent throughout several runs.
"""

from Cell_Via_Utilities import find_representative_vias
from collections import Counter, defaultdict
from Annotation_Helpers import *
from multiprocessing import Pool
from pathlib import Path
import warnings
import pickle
import time

class DataCache:
    def __init__(self, data_path: Path) -> None:
        if (not data_path.exists()):
            raise FileNotFoundError(f"Data path does not exist: {data_path}")

        self.data_path = data_path
        self.cache_dir = Path("./data_cache") / data_path.stem
        self.cache_dir.mkdir(parents=True, exist_ok=True)

        start = time.perf_counter()
        self.sorted_cells = self._retrieve_sorted_cells()
        self.representatives = self._retrieve_representatives()
        self.boxes = self._retrieve_boxes()
        print(f"Retrieving took {time.perf_counter() - start:.4f} seconds")

    # -------- Public API --------
    def update_all(self) -> None:
        start = time.perf_counter()
        self.update_sorted_cells()
        self.update_representatives(list(self.get_sorted_cells().keys()), cell_num=100)
        self.update_boxes()
        warnings.warn("Used cell_num=100 for quicker update. Consider treating these representatives as a draft. \n" \
            "Update representatives with higher cell_num for better accuracies.", UserWarning)
        print(f"Updating took {time.perf_counter() - start:.4f} seconds")

    def update_sorted_cells(self) -> None:
        print("Updating sorted cells...")
        with open(self.data_path, "rb") as f:
            chip_data = pickle.load(f)

        tile_list = list(chip_data.keys())
        all_cells = []
        for tile in tile_list:
            all_cells += chip_data[tile]

        sorted_cells: dict[str, list[Cell]] = defaultdict(list)
        for cell in all_cells:
            sorted_cells[cell["data"]["name"]].append(cell)
        
        self._cache_sorted_cells(sorted_cells)
        print("Updating done.")

    def update_representatives(self, cell_types: list[str], cell_num: int = 1000, replace: bool = False, reset: bool = False) -> None:
        print("Updating representatives...")
        representatives = {} if reset else self.get_representatives()
        cell_types = cell_types if replace else [cell_type for cell_type in cell_types if cell_type not in representatives]

        sorted_cells = self.get_sorted_cells()
        if cell_types:
            with Pool() as pool:
                new_reps = pool.starmap(find_representative_vias, [(sorted_cells[cell_type], cell_num) for cell_type in cell_types])
            new_representatives = dict(zip(cell_types, new_reps))
        else: new_representatives = {}
        
        representatives.update(new_representatives)
        self._cache_representatives(representatives)
        print("Updating done.")

    def update_boxes(self) -> None:
        print("Updating boxes...")
        sorted_cells = self.get_sorted_cells()
        cell_boxes: dict[str, Counter[tuple[float, float]]] = defaultdict(Counter)
        cell_types = list(sorted_cells)
        for cell_type in cell_types:
            for cell in sorted_cells[cell_type]:
                top_left, bottom_right = cell['box']
                box_size = (bottom_right[0] - top_left[0], bottom_right[1] - top_left[1])
                cell_boxes[cell_type][box_size] += 1

        cell_box = {}
        for cell_type, boxes in cell_boxes.items():
            widths, heights = zip(*boxes)
            cell_box[cell_type] = (max(widths), max(heights))

        self._cache_boxes(cell_box)
        print("Updating done.")


    def get_sorted_cells(self) -> dict[str, list[Cell]]:
        return self.sorted_cells
    
    def get_representatives(self) -> dict[str, list[tuple[float, float]]]:
        return self.representatives
    
    def get_boxes(self) -> dict[str, tuple[float, float]]:
        return self.boxes


    def group_cells(self, cell_mapping: dict[str, str]) -> None:
        grouped_cells: dict[str, list[Cell]] = defaultdict(list)
        sorted_cells = self.get_sorted_cells()
        for cell_type, cells in sorted_cells.items():
            for cell in cells:
                cell_map = cell_mapping.get(cell_type, cell_type)
                cell["data"]["name"] = cell_map
                grouped_cells[cell_map].append(cell)
        self._cache_sorted_cells(grouped_cells)


    # -------- Private functions --------
    def _retrieve_sorted_cells(self) -> dict[str, list[Cell]]:
        print("Retrieving sorted cells.")
        path = self.cache_dir / "sorted_cells.pickle"
        if path.exists():
            with open(path, "rb") as f:
                return pickle.load(f)
            
        warnings.warn("No data is cached. Call 'update_sorted_cells()' to cache.", UserWarning)
        return {}

    def _retrieve_representatives(self) -> dict[str, list[tuple[float, float]]]:
        print("Retrieving representatives.")
        path = self.cache_dir / "representatives.pickle"
        if path.exists():
            with open(path, "rb") as f:
                return pickle.load(f)
            
        warnings.warn("No data is cached. Call 'update_representatives()' to cache.", UserWarning)
        return {}

    def _retrieve_boxes(self) -> dict[str, tuple[float, float]]:
        print("Retrieving boxes.")
        path = self.cache_dir / "boxes.pickle"
        if path.exists():
            with open(path, "rb") as f:
                return pickle.load(f)
            
        warnings.warn("No data is cached. Call 'update_boxes()' to cache.", UserWarning)
        return {}


    def _cache_sorted_cells(self, new_cells: dict[str, list[Cell]]) -> None:
        print("Saving new cells")
        self.sorted_cells = new_cells
        path = self.cache_dir / "sorted_cells.pickle"
        with open(path, "wb") as f:
            pickle.dump(new_cells, f)

    def _cache_representatives(self, new_reps: dict[str, list[tuple[float, float]]]) -> None:
        print("Saving new representatives")
        self.representatives = new_reps
        path = self.cache_dir / "representatives.pickle"
        with open(path, "wb") as f:
            pickle.dump(new_reps, f)

    def _cache_boxes(self, new_boxes: dict[str, tuple[float, float]]) -> None:
        print("Saving new boxes")
        self.boxes = new_boxes
        path = self.cache_dir / "boxes.pickle"
        with open(path, "wb") as f:
            pickle.dump(new_boxes, f)