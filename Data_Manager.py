"""
Module: Data Manager
Author: Souren Ishkhanian

As the chip dataset is large, and redunant computations - such as sorting cells or computing representatives - is costly, this module computes
and caches data for each dataset, that remain consistent throughout several runs.
"""

from collections import Counter, defaultdict
from Annotation_Helpers import *
from Cell_Via_Utilities import *
from multiprocessing import Pool
from pathlib import Path
from copy import deepcopy
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
        self.aligned_cells = self._retrieve_aligned_cells()
        print(f"Retrieving took {time.perf_counter() - start:.4f} seconds")

    # -------- Public API --------
    def update_all(self) -> None:
        start = time.perf_counter()
        self.update_sorted_cells()
        self.update_boxes()
        self.update_representatives(list(self.get_sorted_cells().keys()), cell_num=100)
        self.update_aligned_cells(list(self.get_sorted_cells().keys()))
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

    def update_representatives(self, cell_types: list[str], cell_num: Optional[int] = 1000, replace: bool = False, reset: bool = False) -> None:
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

    def update_aligned_cells(self, cell_types: list[str], cell_num: Optional[int] = None, replace: bool = False, reset: bool = False) -> None:
        print("Updating aligned cells...")
        aligned_cells = {} if reset else self.get_aligned_cells()
        cell_types = cell_types if replace else [cell_type for cell_type in cell_types if cell_type not in aligned_cells]

        new_cells = {cell_type : 
                         align_all_cells(
                             self.get_sorted_cells().get(cell_type, [])[:cell_num],
                             self.get_representatives().get(cell_type, None)
                            )
                        for cell_type in cell_types}
        
        aligned_cells.update(new_cells)
        self._cache_aligned_cells(aligned_cells)
        print("Updating done")


    def group_cells(self, cell_mapping: dict[str, str]) -> None:
        grouped_cells: dict[str, list[Cell]] = defaultdict(list)
        sorted_cells = self.get_sorted_cells()
        for cell_type, cells in sorted_cells.items():
            for cell in cells:
                cell_map = cell_mapping.get(cell_type, cell_type)
                cell["data"]["name"] = cell_map
                grouped_cells[cell_map].append(cell)
        self._cache_sorted_cells(grouped_cells)


    def get_sorted_cells(self) -> dict[str, list[Cell]]: return self.sorted_cells
    def get_representatives(self) -> dict[str, list[tuple[float, float]]]: return self.representatives
    def get_boxes(self) -> dict[str, tuple[float, float]]: return self.boxes
    def get_aligned_cells(self) -> dict[str, tuple[list[Cell], list[float]]]: return self.aligned_cells


    # -------- Private functions --------
    def _retrieve_cache(self, name: str) -> dict:
        print(f"Retrieving {name.replace('_', ' ')}")
        path = self.cache_dir / f"{name}.pickle"
        if path.exists():
            with open(path, "rb") as f:
                return pickle.load(f)
        warnings.warn(f"No data is cached. Call 'update_{name}()' to cache.", UserWarning)
        return {}

    def _retrieve_sorted_cells(self): return self._retrieve_cache("sorted_cells")
    def _retrieve_representatives(self): return self._retrieve_cache("representatives")
    def _retrieve_boxes(self): return self._retrieve_cache("boxes")
    def _retrieve_aligned_cells(self): return self._retrieve_cache("aligned_cells")


    def _cache_data(self, name: str, data: dict) -> None:
        print(f"Saving new {name.replace('_', ' ')}")
        setattr(self, name, data)
        path = self.cache_dir / f"{name}.pickle"
        with open(path, "wb") as f:
            pickle.dump(data, f)

    def _cache_sorted_cells(self, new_cells): self._cache_data("sorted_cells", new_cells)
    def _cache_representatives(self, new_reps): self._cache_data("representatives", new_reps)
    def _cache_boxes(self, new_boxes): self._cache_data("boxes", new_boxes)
    def _cache_aligned_cells(self, new_cells): self._cache_data("aligned_cells", new_cells)