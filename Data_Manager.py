"""
Module: Data Manager <br>
Author: Souren Ishkhanian

As the chip dataset is large, and redunant computations - such as sorting cells or computing representatives - is costly, this module computes
and caches data for each dataset, that remain consistent throughout several runs.
"""

import time
import pickle
import warnings
from pathlib import Path
from Cell_Via_Utilities import reset_transform
from Alignment_Utilities import align_cells
from Labeling_Utilities import find_representative_vias
from Annotation_Helpers import *
from multiprocessing import Pool
from collections import Counter, defaultdict

class DataCache:
    """
    Manages cached computations for a chip dataset to avoid redundant processing.

    This class loads, computes, and caches frequently used data such as:
        - Sorted cells by type
        - Representative cell vias
        - Cell bounding boxes
        - Aligned cells

    Cached data is stored in a local directory (`./data_cache/<dataset_name>`) 
    as pickle files and is persistent across multiple runs.

    Attributes:
        data_path (Path): Path to the raw chip dataset.
        cache_dir (Path): Directory where cached data is stored.
        sorted_cells (dict[str, list[Cell]]): Cells sorted by their type/name.
        representatives (dict[str, list[tuple[float, float]]]): Representative via positions for each cell type.
        boxes (dict[str, tuple[float, float]]): Width and height of each cell type's bounding box.
        aligned_cells (dict[str, tuple[list[Cell], list[Optional[float]]]]): Aligned cells and optional alignment errors.
    """

    def __init__(self, data_path: Path) -> None:
        """
        Initialize a DataCache for the given dataset path.

        Args:
            data_path (Path): Path to the raw chip dataset pickle file.

        Raises:
            FileNotFoundError: If `data_path` does not exist.
        
        Behavior:
            - Creates a cache directory for the dataset if it doesn't exist.
            - Retrieves cached data (sorted_cells, representatives, boxes, aligned_cells)
              or initializes empty caches if missing.
            - Prints the time taken to retrieve all caches.
        """
        if (not data_path.exists()):
            raise FileNotFoundError(f"Data path does not exist: {data_path}")

        self.data_path = data_path
        self.cache_dir = Path("./data_cache") / data_path.stem

        if not self.cache_dir.exists():
            warnings.warn(
                "Cache directory created for the first time. "
                "You may want to call 'update_all()' to populate caches."
            )
        self.cache_dir.mkdir(parents=True, exist_ok=True)

        start = time.perf_counter()
        self.sorted_cells = self._retrieve_sorted_cells()
        self.representatives = self._retrieve_representatives()
        self.boxes = self._retrieve_boxes()
        self.aligned_cells = self._retrieve_aligned_cells()
        print(f"Retrieving took {time.perf_counter() - start:.4f} seconds")

    # -------- Public API --------
    def update_all(self) -> None:
        """
        Recompute and update all cached data: sorted cells, boxes, representatives, and aligned cells.
        Uses a reduced cell_num=100 for representatives for speed; treat as draft.
        """
        start = time.perf_counter()
        self.update_sorted_cells()
        self.update_boxes()
        self.update_representatives(list(self.get_sorted_cells().keys()))
        self.update_aligned_cells(list(self.get_sorted_cells().keys()))
        print(f"Total update took {time.perf_counter() - start:.4f} seconds")

    def update_sorted_cells(self) -> None:
        """
        Sort all cells in the dataset by type/name and cache the result.
        """
        print("Updating sorted cells...")
        start = time.perf_counter()
        with open(self.data_path, "rb") as f:
            chip_data = pickle.load(f)

        tile_list = list(chip_data.keys())
        all_cells = []
        for tile in tile_list:
            all_cells += chip_data[tile]

        sorted_cells: dict[str, list[Cell]] = defaultdict(list)
        for cell in all_cells:
            sorted_cells[cell["data"]["name"]].append(reset_transform(cell))
        
        self._cache_sorted_cells(sorted_cells)
        print(f"Updating done. Took {time.perf_counter() - start:.4f} seconds")

    def update_representatives(
            self,
            cell_types: list[str],
            cell_num: Optional[int] = 1000,
            replace: bool = False,
            reset: bool = False
        ) -> None:
        """
        Compute or update representative vias for specified cell types and cache them.

        Args:
            cell_types (list[str]): Cell types to update representatives for.
            cell_num (Optional[int]): Number of cells to use for computing each representative.
            replace (bool): If True, recompute even if representative already exists.
            reset (bool): If True, discard all existing representatives before updating.
        """
        print("Updating representatives...")
        start = time.perf_counter()
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
        print(f"Updating done. Took {time.perf_counter() - start:.4f} seconds")

    def update_boxes(self) -> None:
        """
        Compute and cache the width and height of each cell type's bounding box.
        """
        print("Updating boxes...")
        start = time.perf_counter()
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
        print(f"Updating done. Took {time.perf_counter() - start:.4f} seconds")

    def update_aligned_cells(
            self,
            cell_types: list[str],
            cell_num: Optional[int] = None,
            replace: bool = False,
            reset: bool = False
        ) -> None:
        """
        Align cells to their representatives for the given cell types and cache the result.
        The result contains all cells aligned to their representative and their distances after alignment.

        Args:
            cell_types (list[str]): Cell types to update aligned cells for.
            cell_num (Optional[int]): Number of cells to align per type; if None, use all.
            replace (bool): If True, recompute even if aligned cells already exist.
            reset (bool): If True, discard all existing aligned cells before updating.
        """
        print("Updating aligned cells...")
        start = time.perf_counter()
        aligned_cells = {} if reset else self.get_aligned_cells()
        cell_types = cell_types if replace else [cell_type for cell_type in cell_types if cell_type not in aligned_cells]

        if cell_types:
            with Pool() as pool:
                aligned = pool.starmap(align_cells,
                                       [(self.get_sorted_cells().get(cell_type, [])[:cell_num],
                                         self.get_representatives().get(cell_type, None),
                                         self.get_boxes().get(cell_type, None))
                                         for cell_type in cell_types])
                new_cells = dict(zip(cell_types, aligned))
        else: new_cells = {}
        
        aligned_cells.update(new_cells)
        self._cache_aligned_cells(aligned_cells)
        print(f"Updating done. Took {time.perf_counter() - start:.4f} seconds")
    
    def save_trojans(self, trojans: dict[tuple[int | float, int | float], list[Predicted_Cell]]) -> None:
        """
        Caches cells that were confirmed to be trojans after manual check for a specific box width range.
        These cells can later be retrieved using `get_confirmed_trojans`.

        Args:
            trojans (dict[tuple[int | float, int | float], list[Predicted_Cell]]): Confirmed cells with their actual and predicted label
        """
        current_trojans = self.get_trojans()
        current_trojans.update(trojans)
        path = self.cache_dir / "trojans.pickle"
        with open(path, "wb") as f:
            pickle.dump(current_trojans, f)
    
    def get_trojans(self) -> dict[tuple[int | float, int | float], list[Predicted_Cell]]:
        """Retrieves previously stored confirmed trojans from the cache.

        Returns:
            (dict[tuple[int | float, int | float], list[Predicted_Cell]]): Stored trojans (empty if none were saved)
        """
        path = self.cache_dir / "trojans.pickle"
        if path.exists():
            with open(path, "rb") as f:
                return pickle.load(f)
        return {}


    def group_cells(self, cell_mapping: dict[str, str]) -> None:
        """
        Group cells by a new name mapping and update the cached sorted cells.

        Args:
            cell_mapping (dict[str, str]): Mapping from original cell type to new cell type.
        """
        grouped_cells: dict[str, list[Cell]] = defaultdict(list)
        sorted_cells = self.get_sorted_cells()
        for cell_type, cells in sorted_cells.items():
            for cell in cells:
                cell_map = cell_mapping.get(cell_type, cell_type)
                cell["data"]["name"] = cell_map
                grouped_cells[cell_map].append(cell)
        self._cache_sorted_cells(grouped_cells)


    def get_sorted_cells(self) -> dict[str, list[Cell]]:
        """Return the cached sorted cells."""
        return self.sorted_cells
    def get_representatives(self) -> dict[str, list[tuple[float, float]]]:
        """Return the cached representative vias for each cell type."""
        return self.representatives
    def get_boxes(self) -> dict[str, tuple[float, float]]:
        """Return the cached bounding box dimensions for each cell type."""
        return self.boxes
    def get_aligned_cells(self) -> dict[str, tuple[list[Cell], list[Optional[float]]]]:
        """Return the cached aligned cells and optional alignment metrics."""
        return self.aligned_cells


    # -------- Private functions --------
    def _retrieve_cache(self, name: str) -> dict[str, Any]:
        """
        Load cached data from disk if available.

        Args:
            name (str): Cache key name (e.g. 'sorted_cells', 'boxes').

        Returns:
            dict: Cached data, or an empty dict if no cache exists.
        """
        print(f"Retrieving {name.replace('_', ' ')}")
        path = self.cache_dir / f"{name}.pickle"
        if path.exists():
            with open(path, "rb") as f:
                return pickle.load(f)
        warning = f"No data is cached. Call 'update_{name}()' to cache."
        warnings.warn(warning, UserWarning)
        return {}

    def _retrieve_sorted_cells(self) -> dict[str, list[Cell]]:
        """Retrieve cached sorted cells."""
        return self._retrieve_cache("sorted_cells")
    
    def _retrieve_representatives(self) -> dict[str, list[Point]]:
        """Retrieve cached representative vias."""
        return self._retrieve_cache("representatives")
    
    def _retrieve_boxes(self) -> dict[str, tuple[float, float]]:
        """Retrieve cached cell bounding box sizes."""
        return self._retrieve_cache("boxes")
    
    def _retrieve_aligned_cells(self) -> dict[str, tuple[list[Cell], list[Optional[float]]]]:
        """Retrieve cached aligned cells and alignment scores."""
        return self._retrieve_cache("aligned_cells")


    def _cache_data(self, name: str, data: dict[str, Any]) -> None:
        """
        Persist cached data to disk and update the in-memory attribute.

        Args:
            name (str): Cache key name.
            data (dict): Data to be cached.
        """
        print(f"Saving new {name.replace('_', ' ')}")
        setattr(self, name, data)
        path = self.cache_dir / f"{name}.pickle"
        with open(path, "wb") as f:
            pickle.dump(data, f)

    def _cache_sorted_cells(self, new_cells: dict[str, list[Cell]]) -> None:
        """Cache sorted cells."""
        self._cache_data("sorted_cells", new_cells)

    def _cache_representatives(self, new_reps: dict[str, list[Point]]) -> None:
        """Cache representative vias."""
        self._cache_data("representatives", new_reps)

    def _cache_boxes(self, new_boxes: dict[str, tuple[float, float]]) -> None:
        """Cache cell bounding boxes."""
        self._cache_data("boxes", new_boxes)

    def _cache_aligned_cells(
        self,
        new_cells: dict[str, tuple[list[Cell], list[Optional[float]]]]
    ) -> None:
        """Cache aligned cells and alignment metrics."""
        self._cache_data("aligned_cells", new_cells)
