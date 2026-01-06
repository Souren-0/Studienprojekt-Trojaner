"""
Module: Alignment Utilities <br>
Author: Souren Ishkhanian and Kolja Dorschel

This module provides useful alignment methods for cells and their vias.
"""
import random
import warnings
from tqdm import tqdm
from copy import deepcopy
from collections import Counter
from multiprocessing import Pool
from scipy.spatial import cKDTree as KDTree # type: ignore
# (alternatively use KDTree from sklearn.neighbors, not C optimized)
from Cell_Via_Utilities import diff_pts, add_pts, reset_transform
from Distance_Measures import chamfer
from Annotation_Helpers import *


def sample(lst: list, k: Optional[int]) -> list:
    """
    Return a random ordered subset of elements from a list.

    Args:
        lst (list): Input list.
        k (Optional[int]): Number of elements to sample. If None or >= len(lst), returns a copy of lst.
    """
    if k is None or k >= len(lst):
        return lst[:]
    indices = sorted(random.sample(range(len(lst)), k))
    return [lst[i] for i in indices]


def index(lst: list, idx: list[int]) -> list:
    """
    Index a list with fallback for out-of-range indices.

    Args:
        lst (list): Input list.
        idx (list[int]): Indices to retrieve.
    """
    if not lst: raise ValueError("`lst` must contain at least one element")
    ret = []
    for i in idx:
        if i < len(lst):
            ret.append(lst[i])
        else:
            ret.append(lst[0])
    return ret


def align_to_point(vias : list[Point], point : Point) -> list[Point]:
    """
    Translate a list of points so that each point is shifted relative to a given reference point.

    Args:
        vias (list[Point]): List of points to align.
        point (Point): Reference point to align to.

    Returns:
        list[Point]: New list of points aligned to the reference point.
    """
    return [diff_pts(point, v) for v in vias]


def align_vias_brute_force(cell1_vias : list[Point], cell2_vias : list[Point], itr_count: Optional[int] = 5) -> tuple[list[Point], list[Point]]:
    """
    [DEPRECATED] Use `align_vias` instead.

    Align two sets of vias using a full brute-force approach. Slower and less flexible than `align_vias`.

    Args:
        cell1_vias (list[Point]): Vias of the first cell.
        cell2_vias (list[Point]): Vias of the second cell.
        itr_count (Optional[int]): Number of points to attempt alignment from `cell1_vias`.

    Returns:
        (list[Point], list[Point]): Aligned vias for `cell1_vias` and `cell2_vias`.
    """
    if not (cell1_vias and cell2_vias): return cell1_vias, cell2_vias
    cell1, cell2 = deepcopy(cell1_vias), deepcopy(cell2_vias)
    min_scores = []

    # Align itr_count points to (0, 0)
    alignments_p1 = [(p1, align_to_point(cell1, p1)) for p1 in sample(cell1_vias, itr_count)]

    # Make a list where every point in cell2 is aligned to (0, 0) once.
    alignments_p2 = [(p2, align_to_point(cell2, p2)) for p2 in cell2]

    for p1, p1_aligned in alignments_p1:
        # Score all the alignments bettween cell1 and cell2
        p1_points = KDTree(p1_aligned)
        fitting_scores = [(p1, p2, sum(p1_points.query(p2_points)[0])) for p2, p2_points in alignments_p2]
        min_scores.append(min(fitting_scores, key = lambda x : x[2]))

    p1_alignment, p2_alignment, score = min(min_scores, key = lambda x : x[2])
    aligned_cell1, aligned_cell2 = [diff_pts(p1, p1_alignment) for p1 in cell1], [diff_pts(p2, p2_alignment) for p2 in cell2]
    aligned_cell1, aligned_cell2 = [add_pts(p1, p2_alignment) for p1 in aligned_cell1], [add_pts(p2, p1_alignment) for p2 in aligned_cell2]
    return aligned_cell1, aligned_cell2


def align_vias(
        cell1_vias : list[Point],
        cell2_vias : list[Point],
        itr_count: Optional[int] = None,
        k: int = 3
    ) -> tuple[list[Point], list[Point]]:
    """
    Align two sets of vias (points) as closely as possible.

    Each point in `cell1_vias` is aligned to `cell2_vias` using a heuristic method
    that approximates an optimal alignment without trying every permutation. 
    The alignment tries to minimize the summed Euclidean distances between points
    but may not find the optimal solution.
    For `itr_count = None` and `k = len(cell2_vias)` the function finds the best
    solution with a terrible time-accuracy tradeoff.

    Args:
        cell1_vias (list[Point]): Vias of the first cell.
        cell2_vias (list[Point]): Vias of the second cell.
        itr_count (Optional[int]): Number of points from `cell1_vias` to attempt aligning; 
            use None to consider all points (slower but more accurate).
        k (int): Number of nearest neighbors to consider when scoring alignments (default 3).

    Returns:
        (list[Point], list[Point]): `cell1_vias` aligned to `cell2_vias` and vice versa.
    """
    if not (cell1_vias and cell2_vias): return cell1_vias, cell2_vias
    cell1_vias, cell2_vias = deepcopy(cell1_vias), deepcopy(cell2_vias) # type: ignore
    
    # Align each point to (0,0) once
    alignments_p1 = [(p1, align_to_point(cell1_vias, p1)) for p1 in sample(cell1_vias, itr_count)]
    alignments_p2 = [align_to_point(cell2_vias, p2) for p2 in cell2_vias]
    cell2_vias_KD = KDTree(cell2_vias)
    min_scores = []

    for p1, p1_aligned in alignments_p1:
        # Score all the alignments bettween cell1 and cell2
        p1_points = KDTree(p1_aligned)
        p1_NN = cell2_vias_KD.query(p1, k=k)[1]
        # If k is too large, query will return out of bound values, we map them to p1_NN[0]
        p2_aligned = zip(index(cell2_vias, p1_NN), index(alignments_p2, p1_NN))

        fitting_scores = [(p1, p2, sum(p1_points.query(p2_points)[0])) for p2, p2_points in p2_aligned]
        min_scores.append(min(fitting_scores, key = lambda x : x[2]))

    p1_alignment, p2_alignment, _ = min(min_scores, key = lambda x : x[2])
    aligned_cell1, aligned_cell2 = [diff_pts(p1, p1_alignment) for p1 in cell1_vias], [diff_pts(p2, p2_alignment) for p2 in cell2_vias]
    aligned_cell1, aligned_cell2 = [add_pts(p1, p2_alignment) for p1 in aligned_cell1], [add_pts(p2, p1_alignment) for p2 in aligned_cell2]
    return aligned_cell1, aligned_cell2


def align_cells(
        cells: list[Cell],
        vias: list[Point],
        box: Optional[tuple[int | float, int | float]] = None,
        itr_count: Optional[int] = None,
        distance_measure: Callable[[list[Point], list[Point]], float] = chamfer
    ) -> tuple[list[Cell], list[Optional[float]]]:
    """
    Align a list of cells to a representative set of vias.

    Each cell's vias are aligned to the reference `vias`, optionally filtered by
    a bounding `box`. Returns both the aligned cells and their distances to the
    reference. Distance is computed using `distance_measure`.

    Args:
        cells (list[Cell]): Cells to align.
        vias (list[Point] | None): Reference vias to align to. If None, a draft representative is computed.
        box (Optional[tuple[int | float, int | float]]): Optional bounding box (width, height) to filter aligned vias.
        itr_count (Optional[int]): Number of points to use for alignment heuristics.
        distance_measure (Callable[[list[Point], list[Point]], float]): Function to compute distance after alignment.

    Returns:
        (list[Cell], list[Optional[float]]): Aligned cells and corresponding distances. Distance is `None` if cell has no vias.
    """
    if not vias: return (cells, [None] * len(cells))

    cells = deepcopy(cells)
    # Set vias to a draft representative
    aligned_cells = []
    for cell in tqdm(cells):
        new_vias = align_vias(cell["vias"], vias, itr_count=itr_count)[0]
        if box:
            new_vias = [point for point in new_vias if 0 <= point[0] <= box[0] and 0 <= point[1] <= box[1]]
        cell["vias"] = new_vias
        if cell["vias"]:
            aligned_cells.append((cell, distance_measure(cell["vias"], vias)))
        else:
            warnings.warn("Cannot compute distance of cells with no vias, setting distance to None...")
            aligned_cells.append((cell, None)) # type: ignore

    cells_aligned, distances = zip(*aligned_cells) if aligned_cells else ([], [])
    return list(cells_aligned), list(distances)


def get_aligned_vias(
        cells: Iterable[Cell],
        num_cells: Optional[int] = 100,
        alignment_itr: Optional[int] = 5,
        multiprocess: bool = False
    ) -> tuple[int, list[Point]]:
    """
    Align vias from multiple instances of the same cell type into a common coordinate frame.

    All cells are first normalized to a consistent orientation. A representative starting
    cell is chosen based on the most common via count, and vias from other cells are aligned
    to it and aggregated.

    Args:
        cells (Iterable[Cell]): Collection of cells of the same type.
        num_cells (Optional[int]): Maximum number of cells to use for alignment.
        alignment_itr (Optional[int]): Iterations passed to the alignment routine.
        multiprocess (bool): Whether to parallelize alignment across processes.

    Returns:
        (int, list[Point]): The inferred via count and the list of all aligned vias.
    """
    if not cells: raise ValueError("Cells cannot be empty. Please provide cells to align.")

    # Transform all the cells to the same orientation
    cells = [reset_transform(c, "y") for c in cells]
    cell_type = cells[0]["data"]["name"]

    # Majority vote on the number of vias in the cell
    via_counter = Counter([len(cell["vias"]) for cell in cells])
    via_count = via_counter.most_common(1)[0][0]

    # Choose a good starting cell
    start_cell = None
    for i, cell in enumerate(cells):
        if len(cell["vias"]) == via_count:
            start_cell = cells.pop(i)
            break
    assert start_cell is not None, "Something went wrong. Start cell should not be None."

    # Extract all aligned vias from the cells
    all_vias = deepcopy(start_cell["vias"])
    
    args = [(start_cell["vias"], cell["vias"], alignment_itr) for cell in sample(cells, num_cells)]
    if multiprocess:
        with Pool() as pool:
            results = pool.starmap(align_vias, args)
    else:
        results = [align_vias(*arg) for arg in tqdm(args, desc=cell_type)]

    for _, vias_p2 in results:
        all_vias += vias_p2
    return via_count, all_vias