"""
Module: Cell Via Utilities
Author: Lukas Plätz
Modified by: Souren Ishkhanian

This module contains functions for manipulating and aligning cell structures, including rotation, reflection, distance calculations, and via alignment.
"""

import math
import numpy as np
from tqdm import tqdm
from copy import deepcopy
from collections import Counter
from Annotation_Helpers import *
from typing import Iterable, Optional

from sklearn.cluster import KMeans
from scipy.stats import rayleigh
from scipy.spatial import cKDTree as KDTree # pyright: ignore[reportAttributeAccessIssue] (alternatively use KDTree from sklearn.neighbors, not C optimized)

from multiprocessing import Pool
from Visualizer import Visualizer
from Distance_Measures import chamfer


def euk_dist_sq(p1 : Point, p2 : Point) -> float:
    return pow(p1[0] - p2[0], 2) + pow(p1[1] - p2[1], 2)


def euk_dist(p1 : Point, p2 : Point) -> float:
    return math.sqrt(euk_dist_sq(p1, p2))


def add_pts(p1 : Point, p2 : Point) -> Point:
    return p1[0] + p2[0], p1[1] + p2[1]


def diff_pts(p1 : Point, p2 : Point) -> Point:
    return p1[0] - p2[0], p1[1] - p2[1]


def rotate_cell_180(cell : Cell) -> Cell:
    cell = deepcopy(cell)
    box = cell["box"]
    box_width = box[1][0] - box[0][0]
    box_height = box[1][1] - box[0][1]
    # reflect y and reflect x == rotate 180
    cell["vias"] = [(v[0], box_height - v[1]) for v in cell["vias"]]
    cell["vias"] = [(box_width - v[0], v[1]) for v in cell["vias"]]
    return cell


def reflect_cell(cell : Cell, axis : str) -> Cell:
    axis = axis.lower()
    assert axis in ["x", "y"], "Invalid axis given."

    cell = deepcopy(cell)
    box = cell["box"]
    box_width = box[1][0] - box[0][0]
    box_height = box[1][1] - box[0][1]
    if axis == "y":
        cell["vias"] = [(v[0], box_height - v[1]) for v in cell["vias"]] # powerline direction Y
    else:
        cell["vias"] = [(box_width - v[0], v[1]) for v in cell["vias"]] # powerline direction X_filtered
    return cell


def reset_transform(cell : Cell, powerline_direction : str = 'y') -> Cell:
    """Transform the cell back to rotation = 0 and no reflection."""

    if cell["data"]["rotation"] == 180 and cell["data"]["reflection"]:
        cell = reflect_cell(cell, "x" if powerline_direction == "y" else "y")
    elif cell["data"]["rotation"] == 180:
        cell =  rotate_cell_180(cell)
    elif cell["data"]["reflection"]:
        cell = reflect_cell(cell, powerline_direction)
    
    cell["data"]["rotation"] = 0
    cell["data"]["reflection"] = False
    return cell


def score_matching(pts1 : list[Point], pts2 : list[Point]) -> float:
    """Try to match the set of points 1 to points 2. The points of set 1 are matched to the closest
    point in set 2. The resulting score is the sum of the distances"""

    min_dists = []
    for p1 in pts1:
        dists = [euk_dist_sq(p1, p2) for p2 in pts2]
        min_dists.append(min(dists))
    return sum(min_dists)


def align_to_point(vias : list[Point], point : Point) -> list[Point]:
    return [diff_pts(point, v) for v in vias]


def align_cells(cell1_vias : list[Point], cell2_vias : list[Point], itr_count: Optional[int] = 5) -> tuple[list[Point], list[Point]]:
    """ Given to sets of vias (points) try to align them as well as possible. For the optimal results all points
    would have to be tested but in all cases 5 works well enough. To use all points set itr_count to `None`."""
    if not (cell1_vias and cell2_vias): return cell1_vias, cell2_vias
    cell1, cell2 = deepcopy(cell1_vias), deepcopy(cell2_vias)
    min_scores = []
    
    for p1 in cell1[:itr_count]:
        # Set the chosen via in cell1 to (0, 0).
        p1_aligned_cell1 = align_to_point(cell1, p1)

        # Make a list where every point in cell2 is aligned to (0, 0) once.
        alignments_p2 = [(p2, align_to_point(cell2, p2)) for p2 in cell2]

        # Score all the alignments bettween cell1 and cell2
        p1_points = KDTree(p1_aligned_cell1)
        fitting_scores = [(p1, p2, sum(p1_points.query(p2_points)[0])) for p2, p2_points in alignments_p2]

        min_scores.append(min(fitting_scores, key = lambda x : x[2]))
    p1_alignment, p2_alignment, score = min(min_scores, key = lambda x : x[2])
    aligned_cell1, aligned_cell2 = [diff_pts(p1, p1_alignment) for p1 in cell1], [diff_pts(p2, p2_alignment) for p2 in cell2]
    aligned_cell1, aligned_cell2 = [add_pts(p1, p2_alignment) for p1 in aligned_cell1], [add_pts(p2, p1_alignment) for p2 in aligned_cell2]
    return aligned_cell1, aligned_cell2


def get_aligned_vias(cells: Iterable[Cell],
                     num_cells: Optional[int] = 100,
                     alignment_itr: Optional[int] = 5,
                     multiprocess: bool = False) -> tuple[int, list[Point]]:
    """Choose a cell type and get a list of all the aligned vias from a chosen number of instances. Nested multiprocessing is not recommended."""
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

    if multiprocess:
        with Pool() as pool:
            results = pool.starmap(align_cells, [(start_cell["vias"], cell["vias"], alignment_itr) for cell in cells[:num_cells]])
    else:
        results = [align_cells(start_cell["vias"], cell["vias"], itr_count=alignment_itr) for cell in tqdm(cells[:num_cells], desc=cell_type)]

    for _, vias_p2 in results:
        all_vias += vias_p2
    return via_count, all_vias


def find_representative_vias(cells: Iterable[Cell],
                             num_cells : Optional[int] = 1000,
                             alignment_itr : Optional[int] = 50,
                             plot: bool = False,
                             multiprocess: bool = False) -> list[Point]:
    via_count, cell_vias = get_aligned_vias(cells, num_cells=num_cells, alignment_itr=alignment_itr, multiprocess=multiprocess)

    cells = list(cells)
    filtered_vias = np.array(deepcopy(cell_vias))
    assert cells, "Something went wrong. Cells cannot be empty. 'get_aligned_vias' should've catched that error."

    if via_count < 1: representative = []
    elif len(cells) == 1: representative = cells[0]["vias"]
    else:
        # Predict representative's via count using the most frequent count (approach can be adapted)
        kmeans = KMeans(n_clusters=via_count)
        kmeans.fit(cell_vias)
        # Remove outliers using confidence intervals
        for _ in range(2):
            # Compute the distances from each point to their label point
            distances = np.linalg.norm(filtered_vias - kmeans.cluster_centers_[kmeans.labels_], axis=1)

            # Fit the distances to a rayleigh distribution
            distances[0] += np.finfo(float).eps # Avoids overflow warning when all distances are of equal value
            loc, scale = rayleigh.fit(distances + np.finfo(float).eps) # Avoids divison by 0 when a distance is 0
            
            # Filter points outside of threshold. Filter option can be adapted
            threshold = 0.99
            filtered_vias = filtered_vias[distances < rayleigh.ppf(threshold, loc=loc, scale=scale)] 
            kmeans.fit(filtered_vias)
        representative = [tuple(pt) for pt in kmeans.cluster_centers_]
    
    if plot:
        print("Plotting...")
        visualizer = Visualizer(cells[0]["box"], cells, representative, list(filtered_vias))
        visualizer.display_all()
    return representative


def assign_cell_type(cell: Cell, representatives: dict[str, list[Point]]) -> str:
    #TODO: Align cells with reps first
    cell_vias = cell["vias"]
    dists = [(representative, chamfer(cell_vias, representatives[representative])) for representative in representatives]
    return min(dists, key=lambda x: x[1])[0]