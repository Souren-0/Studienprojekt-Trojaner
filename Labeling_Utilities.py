"""
Module: Representative & Classification Utilities <br>
Authors: Souren Ishkhanian and Kolja Dorschel

This module provides utilities for computing representative (idealized) via
configurations from noisy cell instances, measuring distances to known
representatives, assigning cell types, and filtering or predicting anomalous
cells (e.g. hardware trojans) based on statistical confidence.
"""

import numpy as np
from scipy.stats import rayleigh
from sklearn.cluster import KMeans
from Alignment_Utilities import *
from Annotation_Helpers import *
from Distance_Measures import jaccard


def find_representative_vias(
        cells: list[Cell],
        num_cells : int = 1000,
        alignment_itr : int = 50,
        filter_itr: int = 2,
        filter_threshold: float = 0.995,
        multiprocess: bool = False
    ) -> list[Point]:
    """
    Compute a representative (idealized) set of vias for a given cell type.

    The representative is obtained by aligning vias from multiple cell
    instances, clustering them using k-means, and iteratively filtering
    outliers based on a Rayleigh-distributed distance model.

    Args:
        cells (list[Cell]): Cells of the same type.
        num_cells (int): Maximum number of cells to use for alignment.
        alignment_itr (int): Iterations used during via alignment.
        filter_itr (int): Number of outlier-filtering iterations.
        filter_threshold (float): Confidence threshold for outlier removal.
        multiprocess (bool): Whether to use multiprocessing during alignment.

    Returns:
        list[Point]: Representative via positions ("perfect cell").
    """
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
        for _ in range(filter_itr):
            # Compute the distances from each point to their label point
            distances = np.linalg.norm(filtered_vias - kmeans.cluster_centers_[kmeans.labels_], axis=1)

            # Fit the distances to a rayleigh distribution
            distances[0] += np.finfo(float).eps # Avoids overflow warning when all distances are of equal value
            loc, scale = rayleigh.fit(distances + np.finfo(float).eps) # Avoids divison by 0 when a distance is 0
            
            # Filter points outside of threshold. Filter option can be adapted
            filtered_vias = filtered_vias[distances < rayleigh.ppf(filter_threshold, loc=loc, scale=scale)] 
            kmeans.fit(filtered_vias)
        representative = [tuple(point) for point in kmeans.cluster_centers_]
    return representative


def cell_to_labels_dists(cell: Cell, representatives: dict[str, list[Point]]) -> dict[str, float]:
    """
    Compute distances between a cell and multiple representative labels.

    The cell vias are aligned to each representative before computing
    the Chamfer distance.

    Args:
        cell (Cell): Cell to evaluate.
        representatives (dict[str, list[Point]]): Mapping of labels to representative vias.

    Returns:
        dict[str, float]: Distance from the cell to each representative label.
    """
    cell_vias = cell["vias"]
    dists = {}
    for representative, vias in representatives.items():
        if not vias:
            dists[representative] = float('inf')
            continue
        aligned_cell = align_vias(cell_vias, vias, itr_count=None)[0]
        dists[representative] = chamfer(aligned_cell, vias)
    return dists


def assign_cell_type(cell: Cell, representatives: dict[str, list[Point]], bias_strength: float = 0.5) -> str:
    """
    Assign a cell to the most likely representative label.

    Optionally biases the assignment toward the cell's original label by
    requiring competing labels to be sufficiently better.

    Args:
        cell (Cell): Cell to classify.
        representatives (dict[str, list[Point]]): Known representative vias.
        bias_strength (float): Strength of bias toward the original label.

    Returns:
        str: Assigned cell type label.
    """
    dists = cell_to_labels_dists(cell, representatives)
    favored_label = cell["data"]["name"]
    if favored_label not in dists:
        warnings.warn("Cannot label biased. Cell type not in representatives")
    thr = dists[favored_label] if favored_label in dists else float('inf')

    current_label = (favored_label, thr)
    for representative, dist in dists.items():
        if dist < thr * (1 - bias_strength) and dist < current_label[1]:
            current_label = (representative, dist)
    return current_label[0]


def filter_cells(cells: list[Cell], dists: list[Optional[float]], confidence_threshold: float = 0.9999) -> list[Cell]:
    """
    Filter cells based on statistical confidence of their distance scores.

    Distances are modeled using a Rayleigh distribution, and cells exceeding
    the specified confidence threshold are returned.

    Args:
        cells (list[Cell]): Cells corresponding to the distances.
        dists (list[float]): Distance values for each cell.
        confidence_threshold (float): Rayleigh confidence threshold.

    Returns:
        list[Cell]: Cells classified as outliers.
    """
    cell_dists = [(c, d) for c, d in zip(cells, dists) if d is not None]
    
    if len(cell_dists) < 2:
        return []
    
    cells_np, dists_np = zip(*cell_dists)
    cells_np, dists_np = np.array(cells_np), (np.array(dists_np) + np.finfo(float).eps) # type: ignore
    dists_np[0] += np.finfo(float).eps # type: ignore

    loc, scale = rayleigh.fit(dists_np)

    mask = dists_np > rayleigh.ppf(confidence_threshold, loc=loc, scale=scale)
    return list(cells_np[mask])


def predict_trojans(cells: list[Cell], dists: list[Optional[float]], representatives: dict[str, list[Point]], bias_strength: float = 0.5) -> list[Predicted_Cell]:
    """
    Predict potential trojan cells by combining distance-based filtering
    and representative-based reclassification.

    A cell is flagged if it is a statistical outlier and its predicted
    label differs from its original label.

    Args:
        cells (list[Cell]): Cells to analyze.
        dists (list[Optional[float]]): Precomputed distance scores.
        representatives (dict[str, list[Point]]): Known representative vias.

    Returns:
        list[Predicted_Cell]: Cells predicted to be trojans with labels.
    """
    ret = []
    possible_trojans = filter_cells(cells, dists, 0.9)
    for trojan in tqdm(possible_trojans):
        predicted = assign_cell_type(trojan, representatives, bias_strength=bias_strength)
        actual = trojan["data"]["name"]
        if predicted != actual:
            ret.append(Predicted_Cell(cell=trojan, actual=actual, predicted=predicted))
    return ret


def prune_predicted_trojans(cells: list[Predicted_Cell], thr: int = 6) -> list[Predicted_Cell]:
    def get_key(cell: Predicted_Cell) -> tuple[str, str]:
        return (cell['actual'], cell['predicted'])
    
    confusion_counts: Counter[tuple[str, str]] = Counter()
    for cell in cells:
        confusion_counts[get_key(cell)] += 1
    
    ret: list[Predicted_Cell] = []
    for cell in cells:
        if confusion_counts[get_key(cell)] <= thr: ret.append(cell)

    return ret


def confirm_predicted_trojans(cells: list[Predicted_Cell], dists: dict[str, list[float]], representatives: dict[str, list[Point]]) -> list[Predicted_Cell]:
    def get_jaccard_radius(dists: list[float], thr: float = 0.5) -> float:
            if not dists: return 1
            dists_np = np.array(dists)
            dists_np[0] += np.finfo(float).eps # Avoids overflow warning when all distances are of equal value
            loc, scale = rayleigh.fit(dists_np + np.finfo(float).eps) # Avoids divison by 0 when a distance is 0
            return max(1, float(rayleigh.ppf(thr, loc=loc, scale=scale)))
    
    ret: list[Predicted_Cell] = []
    radii: dict[str, float] = {}
    for cell in cells:
        cell_vias = cell["cell"]["vias"]
        actual_vias = align_vias(representatives[cell["actual"]], cell_vias)[0]
        actual_radius = radii.setdefault(cell["actual"], get_jaccard_radius(dists[cell["actual"]]))
        predicted_vias = align_vias(representatives[cell["predicted"]], cell_vias)[0]
        predicted_radius = radii.setdefault(cell["predicted"], get_jaccard_radius(dists[cell["predicted"]]))
        if jaccard(cell_vias, predicted_vias, predicted_radius) <= jaccard(cell_vias, actual_vias, actual_radius):
            ret.append(cell)

    return ret
