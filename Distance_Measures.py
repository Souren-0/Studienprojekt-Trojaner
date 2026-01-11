"""
Module: Distance Measures <br>
Author: Souren Ishkhanian

This module contains the following distance measures between point clouds:
- Hausdorff
- Chamfer
- Earth Mover
- Jaccard

First three measures can be a normalized in range [0, 1] if a box (width, height) is provided
"""

import ot
import torch
import numpy as np
from shapely.geometry import Point
from shapely.ops import unary_union
from scipy.spatial.distance import directed_hausdorff

from typing import Any, Optional, Sequence
from numpy.typing import ArrayLike


def _pad(arr1: Sequence[Any], arr2: Sequence[Any], pad_pattern: Any) -> tuple[list, list]:
    """
    Pad two sequences to equal length with a specified pattern.

    Args:
        arr1 (Sequence[Any]): First sequence.
        arr2 (Sequence[Any]): Second sequence.
        pad_pattern (Any): Value used for padding.

    Returns:
        tuple[list, list]: Padded sequences of equal length.
    """
    arr1_len, arr2_len = len(arr1), len(arr2)
    max_len = max(arr1_len, arr2_len)
    arr1 = list(arr1) + [pad_pattern] * (max_len - arr1_len)
    arr2 = list(arr2) + [pad_pattern] * (max_len - arr2_len)
    return arr1, arr2


def _box_diagonal(width: float, height: float) -> float:
    """
    Compute the diagonal length of a box.

    Args:
        width (float): Width of the box.
        height (float): Height of the box.

    Returns:
        float: Euclidean diagonal of the box.
    """
    return (width**2 + height**2)**0.5


def _distance_matrix(points1: Sequence[tuple[float, float]], points2: Sequence[tuple[float, float]]) -> np.ndarray:
    """
    Compute the pairwise Euclidean distance matrix between two sets of points.

    Args:
        points1 (Sequence[tuple[float, float]]): First point set.
        points2 (Sequence[tuple[float, float]]): Second point set.

    Returns:
        np.ndarray: Distance matrix of shape (len(points1), len(points2)).
    """
    return np.array([[np.linalg.norm(np.subtract(p1, p2)) for p2 in points2] for p1 in points1])


def _hausdorff(points1: ArrayLike, points2: ArrayLike) -> float:
    """
    Compute the Hausdorff distance between two point sets.

    Args:
        points1 (ArrayLike): First point set.
        points2 (ArrayLike): Second point set.

    Returns:
        float: Hausdorff distance.
    """
    return max(directed_hausdorff(points1, points2)[0], directed_hausdorff(points2, points1)[0])


def hausdorff(points1: ArrayLike, points2: ArrayLike, box: Optional[tuple[float, float]] = None) -> float:
    """
    Compute the (optionally normalized) Hausdorff distance between two point sets.

    Args:
        points1 (ArrayLike): First point set.
        points2 (ArrayLike): Second point set.
        box (Optional[tuple[float, float]]): If provided, normalize by box diagonal.

    Returns:
        float: Normalized Hausdorff distance in [0, 1].
    """
    norm = _box_diagonal(*box) if box else 1
    return _hausdorff(points1, points2) / norm


def chamfer(points1: ArrayLike, points2: ArrayLike, box: Optional[tuple[float, float]] = None, directed: bool = False) -> float:
    """
    Compute the Chamfer distance between two point sets.

    Args:
        points1 (ArrayLike): First point set.
        points2 (ArrayLike): Second point set.
        box (Optional[tuple[float, float]]): If provided, normalize by box diagonal.
        directed (bool): If True, compute directed Chamfer distance; else, average both directions.

    Returns:
        float: Normalized Chamfer distance in [0, 1].
    """
    norm = _box_diagonal(*box) if box else 1

    dist = torch.cdist(
        torch.tensor(points1, dtype=torch.float32).unsqueeze(0),
        torch.tensor(points2, dtype=torch.float32).unsqueeze(0))
    
    if directed:
        chamfer = dist.min(2)[0].mean() 
    else:
        chamfer = (dist.min(2)[0].mean() + dist.min(1)[0].mean()) / 2 
    return chamfer.item() / norm


def emd(points1: Sequence[tuple[float, float]], points2: Sequence[tuple[float, float]], box: Optional[tuple[float, float]] = None) -> float:
    """
    Compute the Earth Mover's Distance (EMD) between two point sets.

    Args:
        points1 (Sequence[tuple[float, float]]): First point set.
        points2 (Sequence[tuple[float, float]]): Second point set.
        box (Optional[tuple[float, float]]): If provided, normalize by box diagonal.

    Returns:
        float: Normalized EMD.
    """
    norm = _box_diagonal(*box) if box else 1

    points1_weights = np.array([1 / len(points1) for _ in range(len(points1))])
    points2_weights = np.array([1 / len(points2) for _ in range(len(points2))])

    vias = _pad(points1, points2, (0, 0))
    weights = _pad(points1_weights, points2_weights, 0) # type: ignore
    dist_matrix = _distance_matrix(*vias)
    flow_matrix = ot.emd(*weights, dist_matrix)
    return np.sum(dist_matrix * flow_matrix) / norm


def jaccard(points1: Sequence[tuple[float, float]], points2: Sequence[tuple[float, float]], radius: float) -> float:
    """
    Compute the Jaccard similarity between two point sets using disks of a given radius.

    Args:
        points1 (Sequence[tuple[float, float]]): First point set.
        points2 (Sequence[tuple[float, float]]): Second point set.
        radius (float): Radius of disks around points.

    Returns:
        float: 1 minus the Jaccard index (distance-like measure).
    """
    disks1 = unary_union([Point(x, y).buffer(radius) for x, y in points1])
    disks2 = unary_union([Point(x, y).buffer(radius) for x, y in points2])
    
    return 1 - (disks1.intersection(disks2).area / disks1.union(disks2).area)


def directed_jaccard(points1: Sequence[tuple[float, float]], points2: Sequence[tuple[float, float]], radius: float) -> float:
    """
    Compute the directed Jaccard from `points1` to `points2` using disks of a given radius.
    This way, outliers of `points1` are mostly ignored because only the disks of `points2` need to be "filled".

    Args:
        points1 (Sequence[tuple[float, float]]): First point set.
        points2 (Sequence[tuple[float, float]]): Second point set.
        radius (float): Radius of disks around points.

    Returns:
        float: 1 minus the Jaccard index (distance-like measure).
    """
    disks1 = unary_union([Point(x, y).buffer(radius) for x, y in points1])
    disks2 = unary_union([Point(x, y).buffer(radius) for x, y in points2])

    return 1 - (disks1.intersection(disks2).area / disks2.area)
