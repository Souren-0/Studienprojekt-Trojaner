"""
Module: Distance Measures
Author: Souren Ishkhanian

This module contains the following distance measures between point clouds:
- Hausdorff
- Chamfer
- Earth Mover
- Jaccard

First three measures can be a normalized in range [0, 1] if a box (width, height) is provided
"""

import numpy as np
import torch
import ot
from shapely.geometry import Point
from shapely.ops import unary_union
from scipy.spatial.distance import directed_hausdorff

from typing import Any, Optional, Sequence
from numpy.typing import ArrayLike


def _pad(arr1: Sequence[Any], arr2: Sequence[Any], pad_pattern: Any) -> tuple[list, list]:
    arr1_len, arr2_len = len(arr1), len(arr2)
    max_len = max(arr1_len, arr2_len)
    arr1 = list(arr1) + [pad_pattern] * (max_len - arr1_len)
    arr2 = list(arr2) + [pad_pattern] * (max_len - arr2_len)
    return arr1, arr2


def _box_diagonal(width: float, height: float) -> float:
    return (width**2 + height**2)**0.5


def _distance_matrix(points1: Sequence[tuple[float, float]], points2: Sequence[tuple[float, float]]) -> np.ndarray:
    return np.array([[np.linalg.norm(np.subtract(p1, p2)) for p2 in points2] for p1 in points1])


def _hausdorff(points1: ArrayLike, points2: ArrayLike) -> float:
    return max(directed_hausdorff(points1, points2)[0], directed_hausdorff(points2, points1)[0])


def hausdorff(points1: ArrayLike, points2: ArrayLike, box: Optional[tuple[float, float]] = None) -> float:
    norm = _box_diagonal(*box) if box else 1
    return _hausdorff(points1, points2) / norm


def chamfer(points1: ArrayLike, points2: ArrayLike, box: Optional[tuple[float, float]] = None, directed: bool = True) -> float:
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
    norm = _box_diagonal(*box) if box else 1

    points1_weights = np.array([1 / len(points1) for _ in range(len(points1))])
    points2_weights = np.array([1 / len(points2) for _ in range(len(points2))])

    vias = _pad(points1, points2, (0, 0))
    weights = _pad(points1_weights, points2_weights, 0) # type: ignore
    dist_matrix = _distance_matrix(*vias)
    flow_matrix = ot.emd(*weights, dist_matrix)
    return np.sum(dist_matrix * flow_matrix) / norm


def jaccard(points1: Sequence[tuple[float, float]], points2: Sequence[tuple[float, float]], radius: float) -> float:
    disks1 = unary_union([Point(x, y).buffer(radius) for x, y in points1])
    disks2 = unary_union([Point(x, y).buffer(radius) for x, y in points2])
    
    return 1 - (disks1.intersection(disks2).area / disks1.union(disks2).area)
