"""
Module: Cell Via Utilities <br>
Author: Kolja Dorschel <br>
Modified by: Souren Ishkhanian

This module contains functions for manipulating and aligning cell structures, including rotation, reflection, distance calculations, and via alignment.
"""

import math
from copy import deepcopy
from Annotation_Helpers import *


def euk_dist_sq(p1 : Point, p2 : Point) -> float:
    """Compute squared Euclidean distance between two points."""
    return pow(p1[0] - p2[0], 2) + pow(p1[1] - p2[1], 2)


def euk_dist(p1 : Point, p2 : Point) -> float:
    """Compute Euclidean distance between two points."""
    return math.sqrt(euk_dist_sq(p1, p2))


def add_pts(p1 : Point, p2 : Point) -> Point:
    """Add two points component-wise."""
    return p1[0] + p2[0], p1[1] + p2[1]


def diff_pts(p1 : Point, p2 : Point) -> Point:
    """Subtract p2 from p1 component-wise."""
    return p1[0] - p2[0], p1[1] - p2[1]


def rotate_cell_180(cell : Cell) -> Cell:
    """Rotate a cell by 180 degrees around the center of its bounding box."""
    cell = deepcopy(cell)
    box = cell["box"]
    box_width = box[1][0] - box[0][0]
    box_height = box[1][1] - box[0][1]
    # reflect y and reflect x == rotate 180
    cell["vias"] = [(v[0], box_height - v[1]) for v in cell["vias"]]
    cell["vias"] = [(box_width - v[0], v[1]) for v in cell["vias"]]
    return cell


def reflect_cell(cell : Cell, axis : str) -> Cell:
    """
    Reflect a cell's vias across the specified axis within its bounding box.

    Args:
        cell (Cell): Cell to reflect.
        axis (str): Axis to reflect across. Must be either `"x"` or `"y"`.

    Raises:
        AssertionError: If `axis` is not `"x"` or `"y"`.
    """
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
