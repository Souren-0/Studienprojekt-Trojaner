"""
Module: Visualizer
Author: Souren Ishkhanian

This module helps visualizing cells by plotting their vias in their box.
"""

from typing import Optional
from matplotlib.patches import Rectangle
from Annotation_Helpers import *
import matplotlib.pyplot as plt
import seaborn as sns

class Visualizer:
    def __init__(self, box: tuple[Point, Point], cells: list[Cell] = [], representative: Optional[list[Point]] = None) -> None:
        self.cells: list[Cell] = cells
        self.representative = representative

        top_left, bottom_right = box
        width = bottom_right[0] - top_left[0]
        height = bottom_right[1] - top_left[1]
        self.box = Rectangle((0, 0), width, height, facecolor='none', edgecolor='black')

        self.fig, self.ax = plt.subplots(figsize=(12,9))
        self.ax.add_patch(self.box)
        self.ax.set_xlim(-0.25*width, 1.25*width)
        self.ax.set_ylim(-0.25*height, 1.25*height)

    def add_cells(self, cell: list[Cell]) -> None:
        self.cells.extend(cell)

    def display_cell(self, cell: Cell) -> None:
        vias = cell.get('vias', [])
        if vias:
            x_vals, y_vals = zip(*vias)
            sns.scatterplot(x=x_vals, y=y_vals, s=5, color='blue', edgecolor='none', alpha=0.5, ax=self.ax)

    def display_all(self) -> None:
        for cell in self.cells:
            self.display_cell(cell)

        if self.representative:
            rep_x, rep_y = zip(*self.representative)
            sns.scatterplot(x=rep_x, y=rep_y, s=20, color='red', edgecolor='none', ax=self.ax)

        plt.show()
        