"""
Module: Visualizer
Author: Souren Ishkhanian

This module helps visualizing cells by plotting their vias in their box.
"""

from typing import Optional, cast
from matplotlib.patches import Rectangle
from Annotation_Helpers import *
import matplotlib.pyplot as plt
import seaborn as sns

plt.style.use('dark_background')

colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
n = len(colors)
colors = {
    'line_color' : colors [1 % n],
    'via_color' : colors[5 % n],
    'cell_color' : colors[0 % n],
    'representative_color' : "red"
}

class Visualizer:
    def __init__(self,
                 box: tuple[float, float] | tuple[Point, Point],
                 cells: list[Cell] = [],
                 representative: Optional[list[Point]] = None,
                 vias: Optional[list[Point]] = None) -> None:
        self.cells: list[Cell] = cells
        self.representative = representative
        self.vias = vias

        if isinstance(box[0], (int, float)):
            width, height = cast(tuple[float, float], box)
        else:
            box = cast(tuple[Point, Point], box)
            top_left, bottom_right = box
            width = bottom_right[0] - top_left[0]
            height = bottom_right[1] - top_left[1]

        self.box = Rectangle((0, 0), width, height, facecolor='none', edgecolor=colors['line_color'])

        self.fig, self.ax = plt.subplots(figsize=(12,9))
        self.ax.add_patch(self.box)
        self.ax.set_xlim(-0.25*width, 1.25*width)
        self.ax.set_ylim(-0.25*height, 1.25*height)


    def add_cells(self, cell: list[Cell]) -> None:
        self.cells.extend(cell)


    def display_all(self, legend: bool = True) -> None:
        self.display_cells()
        self.display_vias()
        self.display_representative()
        self.show(legend)
    

    def display_vias(self) -> None:
        x_val, y_val = zip(*self.vias) if self.vias else (None, None)
        sns.scatterplot(x=x_val, y=y_val, s=3, color=colors['via_color'], edgecolor='none', label="Additional vias", ax=self.ax)


    def _display_cell(self, cell: Cell) -> None:
        vias = cell.get('vias', [])
        if vias:
            x_vals, y_vals = zip(*vias)
            sns.scatterplot(x=x_vals, y=y_vals, s=5, color=colors['cell_color'], edgecolor='none', alpha=0.3, ax=self.ax)


    def display_cells(self) -> None:
        if not self.cells:
            return

        # First cell with legend
        first_cell = self.cells[0]
        vias = first_cell.get('vias', [])
        if vias:
            x_vals, y_vals = zip(*vias)
            sns.scatterplot(x=x_vals, y=y_vals, s=5, color=colors['cell_color'], edgecolor='none', alpha=0.3, label='All cell vias', ax=self.ax)
            
        for cell in self.cells:
            self._display_cell(cell)


    def display_representative(self) -> None:
        if self.representative:
            rep_x, rep_y = zip(*self.representative)
            sns.scatterplot(x=rep_x, y=rep_y, s=10, color=colors['representative_color'], edgecolor='none', label="Representative", ax=self.ax)
    

    def show(self, legend: bool = True) -> None:
        if legend: self.ax.legend()
        plt.show()
