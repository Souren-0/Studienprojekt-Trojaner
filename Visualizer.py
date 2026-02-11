"""
Module: Visualizer

Author: Souren Ishkhanian

This module helps visualizing cells by plotting their vias in their box.
"""

import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from Annotation_Helpers import *
from typing import cast

plt.style.use('dark_background')

colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
n = len(colors)
colors = {
    'line_color' : colors[1 % n],
    'via_color' : colors[5 % n],
    'cell_color' : colors[0 % n],
    'representative_color' : "red"
}

class Visualizer:
    """
    A class for visualizing cells, their vias, and the representative "perfect cell" within a bounding box.

    Attributes:
        cells (list[Cell]): List of cells to visualize.
        representative (Optional[list[Point]]): Points representing the ideal via positions of a cell.
        vias (Optional[list[Point]]): Additional vias to display.
        box (Rectangle): Matplotlib rectangle representing the plotting area.
        fig (Figure): Matplotlib figure.
        ax (Axes): Matplotlib axes.
    """

    def __init__(self,
                 cells: Optional[list[Cell]] = None,
                 representative: Optional[list[Point]] = None,
                 vias: Optional[list[Point]] = None,
                 box: Optional[tuple[float, float] | tuple[Point, Point]] = None) -> None:
        """
        Initialize a Visualizer with optional cells, representative vias, additional vias, and a box.

        Args:
            cells (Optional[list[Cell]]): Initial list of cells.
            representative (Optional[list[Point]]): Points showing the ideal via positions of a cell.
            vias (Optional[list[Point]]): Additional vias to plot.
            box (Optional[tuple]): Either (width, height) or ((top_left), (bottom_right)) defining the plotting box.
        """
        self.cells: list[Cell] = cells or []
        self.representative = representative
        self.vias = vias

        if box and isinstance(box[0], (int, float)):
            width, height = cast(tuple[float, float], box)
        else:
            box = cast(tuple[Point, Point], box) if box \
                else (cells[0]['box'] if cells else ((0,0), (0,0)))
            top_left, bottom_right = box
            width = bottom_right[0] - top_left[0]
            height = bottom_right[1] - top_left[1]

        self.box = Rectangle((0, 0), width, height, facecolor='none', edgecolor=colors['line_color'])

        self.fig, self.ax = plt.subplots(figsize=(12,9))
        self.ax.add_patch(self.box)
        self.ax.set_xlim(-0.25*width, 1.25*width)
        self.ax.set_ylim(-0.25*height, 1.25*height)


    def add_cells(self, cell: list[Cell]) -> None:
        """
        Add more cells to the visualizer.

        Args:
            cell (list[Cell]): Cells to add.
        """
        self.cells.extend(cell)


    def display_all(self, legend: bool = True) -> None:
        """
        Plot all elements: cells, vias, and representative points, then show the figure.

        Args:
            legend (bool): Whether to display a legend.
        """
        self.display_cells()
        self.display_vias()
        self.display_representative()
        self.show(legend)
    

    def display_vias(self) -> None:
        """
        Plot the additional vias on the figure.
        """
        x_val, y_val = zip(*self.vias) if self.vias else ([], [])
        sns.scatterplot(x=x_val, y=y_val, s=3, color=colors['via_color'], edgecolor='none', label="Additional vias", ax=self.ax)


    def _display_cell(self, cell: Cell) -> None:
        """
        Plot the vias of a single cell. Private helper method.

        Args:
            cell (Cell): The cell whose vias are to be plotted.
        """
        vias = cell.get('vias', [])
        if vias:
            x_vals, y_vals = zip(*vias)
            sns.scatterplot(x=x_vals, y=y_vals, s=5, color=colors['cell_color'], edgecolor='none', alpha=0.3, ax=self.ax)


    def display_cells(self) -> None:
        """
        Plot all cells' vias on the figure. The first cell's vias are labeled for the legend.
        """
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
        """
        Plot the representative vias on the figure, if any.
        """
        if self.representative:
            rep_x, rep_y = zip(*self.representative)
            sns.scatterplot(x=rep_x, y=rep_y, s=10, color=colors['representative_color'], edgecolor='none', label="Representative", ax=self.ax)
    

    def show(self, legend: bool = True) -> None:
        """
        Render the figure. Optionally display a legend.

        Args:
            legend (bool): Whether to show a legend.
        """
        if legend: self.ax.legend()
        plt.show()
