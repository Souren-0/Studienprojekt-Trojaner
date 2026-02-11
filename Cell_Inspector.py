"""
Module: Cell Inspector

Author: Souren Ishkhanian

This module helps manually checking for trojans.
"""

import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.backend_bases import Event, KeyEvent
from Alignment_Utilities import align_vias
from Annotation_Helpers import *
from typing import cast

plt.style.use('dark_background')

colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
n = len(colors)
colors = {
    'line_color' : colors[1 % n],
    'cell_color' : colors[0 % n],
    'representative_color' : "red"
}

class Cell_Inspector:
    """
    A class for interactively inspecting predicted cells and identifying trojans.

    Attributes:
        predicted_cells (list[Predicted_Cell]): Cells to be manually inspected.
        trojans (list[Predicted_Cell]): Cells confirmed as trojans.
        representatives (dict[str, list[Point]]): Representative "perfect cell" templates keyed by label.
        distance (Callable[[list[Point], list[Point]], float]): Function to compute distance between vias.
        rep_mode (bool): Whether to use actual or predicted label for distance comparison.
        i (int): Index of the current cell being inspected.
        fig (Figure): Matplotlib figure.
        ax (Axes): Matplotlib axes.
    """

    def __init__(self,
                 predicted_cells: list[Predicted_Cell],
                 trojans: list[Predicted_Cell],
                 representatives: dict[str, list[Point]],
                 distance_measure: Callable[[list[Point], list[Point]], float]) -> None:
        """
        Initialize the Cell_Inspector with predicted cells, mutable list to save trojans, representatives, and a distance function.

        Args:
            predicted_cells (list[Predicted_Cell]): Cells to inspect.
            trojans (list[Predicted_Cell]): List to store confirmed trojans.
            representatives (dict[str, list[Point]]): Representative via positions for each label.
            distance_measure (Callable[[list[Point], list[Point]], float]): Function to compute distance between two sets of vias.
        """
        
        self.predicted_cells = predicted_cells
        self.representatives = representatives
        self.rep_mode = True
        self.i = 0
        self.zoom = 0.2

        self.distance = distance_measure
        self.trojans = trojans

        self.fig, self.ax = plt.subplots(figsize=(12,9))
        self._show_current()


    def _get_current(self) -> Optional[Predicted_Cell]:
        """
        Get the current cell being inspected.

        Returns:
            Optional[Predicted_Cell]: The current predicted cell, or None if no cells remain.
        """
        n = len(self.predicted_cells)
        self.i = (self.i % n) if n > 0 else 0
        return self.predicted_cells[self.i] if self.predicted_cells else None

    
    def _confirm_trojan(self, confirmed: bool) -> None:
        """
        Confirm or discard the current cell as a trojan.

        Args:
            confirmed (bool): If True, mark the current cell as a trojan.
        """
        current = self._get_current()
        if confirmed and current:
            self.trojans.append(current)
        del self.predicted_cells[self.i]


    def _show_current(self) -> None:
        """
        Display the current cell, its vias, representative alignments, and distance metrics on the plot.
        """
        self.ax.cla()
        current = self._get_current()
        if current:
            self._make_box(current)
            cell_vias = current["cell"]["vias"]
            actual_label = current["actual"]
            predicted_label = current["predicted"]
            actual_vias = align_vias(self.representatives[current["actual"]], cell_vias)[0]
            predicted_vias = align_vias(self.representatives[current["predicted"]], cell_vias)[0]
            
            palette = [self.ax.get_facecolor(), colors["representative_color"]]
            self._scatter_vias(actual_vias,    f"Actual Label:    {actual_label}", palette[self.rep_mode], size=20, zorder=self.rep_mode)
            self._scatter_vias(predicted_vias, f"Predicted Label: {predicted_label}", palette[1-self.rep_mode], size=20, zorder=(1-self.rep_mode))
            self._scatter_vias(cell_vias, "Cell vias", colors["cell_color"])

            self.ax.text(0.5, 0.9, f"Showing {self.i + 1}/{len(self.predicted_cells)} possible trojans", transform=self.ax.transAxes, ha='center')
            self.ax.text(0.5, 0.8, f"Distance to label: {self.distance(cell_vias, actual_vias if self.rep_mode else predicted_vias)}", transform=self.ax.transAxes, ha='center')
        else:
            self.ax.text(
                0.5, 0.5, 
                "All trojans have been checked",
                ha='center', va='center',
                fontsize=16,
                color='red',
                transform=self.ax.transAxes  # coordinates relative to axes (0–1)
            )
            self.ax.set_xticks([])
            self.ax.set_yticks([])
        plt.draw()


    def _make_box(self, current: Predicted_Cell) -> None:
        """
        Draw a rectangle representing the cell's bounding box.

        Args:
            current (Predicted_Cell): The cell whose box will be drawn.
        """
        top_left, bottom_right = current["cell"]["box"]
        width = bottom_right[0] - top_left[0]
        height = bottom_right[1] - top_left[1]
        box = Rectangle((0, 0), width, height, facecolor='none', edgecolor=colors['line_color'])
        self.ax.add_patch(box)
        self.ax.set_xlim(-self.zoom*width, (1+self.zoom)*width)
        self.ax.set_ylim(-self.zoom*height, (1+self.zoom)*height)


    def _scatter_vias(self, vias: list[Point], description: str, color: str, size: int = 10, zorder: int = 2) -> None:
        """
        Plot vias as a scatter plot.

        Args:
            vias (list[Point]): List of points to plot.
            description (str): Label for the scatter points.
            color (str): Color for the points.
            size (int): Marker size.
            zorder (int): Layer order for plotting.
        """
        x_val, y_val = zip(*vias) if vias else ([], [])
        sns.scatterplot(x=x_val, y=y_val, s=size, color=color, edgecolor='none', label=description, zorder=zorder, ax=self.ax)


    def _on_key(self, event: Event) -> None:
        """
        Handle keyboard input to navigate cells, toggle representative mode, or confirm trojans.
        Expects a KeyEvent

        Args:
            event (Event): The key event triggered by matplotlib.
        """
        if self._get_current() is None or not hasattr(event, 'key'): return
        else: event = cast(KeyEvent, event)

        allowed = ["left","right", "up", "down", "p", "m", "a", "d"]
        if event.key not in allowed: return
        elif event.key == "left" or event.key == "right":
            self.i += 1 - 2*(event.key == "left")
        elif event.key == "up" or event.key == "down":
            self.rep_mode = not self.rep_mode
        elif event.key == "p" or event.key == "m":
            self.zoom += 0.1 - 0.2*(event.key == "p")
        elif event.key == "a" or event.key == "d":
            self._confirm_trojan(event.key == "a")
        
        self._show_current()


    def start_interactive(self) -> None:
        """
        Start the interactive inspection session.

        Connects keyboard events to the inspector.
        Keyboard controls:
        - left / right: navigate between cells
        - up / down: toggle representative mode
        - p / m: zoom in (plus) / zoom out (minus)
        - a: confirm cell as trojan
        - d: discard cell
        """
        self.fig.canvas.mpl_connect('key_press_event', self._on_key)
        plt.show()
