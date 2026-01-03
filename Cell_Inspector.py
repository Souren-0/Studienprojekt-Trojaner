"""
Module: Cell_Inspector
Author: Souren Ishkhanian

This module helps manually checking for trojans.
"""

from typing import Any, Callable
from matplotlib.patches import Rectangle
from Annotation_Helpers import *
import matplotlib.pyplot as plt
import seaborn as sns
from Cell_Via_Utilities import align_vias

plt.style.use('dark_background')

colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
n = len(colors)
colors = {
    'line_color' : colors[1 % n],
    'cell_color' : colors[0 % n],
    'representative_color' : "red"
}

class Cell_Inspector:
    def __init__(self,
                 predicted_cells: list[Predicted_Cell],
                 trojans: list[Predicted_Cell],
                 representatives: dict[str, list[Point]],
                 distance_measure: Callable[[list[Point], list[Point]], float]) -> None:
        
        self.predicted_cells = predicted_cells
        self.representatives = representatives
        self.rep_mode = True
        self.i = 0

        self.distance = distance_measure
        self.trojans = trojans

        self.fig, self.ax = plt.subplots(figsize=(12,9))
        self._show_current()


    def _get_current(self):
        n = len(self.predicted_cells)
        self.i = (self.i % n) if n > 0 else 0
        return self.predicted_cells[self.i] if self.predicted_cells else None

    
    def _confirm_trojan(self, confirmed: bool):
        current = self._get_current()
        if confirmed and current:
            self.trojans.append(current)
        del self.predicted_cells[self.i]


    def _show_current(self):
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


    def _make_box(self, current):
        top_left, bottom_right = current["cell"]["box"]
        width = bottom_right[0] - top_left[0]
        height = bottom_right[1] - top_left[1]
        box = Rectangle((0, 0), width, height, facecolor='none', edgecolor=colors['line_color'])
        self.ax.add_patch(box)
        self.ax.set_xlim(-0.25*width, 1.25*width)
        self.ax.set_ylim(-0.25*height, 1.25*height)


    def _scatter_vias(self, vias: list[Point], description: str, color: str, size: int = 10, zorder: int = 2):
        x_val, y_val = zip(*vias) if vias else (None, None)
        sns.scatterplot(x=x_val, y=y_val, s=size, color=color, edgecolor='none', label=description, zorder=zorder, ax=self.ax)


    def _on_key(self, event):
        if self._get_current() is None: return

        if event.key == "right":
            self.i += 1
        elif event.key == "left":
            self.i -= 1
        elif event.key == "up" or event.key == "down":
            self.rep_mode = not self.rep_mode
        elif event.key == "a" or event.key == "d":
            self._confirm_trojan(event.key == "a")
        
        self._show_current()


    def start_interactive(self):
        self.fig.canvas.mpl_connect('key_press_event', self._on_key)
        plt.show()
