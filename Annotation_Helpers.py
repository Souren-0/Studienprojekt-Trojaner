from typing import TypedDict, TypeAlias, Any, Callable, Iterable, Optional
from pathlib import Path


Point: TypeAlias = tuple[int | float, int | float]


class Cell_Data(TypedDict):
    name: str
    rotation: int | float
    magnification: int | float
    reflection: bool


class Cell(TypedDict):
    data: Cell_Data
    box: tuple[Point, Point]
    image: str | Path
    vias: list[Point]


class Predicted_Cell(TypedDict):
    cell: Cell
    actual: str
    predicted: str