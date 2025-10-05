from typing import TypedDict, TypeAlias
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