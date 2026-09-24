"""Spinniped: a finite-element toolbox for rotor and shaft dynamics."""

__version__ = "0.1.0"

from .builder import BuiltModel, ModelBuilder, build_model
from .records import (
    BearingElement, BearingProperty, CoordinateSystem, DiskElement,
    DiskProperty, Grid, LumpedMassElement, LumpedMassProperty, Material,
    ModelDefinition, RandomDistribution, ShaftElement, ShaftProperty,
)
from .solver import Solver
from .plotting import plot_campbell

__all__ = [
    "BearingElement", "BearingProperty", "BuiltModel", "CoordinateSystem",
    "DiskElement", "DiskProperty", "Grid", "LumpedMassElement",
    "LumpedMassProperty", "Material", "ModelBuilder", "ModelDefinition",
    "RandomDistribution",
    "ShaftElement", "ShaftProperty", "Solver", "build_model",
    "plot_campbell",
]
