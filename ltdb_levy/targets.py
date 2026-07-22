"""C-compatible target dimensions, torus geometry, and membership tests."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Tuple

import numpy as np

from .errors import GeometryError

SUPPORTED_SHAPES: Tuple[str, ...] = ("Ball", "Disk", "Line")


def minimum_image(delta: np.ndarray, side_length: float) -> np.ndarray:
    """Map signed coordinate differences to the torus minimum image."""
    values = np.asarray(delta, dtype=float)
    if side_length <= 0.0:
        raise ValueError("side_length must be positive")
    return np.mod(values + side_length / 2.0, side_length) - side_length / 2.0


def wrap_positions(positions: np.ndarray, side_length: float) -> np.ndarray:
    """Wrap coordinates into the half-open box ``[0, side_length)``."""
    if side_length <= 0.0:
        raise ValueError("side_length must be positive")
    return np.mod(np.asarray(positions, dtype=float), side_length)


def dimension_from_projected_area(
    area: float, shape: str, detection_radius: float = 1.0
) -> float:
    """Invert the C simulator's face-on maximum-projection formulas.

    ``dimension`` is a diameter for Ball/Disk and the centre-line length for
    Line.  At ``d=1`` these reduce exactly to
    ``2*sqrt(A/pi)-2`` and ``(A-pi)/2`` from ``func.c``.
    """
    if shape not in SUPPORTED_SHAPES:
        raise GeometryError("Unsupported target shape {!r}".format(shape))
    if detection_radius <= 0.0:
        raise GeometryError("detection_radius must be positive")
    minimum = np.pi * detection_radius ** 2
    if not np.isfinite(area) or area <= minimum:
        raise GeometryError(
            "Projected area must exceed pi*d^2 ({:.6g}), got {}".format(minimum, area)
        )
    if shape in ("Ball", "Disk"):
        return float(2.0 * (np.sqrt(area / np.pi) - detection_radius))
    return float((area - np.pi * detection_radius ** 2) / (2.0 * detection_radius))


def projected_area(
    dimension: float, shape: str, detection_radius: float = 1.0
) -> float:
    """Return the face-on projected area of a detection neighbourhood."""
    if dimension < 0.0:
        raise GeometryError("target dimension must be non-negative")
    if shape in ("Ball", "Disk"):
        return float(np.pi * (dimension / 2.0 + detection_radius) ** 2)
    if shape == "Line":
        return float(
            2.0 * detection_radius * dimension + np.pi * detection_radius ** 2
        )
    raise GeometryError("Unsupported target shape {!r}".format(shape))


def dimension_from_neighbourhood_volume(
    volume: float, shape: str, detection_radius: float = 1.0
) -> float:
    """Invert the detection-neighbourhood volumes used by :class:`Target`.

    The nominal Disk and Line targets are lower-dimensional objects, so an
    intrinsic three-dimensional target volume is not common to all three
    shapes.  This function instead matches the volume within detection radius
    ``d`` of the target, which is exactly the region tested by ``contains``.
    """
    if shape not in SUPPORTED_SHAPES:
        raise GeometryError("Unsupported target shape {!r}".format(shape))
    if detection_radius <= 0.0:
        raise GeometryError("detection_radius must be positive")
    if not np.isfinite(volume) or volume <= 0.0:
        raise GeometryError(
            "Detection-neighbourhood volume must be positive, got {}".format(
                volume
            )
        )

    d = float(detection_radius)
    if shape == "Ball":
        radius = (3.0 * volume / (4.0 * np.pi)) ** (1.0 / 3.0)
        dimension = 2.0 * (radius - d)
    elif shape == "Disk":
        radius = np.sqrt(volume / (2.0 * d * np.pi))
        dimension = 2.0 * (radius - d)
    else:
        dimension = (
            volume - 4.0 * np.pi * d ** 3 / 3.0
        ) / (np.pi * d * d)

    tolerance = 1e-12 * max(1.0, d)
    if dimension < -tolerance:
        minimum = (
            2.0 * np.pi * d ** 3
            if shape == "Disk"
            else 4.0 * np.pi * d ** 3 / 3.0
        )
        raise GeometryError(
            "Volume {:.6g} is below the {:.6g} minimum for a non-negative "
            "{} target at detection radius {:.6g}".format(
                volume, minimum, shape, d
            )
        )
    return float(max(0.0, dimension))


@dataclass(frozen=True)
class Target:
    """One target centred in a periodic cubic domain."""

    shape: str
    dimension: float
    side_length: float
    detection_radius: float = 1.0
    center: Tuple[float, float, float] = None  # type: ignore[assignment]

    def __post_init__(self) -> None:
        if self.shape not in SUPPORTED_SHAPES:
            raise GeometryError("Unsupported target shape {!r}".format(self.shape))
        if self.dimension < 0.0:
            raise GeometryError("target dimension must be non-negative")
        if self.side_length <= 0.0 or self.detection_radius <= 0.0:
            raise GeometryError("side_length and detection_radius must be positive")
        if self.center is None:
            middle = self.side_length / 2.0
            object.__setattr__(self, "center", (middle, middle, middle))
        if len(self.center) != 3:
            raise GeometryError("target center must have three coordinates")

    def contains(self, positions: np.ndarray) -> np.ndarray:
        """Vectorized C-compatible membership, including boundary points."""
        points = np.asarray(positions, dtype=float)
        if points.shape[-1:] != (3,):
            raise ValueError("positions must end in a coordinate axis of length 3")
        delta = minimum_image(points - np.asarray(self.center), self.side_length)
        d = self.detection_radius
        if self.shape == "Ball":
            radius = self.dimension / 2.0 + d
            return np.sum(delta * delta, axis=-1) <= radius * radius
        if self.shape == "Disk":
            radius = self.dimension / 2.0 + d
            radial_squared = delta[..., 0] ** 2 + delta[..., 1] ** 2
            return (radial_squared <= radius * radius) & (np.abs(delta[..., 2]) <= d)
        # Capsule around a line segment oriented along y.
        residual_y = np.maximum(np.abs(delta[..., 1]) - self.dimension / 2.0, 0.0)
        distance_squared = (
            delta[..., 0] ** 2 + residual_y ** 2 + delta[..., 2] ** 2
        )
        return distance_squared <= d * d

    def neighbourhood_volume(self) -> float:
        """Closed-form Euclidean volume used by the analytic one-step check.

        The result assumes the neighbourhood does not overlap itself through
        periodic boundaries. Tests enforce that precondition.
        """
        d = self.detection_radius
        if self.shape == "Ball":
            radius = self.dimension / 2.0 + d
            return float(4.0 * np.pi * radius ** 3 / 3.0)
        if self.shape == "Disk":
            radius = self.dimension / 2.0 + d
            return float(2.0 * d * np.pi * radius ** 2)
        return float(np.pi * d * d * self.dimension + 4.0 * np.pi * d ** 3 / 3.0)
