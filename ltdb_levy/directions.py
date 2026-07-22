"""Direction sampling and reproducible Haar-uniform three-dimensional rotations."""

from __future__ import annotations

import numpy as np


def isotropic_directions(size: int, rng: np.random.Generator) -> np.ndarray:
    """Draw unit vectors uniformly on the sphere.

    Uniform azimuth and uniform ``cos(phi)`` are required; drawing ``phi``
    uniformly would over-sample both poles (the bug fixed in ``func.c``).
    """
    if size < 0:
        raise ValueError("size must be non-negative")
    azimuth = rng.uniform(0.0, 2.0 * np.pi, size=size)
    z = rng.uniform(-1.0, 1.0, size=size)
    radial = np.sqrt(np.maximum(0.0, 1.0 - z * z))
    return np.column_stack((radial * np.cos(azimuth), radial * np.sin(azimuth), z))


def random_quaternions(size: int, rng: np.random.Generator) -> np.ndarray:
    """Draw Haar-uniform unit quaternions using Shoemake's construction.

    The returned component order is ``(x, y, z, w)``.
    """
    if size < 0:
        raise ValueError("size must be non-negative")
    u1 = rng.random(size)
    u2 = rng.random(size)
    u3 = rng.random(size)
    a = np.sqrt(1.0 - u1)
    b = np.sqrt(u1)
    return np.column_stack(
        (
            a * np.sin(2.0 * np.pi * u2),
            a * np.cos(2.0 * np.pi * u2),
            b * np.sin(2.0 * np.pi * u3),
            b * np.cos(2.0 * np.pi * u3),
        )
    )


def rotate_vectors(vectors: np.ndarray, quaternions: np.ndarray) -> np.ndarray:
    """Rotate vectors by unit quaternions without constructing rotation matrices.

    Either one quaternion may be supplied for a whole block, or one quaternion
    per vector.
    """
    vectors = np.asarray(vectors, dtype=float)
    quaternions = np.asarray(quaternions, dtype=float)
    if vectors.ndim != 2 or vectors.shape[1] != 3:
        raise ValueError("vectors must have shape (n, 3)")
    if quaternions.ndim == 1:
        quaternions = quaternions.reshape(1, 4)
    if quaternions.ndim != 2 or quaternions.shape[1] != 4:
        raise ValueError("quaternions must have shape (1, 4) or (n, 4)")
    if len(quaternions) not in (1, len(vectors)):
        raise ValueError("Provide one quaternion or one quaternion per vector")

    q_xyz = quaternions[:, :3]
    q_w = quaternions[:, 3:4]
    if len(quaternions) == 1 and len(vectors) != 1:
        q_xyz = np.broadcast_to(q_xyz, vectors.shape)
        q_w = np.broadcast_to(q_w, (len(vectors), 1))
    t = 2.0 * np.cross(q_xyz, vectors)
    return vectors + q_w * t + np.cross(q_xyz, t)


def independently_rotate_vectors(
    vectors: np.ndarray, rng: np.random.Generator
) -> np.ndarray:
    """Apply an independent random rotation to each vector."""
    vectors = np.asarray(vectors, dtype=float)
    return rotate_vectors(vectors, random_quaternions(len(vectors), rng))

