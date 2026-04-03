"""
functions.py
------------
Shared utility functions for computing projected halo ellipticities from
AbacusSummit shape-tensor data.

The pipeline implemented here follows a three-step approach:
    1. `population_3D`  — draw galaxy axis ratios from a multivariate normal
                          distribution and scale them by the halo semi-axes.
    2. `projection`     — project the 3D ellipsoid onto the sky plane and
                          compute the projected semi-axes and position angle.
    3. `simulator`      — top-level wrapper that chains steps 1 and 2 and
                          returns the complex ellipticity components (e1, e2).

Helper functions:
    `e_complex`         — converts semi-axes and position angle to (e1, e2).
    `format_ellipsoid`  — packages eigenvectors/eigenvalues into an Astropy
                          Table with the AbacusSummit column naming convention.

References
----------
(1) Zheng et al. 2007, ApJ, 667, 760  — HOD model
(2) Schneider & Bartelmann 1995 (or equivalent) — projected ellipsoid algebra
"""

import treecorr
import IACorr
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import time as time
import pandas as pd
from scipy import stats
import math
from astropy.table import vstack, Table
from astropy import units as u


# ---------------------------------------------------------------------------
# Ellipticity helpers
# ---------------------------------------------------------------------------

def e_complex(a, b, r):
    """
    Convert ellipse parameters to complex ellipticity components.

    Parameters
    ----------
    a : array_like
        Semi-major axis length.
    b : array_like
        Semi-minor axis length (b <= a).
    r : array_like
        Position angle of the major axis [radians].

    Returns
    -------
    e1, e2 : ndarray
        Real and imaginary parts of the complex ellipticity
        e = |e| * exp(2i * r), where |e| = (a - b) / (a + b).
    """
    abs_e = (1 - (b / a)) / (1 + (b / a))  # ellipticity modulus in [0, 1)
    e1 = abs_e * np.cos(2 * r)
    e2 = abs_e * np.sin(2 * r)
    return e1, e2


def format_ellipsoid(eigenvectors, eigenvalues, position=np.asarray([0, 0, 0])):
    """
    Package 3D ellipsoid parameters into an Astropy Table using the AbacusSummit
    column naming convention.

    Eigenvectors and eigenvalues must be ordered from least to greatest
    (minor → intermediate → major axis), matching the AbacusSummit convention.

    Parameters
    ----------
    eigenvectors : array_like, shape (3, 3)
        Rows are the eigenvectors of the shape tensor, ordered
        [minor, intermediate, major].
    eigenvalues : array_like, shape (3,)
        Eigenvalues of the shape tensor (variances), same ordering.
    position : array_like, shape (3,), optional
        Centre-of-mass position (not used internally; kept for API compatibility).

    Returns
    -------
    el : astropy.table.Table
        Table with columns:
        - ``sigman_eigenvecsMin_L2com`` : minor-axis eigenvector
        - ``sigman_eigenvecsMid_L2com`` : intermediate-axis eigenvector
        - ``sigman_eigenvecsMaj_L2com`` : major-axis eigenvector
        - ``sigman_L2com``              : sqrt(eigenvalues) = semi-axes [Mpc/h]
    """
    el = Table()
    el['sigman_eigenvecsMin_L2com'] = eigenvectors[0]  # minor-axis unit vector
    el['sigman_eigenvecsMid_L2com'] = eigenvectors[1]  # intermediate-axis unit vector
    el['sigman_eigenvecsMaj_L2com'] = eigenvectors[2]  # major-axis unit vector
    el['sigman_L2com'] = np.sqrt(eigenvalues)          # semi-axis lengths (std dev)
    return el


# ---------------------------------------------------------------------------
# 3D galaxy population
# ---------------------------------------------------------------------------

def population_3D(mu_tau_A, mu_tau_B, mu_tau_C,
                  sigma_tau_A, sigma_tau_B, sigma_tau_C,
                  r_tau,
                  el,
                  nb_halos):
    """
    Generate 3D galaxy axis ratios by drawing from a correlated multivariate
    normal distribution and rescaling by the halo semi-axes.

    Each galaxy semi-axis is obtained as:
        A_gal = axis_halo[0] * tau_A
        B_gal = axis_halo[1] * tau_B * tau_A
        C_gal = axis_halo[2] * tau_C * tau_A

    Draws are repeated (rejection sampling) until the ordering constraints
    B_gal >= C_gal > 0 and A_gal >= B_gal are satisfied for every halo.

    Parameters
    ----------
    mu_tau_A, mu_tau_B, mu_tau_C : float
        Mean axis-ratio parameters for the A, B, C semi-axes.
    sigma_tau_A, sigma_tau_B, sigma_tau_C : float
        Standard deviations of the axis-ratio distributions.
    r_tau : float
        Correlation coefficient shared across all off-diagonal covariance terms.
    el : astropy.table.Table
        Halo table containing ``sigman_L2com`` and the three eigenvector columns
        (as produced by CompaSOHaloCatalog).
    nb_halos : int
        Number of halos to process.

    Returns
    -------
    evcl : ndarray, shape (nb_halos, 3, 3)
        Eigenvector matrices for each galaxy (rows: major, mid, minor axes).
    evls : ndarray, shape (nb_halos, 3)
        Squared semi-axis lengths (eigenvalues) for each galaxy.
    """
    halos_table2 = el.copy()
    axis_orig = np.array(halos_table2['sigman_L2com'])  # halo semi-axes [Mpc/h]

    valid_axis = np.zeros((nb_halos, 3))  # container for accepted galaxy semi-axes

    # Track which halos still need a valid draw
    pending = np.ones(nb_halos, dtype=bool)

    # Full 3x3 covariance matrix for the correlated axis-ratio draw
    cov = [
        [sigma_tau_A**2,                     r_tau * sigma_tau_A * sigma_tau_B, r_tau * sigma_tau_A * sigma_tau_C],
        [r_tau * sigma_tau_B * sigma_tau_A,  sigma_tau_B**2,                    r_tau * sigma_tau_B * sigma_tau_C],
        [r_tau * sigma_tau_A * sigma_tau_C,  r_tau * sigma_tau_B * sigma_tau_C, sigma_tau_C**2],
    ]

    # Rejection-sampling loop: keep drawing until all halos have valid axes
    while np.any(pending):
        n_pending = pending.sum()

        # Draw correlated tau values for all pending halos simultaneously
        taus = np.random.multivariate_normal(
            mean=[mu_tau_A, mu_tau_B, mu_tau_C],
            cov=cov,
            size=n_pending,
        )

        # Clip axis ratios to [0, 1] to enforce physical bounds
        tau_A = np.clip(taus[:, 0], 0, 1)
        tau_B = np.clip(taus[:, 1], 0, 1)
        tau_C = np.clip(taus[:, 2], 0, 1)

        # Scale halo semi-axes by the drawn axis ratios
        Ag = axis_orig[pending, 0] * tau_A
        Bg = axis_orig[pending, 1] * tau_B * tau_A
        Cg = axis_orig[pending, 2] * tau_C * tau_A

        # Accept only draws that respect the ordering A >= B >= C > 0
        mask = (
            (Bg / Ag <= 1) &
            (Cg / Ag <= 1) &
            (Bg >= Cg) &
            (Bg > 0) &
            (Cg > 0)
        )

        # Store accepted semi-axes at the correct global indices
        accepted_idx = np.where(pending)[0][mask]
        valid_axis[accepted_idx] = np.stack((Ag[mask], Bg[mask], Cg[mask]), axis=1)
        pending[accepted_idx] = False

    # Retrieve halo eigenvectors (shape orientation in 3D)
    eigenvecs_Min = halos_table2['sigman_eigenvecsMin_L2com']
    eigenvecs_Mid = halos_table2['sigman_eigenvecsMid_L2com']
    eigenvecs_Max = halos_table2['sigman_eigenvecsMaj_L2com']

    # Stack eigenvectors into a (nb_halos, 3, 3) array
    eigenvectors = np.stack((eigenvecs_Min, eigenvecs_Mid, eigenvecs_Max), axis=1)

    # Build formatted ellipsoid tables for each galaxy
    ellipses = np.array([
        format_ellipsoid(eigenvectors[i], valid_axis[i]**2)
        for i in range(nb_halos)
    ])

    # Re-order columns: (major, mid, minor) for the projection step
    evcl = np.array([
        ellipses['sigman_eigenvecsMaj_L2com'],
        ellipses['sigman_eigenvecsMid_L2com'],
        ellipses['sigman_eigenvecsMin_L2com'],
    ])
    evcl = np.transpose(evcl, (1, 0, 2))  # shape: (nb_halos, 3, 3)

    evls = ellipses['sigman_L2com']**2  # eigenvalues (squared semi-axes)

    return evcl, evls


# ---------------------------------------------------------------------------
# Projection: 3D ellipsoid -> 2D projected ellipse
# ---------------------------------------------------------------------------

def projection(evcl, evls, p_axis=''):
    """
    Project a 3D ellipsoid onto a 2D plane and compute projected ellipticity.

    The projection algebra follows eq. 23 of the reference (Schneider &
    Bartelmann or equivalent). The line-of-sight direction is chosen via
    `p_axis`.

    Parameters
    ----------
    evcl : ndarray, shape (N, 3, 3)
        Eigenvector matrices for N galaxies (major, mid, minor axes stacked).
    evls : ndarray, shape (N, 3)
        Eigenvalues (squared semi-axes) for N galaxies.
    p_axis : {'x', 'y'}
        Axis along which to project:
        - ``'x'``: project perpendicular to the x-axis (LOS along x).
        - ``'y'``: project perpendicular to the y-axis (LOS along y).

    Returns
    -------
    e1, e2 : ndarray, shape (N,)
        Projected complex ellipticity components.
    """
    # --- Project the shape tensor onto the 2D plane ---
    if p_axis == 'x':  # line-of-sight is the x-axis
        K = np.sum(evcl[:, :, 0][:, :, None] * (evcl / evls[:, None]), axis=1)
        r = evcl[:, :, 2] - evcl[:, :, 0] * K[:, 2][:, None] / K[:, 0][:, None]
        s = evcl[:, :, 1] - evcl[:, :, 0] * K[:, 1][:, None] / K[:, 0][:, None]

    if p_axis == 'y':  # line-of-sight is the y-axis
        K = np.sum(evcl[:, :, 1][:, :, None] * (evcl / evls[:, None]), axis=1)
        r = evcl[:, :, 0] - evcl[:, :, 1] * K[:, 0][:, None] / K[:, 1][:, None]
        s = evcl[:, :, 2] - evcl[:, :, 1] * K[:, 2][:, None] / K[:, 1][:, None]

    # --- Coefficients of the projected ellipse equation (ref. eq. 23) ---
    A1 = np.sum(r**2 / evls, axis=1)
    B1 = np.sum(2 * r * s / evls, axis=1)
    C1 = np.sum(s**2 / evls, axis=1)

    # --- Projected semi-axes and position angle ---
    r_p = np.pi / 2 + np.arctan2(B1, A1 - C1) / 2          # position angle [rad]
    a_p = 1 / np.sqrt((A1 + C1) / 2 + (A1 - C1) / (2 * np.cos(2 * r_p)))  # semi-major
    b_p = 1 / np.sqrt(A1 + C1 - (1 / a_p**2))               # semi-minor

    # --- Convert to complex ellipticity ---
    e1, e2 = e_complex(a_p, b_p, r_p)

    return e1, e2


# ---------------------------------------------------------------------------
# Top-level simulator
# ---------------------------------------------------------------------------

def simulator(theta, el, nb_halos, p_axis='y'):
    """
    Generate projected ellipticity components for a population of halos.

    This is the main entry point used by generate_halo_catalogue.py.
    It chains `population_3D` and `projection` together.

    Parameters
    ----------
    theta : array_like, length 7
        Model parameters in order:
        [mu_tau_A, mu_tau_B, mu_tau_C,
         sigma_tau_A, sigma_tau_B, sigma_tau_C, r_tau].
        Set all sigmas and r_tau to 0 for a deterministic (no-scatter) run.
    el : astropy.table.Table
        Halo table with AbacusSummit shape-tensor columns (see `population_3D`).
    nb_halos : int
        Number of halos to process.
    p_axis : {'x', 'y'}, optional
        Projection axis (default ``'y'`` = line-of-sight along y).

    Returns
    -------
    e1, e2 : ndarray, shape (nb_halos,)
        Projected ellipticity components.
    """
    mu_tau_A, mu_tau_B, mu_tau_C, sigma_tau_A, sigma_tau_B, sigma_tau_C, r_tau = theta

    # Step 1: generate galaxy 3D axis ratios scaled by halo semi-axes
    evcl, evls = population_3D(
        mu_tau_A, mu_tau_B, mu_tau_C,
        sigma_tau_A, sigma_tau_B, sigma_tau_C,
        r_tau, el, nb_halos,
    )

    # Step 2: project onto the sky plane and compute ellipticity components
    e1, e2 = projection(evcl, evls, p_axis=p_axis)

    return e1, e2
