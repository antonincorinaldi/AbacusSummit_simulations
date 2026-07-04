import treecorr
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



# Compute the complex ellipticity e = |e| * exp(2i*r) from semi-axes a, b (a >= b)
# and position angle r [rad]
def e_complex(a,b,r):
    abs_e = (1-(b/a)) / (1+(b/a))
    e1 = abs_e*np.cos(2*r)
    e2 = abs_e*np.sin(2*r)
    return e1, e2


# Package eigenvectors/eigenvalues of a shape tensor into an Astropy Table using
# the AbacusSummit column naming convention (eigenvectors/eigenvalues must be
# ordered from least to greatest: minor, intermediate, major axis)
def format_ellipsoid(eigenvectors, eigenvalues, position = np.asarray([0,0,0])):

    el = Table()
    el['sigman_eigenvecsMin_L2com'] = eigenvectors[0]
    el['sigman_eigenvecsMid_L2com'] = eigenvectors[1]
    el['sigman_eigenvecsMaj_L2com'] = eigenvectors[2]
    
    el['sigman_L2com'] = np.sqrt(eigenvalues)
    
    return el
    
    




# Draw galaxy B/C axis ratios (tau_B, tau_C) from an independent bivariate normal
# distribution and rescale the halo semi-axes to get the galaxy's 3D shape.
# The major axis (A) is kept equal to the halo's major axis (no scatter on tau_A).
# Draws are repeated (rejection sampling) per halo until B <= A and C <= A hold.
def population_3D(
                  mu_tau_B, mu_tau_C,
                  sigma_tau_B, sigma_tau_C,
                  el,
                  nb_halos,
                  ):

    halos_table2 = el.copy()
    axis_orig = np.array(halos_table2['sigman_L2com'])

    valid_axis = np.zeros((nb_halos, 3))

    pending = np.ones(nb_halos, dtype=bool) # Table of True = nb of halos without any central galaxy

    # Diagonal covariance: tau_B and tau_C are drawn independently (no cross-correlation)
    cov = [
           [sigma_tau_B**2, 0],
           [0             , sigma_tau_C**2]]


    while np.any(pending):

        n_pending = pending.sum() # Count the number of True that remain

        taus = np.random.multivariate_normal(
            mean=[mu_tau_B, mu_tau_C],
            cov=cov,
            size=n_pending
        )

        tau_B = np.clip(taus[:, 0], 0, 1.)
        tau_C = np.clip(taus[:, 1], 0, 1.)

        Ag = axis_orig[pending, 0]  # major axis unchanged (tau_A = 1)
        Bg = axis_orig[pending, 1] * tau_B
        Cg = axis_orig[pending, 2] * tau_C

        # Accept only draws respecting the physical ordering A >= B and A >= C
        mask = (
            (Bg <= Ag) &
            (Cg <= Ag)
        )

        accepted_idx = np.where(pending)[0][mask] # Index of the galaxies that respect the conditions of the mask

        valid_axis[accepted_idx] = np.stack(
            (Ag[mask], Bg[mask], Cg[mask]), axis=1 #
        )

        pending[accepted_idx] = False

    # Halo eigenvectors giving the 3D orientation of the shape tensor
    eigenvecs_Min = halos_table2['sigman_eigenvecsMin_L2com']
    eigenvecs_Mid = halos_table2['sigman_eigenvecsMid_L2com']
    eigenvecs_Max = halos_table2['sigman_eigenvecsMaj_L2com']

    eigenvectors = np.stack(
        (eigenvecs_Min, eigenvecs_Mid, eigenvecs_Max),
        axis=1
    )

    # Build one formatted ellipsoid table per galaxy from its accepted axes
    ellipses = np.array([
        format_ellipsoid(eigenvectors[i], valid_axis[i]**2)
        for i in range(nb_halos)
    ])

    # Re-order axes as (major, mid, minor) for the projection step
    evcl = np.array([
        ellipses['sigman_eigenvecsMaj_L2com'],
        ellipses['sigman_eigenvecsMid_L2com'],
        ellipses['sigman_eigenvecsMin_L2com']
    ])
    evcl = np.transpose(evcl, (1, 0, 2))  # shape: (nb_halos, 3, 3)

    evls = ellipses['sigman_L2com']**2  # eigenvalues (squared semi-axes)

    return evcl, evls

    




# Project a 3D ellipsoid (eigenvectors evcl, eigenvalues evls) onto the sky plane
# and return its complex ellipticity (e1, e2). `p_axis` sets the line-of-sight axis:
# 'x' or 'y'. Follows the projected-ellipse algebra of eq. 23 (Schneider & Bartelmann).
def projection (evcl, evls, p_axis=''):

    # Projection 3D => 2D
    if p_axis=='x': # Projection perpendicular to the LOS
        K = np.sum(evcl[:,:,0][:,:,None]*(evcl/evls[:,None]), axis=1)
        r = evcl[:,:,2] - evcl[:,:,0] * K[:,2][:,None] / K[:,0][:,None]
        s = evcl[:,:,1] - evcl[:,:,0] * K[:,1][:,None] / K[:,0][:,None]

    if p_axis=='y': # Projection along the LOS
        K = np.sum(evcl[:,:,1][:,:,None] * (evcl/evls[:,None]), axis=1)
        r = evcl[:,:,0] - evcl[:,:,1] * K[:,0][:,None] / K[:,1][:,None]
        s = evcl[:,:,2] - evcl[:,:,1] * K[:,2][:,None] / K[:,1][:,None]


    # Coefficients A,B,C (eq 23 of (2))
    A1 = np.sum(r**2 / evls, axis=1)
    B1 = np.sum(2*r*s / evls, axis=1)
    C1 = np.sum(s**2 / evls, axis=1)


    # Axis a_p,b_p and orientation angle r_p of the projected galaxy
    r_p = np.pi / 2 + np.arctan2(B1,A1-C1)/2
    a_p = 1/np.sqrt((A1+C1)/2 + (A1-C1)/(2*np.cos(2*r_p)))
    b_p = 1/np.sqrt(A1+C1-(1/a_p**2))


    # Projected ellipticity
    e1, e2 = e_complex(a_p, b_p, r_p)

    return e1, e2




# Top-level wrapper: draw the galaxy population (population_3D) and project it
# (projection) to return the ellipticity components (e1, e2) for nb_halos galaxies.
# theta = [mu_tau_B, mu_tau_C, sigma_tau_B, sigma_tau_C]; set the sigmas to 0
# for a deterministic (no-scatter) run.
def simulator(theta,
                el,
                nb_halos,
                p_axis='y'
               ):

    mu_tau_B, mu_tau_C, sigma_tau_B, sigma_tau_C, = theta

    evcl, evls = population_3D (mu_tau_B, mu_tau_C, sigma_tau_B, sigma_tau_C, el, nb_halos)

    e1, e2 = projection (evcl, evls, p_axis=p_axis)

    return e1, e2


