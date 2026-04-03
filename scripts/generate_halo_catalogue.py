"""
generate_halo_catalogue.py
--------------------------
Builds a halo catalogue matched to a mock galaxy catalogue output by
generate_galaxy_mock.py, then computes halo ellipticities and writes the
result to a FITS file for downstream analysis.

Pipeline overview:
    [1] prepare_sim.py            <- particle subsample preparation (run first)
    [2] generate_galaxy_mock.py   <- HOD galaxy catalog generation
    [3] generate_halo_catalogue.py <- THIS SCRIPT (run after step 2)

Steps performed here:
    1. Load the mock galaxy catalogue produced in step 2.
    2. Load the full AbacusSummit CompaSOHaloCatalog at the same redshift.
    3. Match halos to mock galaxies by halo ID.
    4. Compute ellipticity components (e1, e2) for the matched halos.
    5. Write the enriched halo catalogue to a FITS file.

Usage:
    From the project root:
        python scripts/generate_halo_catalogue.py

    Update the hard-coded file paths if the directory layout changes.
"""

from abacusnbody.hod import prepare_sim
from astropy.io import fits
from astropy.table import Table
from astropy import units as u
from abacusnbody.data.compaso_halo_catalog import CompaSOHaloCatalog
from astropy.cosmology import FlatLambdaCDM
import matplotlib.pyplot as plt
import time
import yaml
import numpy as np
from abacusnbody.hod.abacus_hod import AbacusHOD

import sys
import os

# Add the local functions/ directory (at the project root) to the Python path
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'functions'))
from functions import *


# Redshift of the snapshot to process
z_snapshot = 0.5

# Galaxy tracer type to extract from the mock catalogue
tracer = 'LRG'


# -----------------------------------------------------------------------
# 1. Load mock galaxy catalogue and extract halo IDs
# -----------------------------------------------------------------------
# Open the FITS file produced by generate_galaxy_mock.py and read the
# extension corresponding to the chosen tracer (LRG or ELG).
with fits.open(f"/n17data/corinaldi/mock_galaxy_catalogues/mock_galaxy_catalogue_LRG+ELG_z{z_snapshot}.fits") as hdul:
    mock = hdul[tracer].data

mock = Table(mock)

# Add a uniform weight column (all galaxies weighted equally)
mock['WEIGHT_TOT'] = np.ones(len(mock))

# Extract the halo IDs associated with mock galaxies
halos_ids_mock = np.asarray(mock['id'])


# -----------------------------------------------------------------------
# 2. Load the full AbacusSummit halo catalogue at the target redshift
# -----------------------------------------------------------------------
# Only the fields needed for position, shape, particle count, and ID are
# requested to limit memory usage.
cat = CompaSOHaloCatalog(
    f'/n17data/corinaldi/halos/AbacusSummit_base_c000_ph000/halos/z{z_snapshot}',
    fields=[
        'x_L2com',              # L2-center-of-mass position [Mpc/h]
        'sigman_L2com',         # velocity dispersion tensor (scalar, L2com frame)
        'sigman_eigenvecsMin_L2com',  # eigenvector along minor axis
        'sigman_eigenvecsMid_L2com',  # eigenvector along intermediate axis
        'sigman_eigenvecsMaj_L2com',  # eigenvector along major axis
        'N',                    # number of particles in the halo
        'id',                   # unique halo identifier
    ],
    cleaned=False,
)


# -----------------------------------------------------------------------
# 3. Match halos to mock galaxies by halo ID
# -----------------------------------------------------------------------
# Cast mock halo IDs to the same dtype as the catalogue IDs before masking
halos_ids_mock = halos_ids_mock.astype(cat.halos['id'].dtype)

# Boolean mask: True for halos whose ID appears in the mock galaxy list
mask = np.isin(cat.halos['id'], halos_ids_mock)

# Retain only the matched halos and clear simulation metadata from the header
halos = Table(cat.halos[mask])
halos.meta.clear()


# -----------------------------------------------------------------------
# 4. Compute halo ellipticities
# -----------------------------------------------------------------------
# Extract Cartesian position components for convenience
x = halos['x_L2com'][:, 0]
y = halos['x_L2com'][:, 1]
z = halos['x_L2com'][:, 2]

# Compute ellipticity components (e1, e2) using the shape tensor eigenvectors.
# The first argument is a weight vector; here all components are set to 1
# (isotropic weighting). `el` is the halo table; `nb_halos` sets the sample size.
e1_halos, e2_halos = simulator([1, 1, 1, 0, 0, 0, 0], el=halos, nb_halos=len(halos))

# Add derived columns to the halo table
halos['x'] = x                        # RA-like coordinate [Mpc/h]
halos['y'] = y                        # Dec-like coordinate [Mpc/h]
halos['z'] = z                        # line-of-sight coordinate [Mpc/h]
halos['e1'] = e1_halos                # first ellipticity component
halos['e2'] = e2_halos                # second ellipticity component
halos['weights'] = np.ones(len(halos))  # uniform halo weights


# -----------------------------------------------------------------------
# 5. Write the matched halo catalogue to disk
# -----------------------------------------------------------------------
halos.write(
    f'/n17data/corinaldi/halos/Abacus_HOD_halos/halo_catalogue_{tracer}_z{z_snapshot}.fits',
    format='fits',
    overwrite=True,
)
