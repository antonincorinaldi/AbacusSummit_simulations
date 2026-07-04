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
# Hardcoded path to the functions/ directory; update if the directory layout changes
sys.path.append(os.path.abspath("/n17data/corinaldi/halos/Abacus_HOD/abacusutils/functions"))
from functions import *




z_snapshot= 0.2
tracer = 'LRG'


# -----------------------------------------------------------------------
# 1. Build associated halo catalog for LRG galaxies
# -----------------------------------------------------------------------
# Open the mock galaxy catalogue (from generate_galaxy_mock.py) and read the
# extension corresponding to the chosen tracer
with fits.open(f"/n17data/corinaldi/mock_galaxy_catalogues/mock_galaxy_catalogue_LRG+ELG_z{z_snapshot}.fits") as hdul:
    mock = hdul[tracer].data

mock = Table(mock)
mock['WEIGHT_TOT'] = np.ones(len(mock))

# Halo IDs associated with the mock galaxies
halos_ids_mock = np.asarray(mock['id'])

cat = CompaSOHaloCatalog(
    f'/n17data/corinaldi/halos/AbacusSummit_base_c000_ph000/halos/z{z_snapshot}',
    fields=[
        'x_L2com',
        'sigman_L2com',
        'sigman_eigenvecsMin_L2com',
        'sigman_eigenvecsMid_L2com',
        'sigman_eigenvecsMaj_L2com',
        'N',
        'id',
    ],
    cleaned=False,
)

# Match halos to mock galaxies by ID
halos_ids_mock = halos_ids_mock.astype(cat.halos['id'].dtype)
mask = np.isin(cat.halos['id'], halos_ids_mock)
halos = Table(cat.halos[mask])
halos.meta.clear()


# Selecting only 2_000_000 halos
# Cap the sample size for memory/runtime reasons: subsample randomly without
# replacement if there are more matched halos than nb_halos_max
nb_halos_max = 2_000_000
if len(halos) > nb_halos_max:
    rng_idx = np.random.choice(len(halos), size=nb_halos_max, replace=False)
    halos = halos[rng_idx]



x = halos['x_L2com'][:, 0]
y = halos['x_L2com'][:, 1]
z = halos['x_L2com'][:, 2]



# Compute ellipticity components (e1, e2) with no scatter on the axis ratios
# (theta = [mu_tau_B, mu_tau_C, sigma_tau_B, sigma_tau_C] = [1, 1, 0, 0]),
# i.e. galaxy shapes are set equal to their host halo's shape
e1_halos, e2_halos = simulator([1, 1, 0, 0], el=halos, nb_halos=len(halos))

halos['x'] = x                        # RA-like coordinate [Mpc/h]
halos['y'] = y                        # Dec-like coordinate [Mpc/h]
halos['z'] = z                        # line-of-sight coordinate [Mpc/h]
halos['e1'] = e1_halos                # first ellipticity component
halos['e2'] = e2_halos                # second ellipticity component
halos['weights'] = np.ones(len(halos))  # uniform halo weights



# Write the enriched halo catalogue to disk
halos.write(f'/n17data/corinaldi/halos/Abacus_HOD_halos/halo_catalogue_{tracer}_z{z_snapshot}.fits', format='fits', overwrite=True)










