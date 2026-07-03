# Compute the projected galaxy clustering signal wp(rp) for a mock galaxy
# catalogue populated in a periodic AbacusSummit box (Step 2 output).
# Because the box is periodic, no random catalogue is needed: pycorr uses
# the box periodicity directly (`boxsize` argument) instead of a data/random
# estimator.

from abacusnbody.hod import prepare_sim
from astropy.io import fits
from astropy.table import Table
from astropy import units as u
from abacusnbody.data.compaso_halo_catalog import CompaSOHaloCatalog
from astropy.cosmology import FlatLambdaCDM
import matplotlib.pyplot as plt
import os
import glob
import time
import yaml
import numpy as np
import argparse
from abacusnbody.hod.abacus_hod import AbacusHOD
import treecorr
from tqdm import tqdm
import torch
import pycorr
from pycorr import TwoPointCorrelationFunction




tracer = 'LRG' ; hdu = 1 # hdu=1 for LRGs, hdu=2 for ELGs
z_snapshot = 0.2


# Load the mock galaxy catalogue produced by generate_galaxy_mock.py and
# randomly subsample it down to nb_halos objects to keep the correlation
# function computation tractable.
halo_catalogue = Table.read(f"/n17data/corinaldi/mock_galaxy_catalogues/mock_galaxy_catalogue_LRG+ELG_z{z_snapshot}.fits", hdu=hdu)
nb_halos = 2_000_000; halo_catalogue = halo_catalogue[np.random.choice(len(halo_catalogue), nb_halos, replace=False)]
x = halo_catalogue['x'].value
y = halo_catalogue['y'].value
z = halo_catalogue['z'].value
pos = np.array([x, y, z])


# rp/pi binning: nrp log-spaced radial bins (transverse separation) and npi
# linear bins along the line of sight, integrated up to pi_max to obtain wp.
nrp  = 20
rmin   = 1.
rmax   = 200.
pi_max = 60.
npi    = 10


edges = (np.geomspace(rmin, rmax, nrp), np.linspace(-pi_max, pi_max, npi))

# 'z' is used as the line-of-sight direction since the box is a simulation
# snapshot with no observer position; boxsize enables periodic wrapping.
result = TwoPointCorrelationFunction('rppi', edges, data_positions1=pos, boxsize=2000, engine='corrfunc', nthreads=1, los='z')

# Integrate xi(rp, pi) over pi to get the projected correlation function wp(rp).
rp, wgg = result.get_corr(return_sep=True, mode='wp')

# Save the result for later comparison against the data measurement (run_wgg_data.py).
np.savez(
    f'/n17data/corinaldi/halos/Abacus_HOD/pycorr/results_mock/{tracer}_mock_z{z_snapshot}.npz',
    rp = rp,
    wgg= wgg
)






