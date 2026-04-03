"""
generate_galaxy_mock.py
-----------------------
Generates a mock galaxy catalogue from an AbacusSummit N-body simulation
using the Halo Occupation Distribution (HOD) framework provided by AbacusHOD.

This script must be run *after* prepare_sim.py has written the particle
subsample files for the target snapshot.

Pipeline overview:
    [1] prepare_sim.py            <- particle subsample preparation (run first)
    [2] generate_galaxy_mock.py   <- THIS SCRIPT
    [3] generate_halo_catalogue.py <- matched halo catalogue with ellipticities

Steps performed here:
    1. Load the YAML configuration and extract HOD / clustering parameters.
    2. Create an AbacusHOD object and run the HOD to populate halos with galaxies.
    3. Save the raw mock dictionary as a compressed NumPy archive (.npz).
    4. Reload the archive and convert each tracer (LRG, ELG) to a FITS extension.

Usage:
    From the project root:
        python scripts/generate_galaxy_mock.py

    Update `path2config` if the project is moved to a different machine.
"""

from abacusnbody.hod import prepare_sim
from astropy.io import fits
from astropy.table import Table
from astropy import units as u
from abacusnbody.data.compaso_halo_catalog import CompaSOHaloCatalog
from astropy.cosmology import FlatLambdaCDM
import matplotlib.pyplot as plt
import os
import time
import yaml
import numpy as np

from abacusnbody.hod.abacus_hod import AbacusHOD


# Cosmology matching AbacusSummit base cosmology (Planck 2018, h=1 units)
cosmo = FlatLambdaCDM(H0=100, Om0=0.3158)
h = 1  # dimensionless Hubble parameter (distances already in Mpc/h)

# Redshift of the snapshot to process — must match z_mock in the YAML
z_snapshot = 1.1

# Absolute path to the AbacusHOD configuration file
path2config = "/n17data/corinaldi/halos/Abacus_HOD/abacusutils/abacus_hod.yaml"

# Load the full configuration
config = yaml.safe_load(open(path2config))
sim_params = config['sim_params']
HOD_params = config['HOD_params']
clustering_params = config['clustering_params']

# --- HOD runtime flags ---
want_rsd = HOD_params['want_rsd']        # apply redshift-space distortions (RSD)?
write_to_disk = HOD_params['write_to_disk']  # let AbacusHOD write its own output?

# --- Projected clustering bin edges ---
bin_params = clustering_params['bin_params']

# rpbins: nbins+1 logarithmically-spaced edges => nbins bins in r_perp [Mpc/h]
rpbins = np.logspace(
    bin_params['logmin'], bin_params['logmax'], bin_params['nbins'] + 1
)
# Bin centres (arithmetic mean of edges) for labelling / plotting
rp_centers = 0.5 * (rpbins[1:] + rpbins[:-1])

# Line-of-sight integration limits for xi(rp, pi)
pimax = clustering_params['pimax']            # maximum pi [Mpc/h]
pi_bin_size = clustering_params['pi_bin_size']  # pi bin width [Mpc/h]


# -----------------------------------------------------------------------
# 1. Create AbacusHOD object and generate galaxy mock catalog
# -----------------------------------------------------------------------
# AbacusHOD reads the particle subsamples prepared by prepare_sim.py and
# places central and satellite galaxies into halos according to the HOD.
newBall = AbacusHOD(sim_params, HOD_params, clustering_params)

# Run the HOD: returns a nested dict {tracer_name: {field: array, ...}, ...}
# Nthread=1 keeps memory usage predictable; increase for faster runs on
# machines with many cores.
mock_dict = newBall.run_hod(
    newBall.tracers, want_rsd, write_to_disk=write_to_disk, Nthread=1
)

# Persist the raw mock dictionary to disk before converting to FITS,
# so intermediate results are not lost if the FITS step fails.
np.savez(
    f"/n17data/corinaldi/mock_galaxy_catalogues/mock_dict_LRG+ELG_z{z_snapshot}.npz",
    mock_dict=mock_dict,
    allow_pickle=True,
)


# -----------------------------------------------------------------------
# 2. Convert the mock dictionary to a multi-extension FITS file
# -----------------------------------------------------------------------
# Reload from disk to ensure the saved and served data are identical
mock = np.load(
    f"/n17data/corinaldi/mock_galaxy_catalogues/mock_dict_LRG+ELG_z{z_snapshot}.npz",
    allow_pickle=True,
)
mock = mock['mock_dict'].item()  # unwrap the object array produced by np.savez

# Start the HDU list with a minimal primary (no data, required by FITS standard)
hdul = [fits.PrimaryHDU()]

for tracer in ['LRG', 'ELG']:
    cat = mock[tracer]
    cols = []

    for key, val in cat.items():
        # Skip scalar metadata fields (e.g. Ncent — total central count)
        if np.isscalar(val):
            continue

        arr = np.asarray(val)

        # Choose FITS column format: 'K' (64-bit int) or 'D' (64-bit float)
        if np.issubdtype(arr.dtype, np.integer):
            fmt = 'K'
        else:
            fmt = 'D'

        cols.append(fits.Column(name=key, format=fmt, array=arr))

    # Name each extension after the tracer for easy access (e.g. hdul['LRG'])
    hdu = fits.BinTableHDU.from_columns(cols, name=tracer)

    # Store the number of central galaxies as a header keyword
    if 'Ncent' in cat:
        hdu.header['NCENT'] = cat['Ncent']

    hdul.append(hdu)


# Write the final multi-extension FITS catalogue to disk
fits.HDUList(hdul).writeto(
    f"/n17data/corinaldi/mock_galaxy_catalogues/mock_galaxy_catalogue_LRG+ELG_z{z_snapshot}.fits",
    overwrite=True,
)
