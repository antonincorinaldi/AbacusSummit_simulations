"""
generate_galaxy_mock.py
-----------------------
Main script for generating a galaxy mock catalog from AbacusSummit simulations
using the HOD (Halo Occupation Distribution) framework.

Pipeline:
    1. Load configuration (config/abacus_hod.yaml)
    2. Initialise AbacusHOD with simulation and HOD parameters
    3. Run the HOD -> mock_dict
    4. Save the catalog to a FITS file
    5. Add a unit-weight column to the catalog

Usage:
    From the project root:
        python scripts/generate_galaxy_mock.py
"""

import os

import yaml
import numpy as np
from astropy.table import Table
from abacusnbody.hod.abacus_hod import AbacusHOD

from functions import save_mock_dict_to_fits


# ---------------------------------------------------------------------------
# 1. Load configuration
# ---------------------------------------------------------------------------

# Path to the config file, relative to this script's location
path2config = os.path.join(os.path.dirname(__file__), "..", "config", "abacus_hod.yaml")

config = yaml.safe_load(open(path2config))
sim_params        = config['sim_params']
HOD_params        = config['HOD_params']
clustering_params = config['clustering_params']

# HOD run options
want_rsd      = HOD_params['want_rsd']       # Apply Redshift Space Distortions
write_to_disk = HOD_params['write_to_disk']  # Internal AbacusHOD disk write (disabled here)

# Projected separation bins for the correlation function
bin_params = clustering_params['bin_params']
rpbins     = np.logspace(bin_params['logmin'], bin_params['logmax'], bin_params['nbins'] + 1)
rp_centers = 0.5 * (rpbins[1:] + rpbins[:-1])  # Bin centres

pimax       = clustering_params['pimax']       # Max line-of-sight separation (Mpc/h)
pi_bin_size = clustering_params['pi_bin_size'] # Line-of-sight bin size (Mpc/h)


# ---------------------------------------------------------------------------
# 2. Initialise AbacusHOD
# ---------------------------------------------------------------------------

newBall = AbacusHOD(sim_params, HOD_params, clustering_params)


# ---------------------------------------------------------------------------
# 3. Run the HOD
# ---------------------------------------------------------------------------

# Returns a dict {tracer_name: {field: array, ...}, ...}
# e.g. {'LRG': {'x': array, 'y': array, 'z': array, 'vz': array, ...}}
mock_dict = newBall.run_hod(
    newBall.tracers,
    want_rsd,
    write_to_disk=write_to_disk,
    Nthread=16,
)


# ---------------------------------------------------------------------------
# 4. Save the FITS catalog
# ---------------------------------------------------------------------------

output_path = os.path.join(
    sim_params['output_dir'],
    f"mock_galaxy_catalogue_LRG_z{sim_params['z_mock']:.3f}.fits",
)

save_mock_dict_to_fits(mock_dict, output_path)
print(f"Catalog saved: {output_path}")


# ---------------------------------------------------------------------------
# 5. Add unit weights column
# ---------------------------------------------------------------------------

catalog = Table.read(output_path)
catalog['weights'] = np.ones(len(catalog))
catalog.write(output_path, overwrite=True)
print(f"'weights' column added ({len(catalog)} galaxies).")
