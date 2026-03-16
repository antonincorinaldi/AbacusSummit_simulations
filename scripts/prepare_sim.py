"""
prepare_sim.py
--------------
Preprocessing script for AbacusSummit simulations, to be run before any mock
generation.

This step must be executed once per snapshot before calling
generate_galaxy_mock.py. It reads the raw halo catalogs (.asdf files in
sim_dir) and writes the particle subsample files required by the AbacusHOD
framework to subsample_dir (both paths defined in config/abacus_hod.yaml).

Pipeline overview:
    [1] prepare_sim.py          <- THIS SCRIPT (run first)
    [2] generate_galaxy_mock.py <- HOD galaxy catalog generation

Usage:
    From the project root:
        python scripts/prepare_sim.py

    The path to the configuration file is currently hard-coded in main().
    Update `path2config` if the project is moved to a different machine.
"""

# --- AbacusSummit dependencies ---
from abacusnbody.hod import prepare_sim   # Particle subsample preparation module

# --- Unused imports (kept for potential future extensions) ---
from astropy.io import fits
from astropy.table import Table
from astropy import units as u
from abacusnbody.data.compaso_halo_catalog import CompaSOHaloCatalog
from astropy.cosmology import FlatLambdaCDM
from abacusnbody.hod.abacus_hod import AbacusHOD
import matplotlib.pyplot as plt
import os
import time
import yaml
import numpy as np


def main():
    """
    Main entry point of the script.

    Calls prepare_sim.main() from abacusnbody with the YAML configuration
    file. That function:
      - loads simulation parameters (sim_params) and preparation options
        (prepare_sim) from the YAML;
      - reads the halo catalogs (CompaSOHaloCatalog) from sim_dir;
      - writes the particle subsample files to subsample_dir.

    Notes
    -----
    The path to the configuration file is hard-coded. Update `path2config`
    if the project is moved to a different machine or directory layout.
    """

    # Absolute path to the AbacusHOD configuration file.
    # Update this path if the directory layout changes.
    path2config = "/n17data/corinaldi/halos/Abacus_HOD/abacusutils/abacus_hod.yaml"

    # Run the preparation: read halo catalogs and write particle subsamples
    prepare_sim.main(path2config)


if __name__ == "__main__":
    main()
