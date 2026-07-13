# AbacusSummit simulations

This repository contains scripts to work with AbacusSummit N-body simulations.
The idea is to generate mock galaxy catalogues using the Halo Occupation
Distribution (HOD) model and to build matched halo catalogues with shape
information for weak lensing and intrinsic alignments analyses.


## The HOD model

The HOD (Halo Occupation Distribution) statistically describes how many galaxies
occupy a dark matter halo of a given mass. The model is based on Zheng+2007
with assembly bias extensions:

- **Central galaxies**: the probability that a halo of mass M hosts a central
  follows an error function in `log(M)`, parametrised by `logM_cut` and `sigma`.
- **Satellite galaxies**: their mean count follows a power law of slope `alpha`,
  activated above `kappa * M_cut`.
- **Assembly bias**: parameters `Acent`, `Asat`, `Bcent`, `Bsat` modulate the
  occupation as a function of halo concentration and local environment.

Three tracer types are supported: **LRG**, **ELG**, **QSO**.


## Configuration (`config/abacus_hod.yaml`)

The key parameters to adapt are in the `sim_params` section:

| Parameter       | Description                                         |
|-----------------|-----------------------------------------------------|
| `sim_name`      | AbacusSummit box name                               |
| `sim_dir`       | Directory containing the halo catalogs (`.asdf`)    |
| `output_dir`    | Output directory for the FITS mock catalog          |
| `subsample_dir` | Output directory for particle subsamples            |
| `z_mock`        | Redshift of the simulated snapshot                  |

LRG HOD parameters are found under `HOD_params.LRG_params`.
To enable ELG or QSO tracers, set `tracer_flags.ELG: True` (etc.).

See Yuan et al. 2023 for best-fit HOD parameters for LRGs.


## Dependencies

```
abacusnbody    # AbacusSummit package (HOD + halo catalog reader)
astropy        # FITS I/O and table handling
numpy          # Numerical computing
pyyaml         # YAML configuration file parsing
treecorr       # Two-point correlation functions (used in functions.py)
scipy          # Statistical distributions (used in functions.py)
pandas         # Data table utilities (used in functions.py)
pycorr         # Two-point correlation function estimator (used in scripts/clustering)
tqdm           # Progress bars
torch          # Used by some HOD/clustering scripts
```

Installation follows the AbacusHOD instructions:
https://github.com/abacusorg/abacusutils

`pycorr` (two-point correlation function estimator used in
`scripts/clustering`) is available at:
https://github.com/cosmodesi/pycorr


## Scripts

| Script | Description |
|--------|-------------|
| `scripts/prepare_sim.py` | **Step 1 (required, once per snapshot)** — reads raw halo catalogs and writes particle subsample files to `subsample_dir` |
| `scripts/generate_galaxy_mock.py` | **Step 2** — populates halos with the HOD model and writes a multi-extension FITS catalog (one extension per tracer: LRG, ELG) |
| `scripts/generate_halo_catalogue.py` | **Step 3** — matches halos to mock galaxies by ID, computes ellipticity components (e1, e2), and writes the enriched halo catalog to a FITS file |
| `functions/functions.py` | Shared utility functions used by Step 3: 3D axis-ratio sampling (`population_3D`), ellipsoid projection (`projection`), and the top-level wrapper (`simulator`). `simulator` takes `theta = [mu_tau_B, mu_tau_C, sigma_tau_B, sigma_tau_C]`: galaxy B/C semi-axes are drawn from an independent bivariate normal and rescale the halo's B/C axes, while the major axis (A) is kept equal to the halo's major axis |
| `scripts/clustering/run_wgg_box.py` | **Step 4 (mock)** — computes the projected correlation function wp(rp) for the mock galaxy catalogue in the periodic simulation box |
| `scripts/clustering/run_wgg_data.py` | **Step 4 (data)** — computes wp(rp) with jackknife error bars for a real DESI galaxy sample (data + randoms) |
| `notebooks/clustering.ipynb` | **Step 5** — plots wp(rp) for the LRG and ELG tracers, comparing DESI data against the AbacusHOD mock, using the `.npz` outputs of Step 4 |


## Usage

**Step 1 — Prepare the simulation** (run once per snapshot, before any mock generation):

```bash
python scripts/prepare_sim.py
```

Reads halo catalogs from `sim_dir` and writes particle subsample files to
`subsample_dir`. Required before running the HOD.

**Step 2 — Generate the galaxy mock catalog**:

```bash
python scripts/generate_galaxy_mock.py
```

Populates halos with galaxies using the HOD parameters in `abacus_hod.yaml`.
The output FITS file is written to `output_dir` with the name:
`mock_galaxy_catalogue_LRG+ELG_z<redshift>.fits`

Each tracer (LRG, ELG) is stored in a separate FITS extension with the
number of central galaxies recorded in the `NCENT` header keyword.

**Step 3 — Build the halo catalogue**:

```bash
python scripts/generate_halo_catalogue.py
```

Reads the galaxy mock produced in Step 2, matches the corresponding halos
from the AbacusSummit catalog, computes halo ellipticity components (e1, e2)
from the shape-tensor eigenvectors, and writes the result to:
`halo_catalogue_<tracer>_z<redshift>.fits`

By default `simulator` is called with `theta = [1, 1, 0, 0]` (no scatter),
so galaxy shapes are set equal to their host halo's shape. If the number of
matched halos exceeds `nb_halos_max` (2,000,000), the script randomly
subsamples down to that limit before computing ellipticities.

As with the clustering scripts, the input/output file paths and the path to
`functions/functions.py` are hardcoded and should be adapted to your
environment.

**Step 4 — Galaxy clustering (wp(rp))**:

```bash
python scripts/clustering/run_wgg_box.py    # mock catalogue, periodic box
python scripts/clustering/run_wgg_data.py   # real DESI data + randoms, with jackknife errors
```

Both scripts use `pycorr` to measure the projected two-point correlation
function wp(rp) in `rp`/`pi` bins:

- `run_wgg_box.py` measures wp(rp) directly on the mock galaxy catalogue
  (Step 2 output) inside the periodic simulation box, so no random
  catalogue is required.
- `run_wgg_data.py` measures wp(rp) on a real DESI data + randoms sample,
  converting (RA, DEC, Z) to comoving Cartesian coordinates, and estimates
  error bars via angular jackknife regions built with K-Means.

Results are saved as `.npz` files containing `rp`, `wgg` (and `err_jk` /
`cov_jk` for the data case).

Paths to the input catalogues and output directories in these two scripts
are currently hardcoded (e.g. `/n17data/corinaldi/...`) — update them to
match your own environment before running.

**Step 5 — Visualize the clustering results**:

```bash
jupyter notebook notebooks/clustering.ipynb
```

Loads the `.npz` files produced in Step 4 and plots `rp * wp(rp)` for
data vs. mock:

- LRG tracer, comparing three DESI redshift bins (`0.4 < z < 0.6`,
  `0.6 < z < 0.8`, `0.8 < z < 1.1`) against the matching AbacusHOD mock
  snapshots.
- ELG tracer, comparing the `0.8 < z < 1.6` DESI bin against the mock at
  `z = 1.1`.

As with the clustering scripts, the input file paths are hardcoded and
should be adapted to your environment.
