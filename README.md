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


## Dependencies

```
abacusnbody    # AbacusSummit package (HOD + halo catalog reader)
astropy        # FITS I/O and table handling
numpy          # Numerical computing
pyyaml         # YAML configuration file parsing
treecorr       # Two-point correlation functions (used in functions.py)
IACorr         # Intrinsic-alignment correlation utilities (used in functions.py)
scipy          # Statistical distributions (used in functions.py)
pandas         # Data table utilities (used in functions.py)
```

Installation follows the AbacusHOD instructions:
https://github.com/abacusorg/abacusutils


## Scripts

| Script | Description |
|--------|-------------|
| `scripts/prepare_sim.py` | **Step 1 (required, once per snapshot)** — reads raw halo catalogs and writes particle subsample files to `subsample_dir` |
| `scripts/generate_galaxy_mock.py` | **Step 2** — populates halos with the HOD model and writes a multi-extension FITS catalog (one extension per tracer: LRG, ELG) |
| `scripts/generate_halo_catalogue.py` | **Step 3** — matches halos to mock galaxies by ID, computes ellipticity components (e1, e2), and writes the enriched halo catalog to a FITS file |
| `functions/functions.py` | Shared utility functions used by Step 3: 3D axis-ratio sampling (`population_3D`), ellipsoid projection (`projection`), and the top-level wrapper (`simulator`) |


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
