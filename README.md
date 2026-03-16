# AbacusSummit simulations

This repository contains some scripts and notebooks to work with AbacusSummit N-body simulations. In particular, it is possible
to generate mock galaxy catalogues using the HOD model.



### The HOD model

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
```

Installation follows the AbacusHOD instructions:
https://github.com/abacusorg/abacusutils


## Scripts

| Script | Description |
|--------|-------------|
| `scripts/prepare_sim.py` | **Step 1 (required, once per snapshot)** — reads raw halo catalogs and writes particle subsample files to `subsample_dir` |
| `scripts/generate_galaxy_mock.py` | **Step 2** — populates halos with the HOD model and writes the galaxy catalog to a FITS file |
| `scripts/functions.py` | Shared utility functions (e.g. `save_mock_dict_to_fits`) |

## Usage

**1. Prepare the simulation** (run once per snapshot, before any mock generation):

```bash
python scripts/prepare_sim.py
```

This step reads halo catalogs from `sim_dir` and writes particle subsample
files to `subsample_dir`. It is required before running the HOD.

**2. Generate the galaxy catalog**:

```bash
python scripts/generate_galaxy_mock.py
```

The output catalog is written to `output_dir` with the name:
`mock_galaxy_catalogue_LRG_z<redshift>.fits`


