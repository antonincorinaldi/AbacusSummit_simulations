# Compute the projected galaxy clustering signal wp(rp), with jackknife
# error bars, for a real DESI spectroscopic sample (data + randoms).
# Unlike run_wgg_box.py, angular positions (RA/DEC/Z) are converted to
# comoving Cartesian coordinates and a data-random estimator is used since
# the survey footprint is not periodic.

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
from pycorr import TwoPointCorrelationFunction, KMeansSubsampler



# Flat LCDM cosmology used to convert redshifts to comoving distances.
# H0=100 means distances are expressed in Mpc/h.
cosmo = FlatLambdaCDM(H0=100, Om0=0.3158)
#cosmo = FlatLambdaCDM(H0=100, Om0=0.308)

tracer = 'BGS_BRIGHT'
zmin = 0.1; zmax = 0.3



# Data: load the DESI galaxy sample and restrict it to the redshift range
# of interest, then convert (RA, DEC, comoving distance) to Cartesian
# coordinates (x, y, z). FKP + systematic weights are combined into a
# single per-galaxy weight for the correlation function estimator.
galaxy_sample = Table.read(f'/n17data/corinaldi/DESI/clustering_data/{tracer}/{tracer}_data_NGC+SGC.fits')
galaxy_sample = galaxy_sample[(galaxy_sample['Z'] > zmin) & (galaxy_sample['Z'] < zmax)].copy()
chi = cosmo.comoving_distance(galaxy_sample['Z']).value
x_d = chi * np.cos(galaxy_sample['RA'].value * np.pi / 180) * np.cos(galaxy_sample['DEC'].value * np.pi / 180)
y_d = chi * np.sin(galaxy_sample['RA'].value * np.pi / 180) * np.cos(galaxy_sample['DEC'].value * np.pi / 180)
z_d = chi * np.sin(galaxy_sample['DEC'].value * np.pi / 180)
pos_data = np.array([x_d, y_d, z_d])
w_data = galaxy_sample['WEIGHT'] * galaxy_sample['WEIGHT_FKP']


# Randoms: same redshift cut and coordinate conversion as the data, used to
# estimate and correct for the survey selection function/geometry.
random_galaxy_sample = Table.read(f'/n17data/corinaldi/DESI/randoms/{tracer}/{tracer}_random_NGC+SGC.fits')
random_galaxy_sample = random_galaxy_sample[(random_galaxy_sample['Z'] > zmin) & (random_galaxy_sample['Z'] < zmax)].copy()
chi_r = cosmo.comoving_distance(random_galaxy_sample['Z']).value
x_r = chi_r * np.cos(random_galaxy_sample['RA'].value * np.pi / 180) * np.cos(random_galaxy_sample['DEC'].value * np.pi / 180)
y_r = chi_r * np.sin(random_galaxy_sample['RA'].value * np.pi / 180) * np.cos(random_galaxy_sample['DEC'].value * np.pi / 180)
z_r = chi_r * np.sin(random_galaxy_sample['DEC'].value * np.pi / 180)
pos_randoms = np.array([x_r, y_r, z_r])
w_randoms = random_galaxy_sample['WEIGHT'] * random_galaxy_sample['WEIGHT_FKP']



# rp/pi binning: nbins log-spaced radial bins (transverse separation) and
# npi linear bins on each side of the line of sight, integrated to obtain wp.
nbins  = 20
rmin   = 1.
rmax   = 200.
pi_max = 60.
npi    = 10

rp_edges = np.geomspace(rmin, rmax, nbins + 1)
pi_edges = np.linspace(-pi_max, pi_max, 2 * npi + 1)

# ---- Jackknife: angular regions defined via K-Means clustering ----
# Data and randoms are labelled with the same angular regions so that each
# jackknife realisation can drop one region consistently from both samples.
nsamples = 50

subsampler = KMeansSubsampler(
      'angular',
      positions=pos_data,
      nsamples=nsamples,
      position_type='xyz',
      random_state=42,
  )
samples_data    = subsampler.label(pos_data,     position_type='xyz')
samples_randoms = subsampler.label(pos_randoms,  position_type='xyz')

# ---- Correlation function computation (data-random estimator) ----
# los='midpoint' uses the angular midpoint between each pair as the
# line-of-sight direction, appropriate for a wide-angle survey geometry.
result = TwoPointCorrelationFunction(
      'rppi',
      edges=(rp_edges, pi_edges),
      data_positions1=pos_data,
      data_weights1=w_data,
      data_samples1=samples_data,
      randoms_positions1=pos_randoms,
      randoms_weights1=w_randoms,
      randoms_samples1=samples_randoms,
      los='midpoint',
      engine='corrfunc',
      nthreads=30
  )

# ---- wp + jackknife error bars ----
# Integrate xi(rp, pi) over pi to get wp(rp); wp_jk is the jackknife
# covariance matrix, from which the diagonal gives the 1-sigma errors.
rp, wgg, wp_jk = result.get_corr(return_sep=True, mode='wp')
err_jk  = np.sqrt(np.diag(wp_jk))

np.savez(
      f"/n17data/corinaldi/halos/Abacus_HOD/pycorr/results_data/wgg_{tracer}_z{zmin}_{zmax}_pycorr.npz",
      rp=rp,
      wgg=wgg,
      err_jk=err_jk,
      cov_jk=wp_jk
  )



