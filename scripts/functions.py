"""
functions.py
------------
Utility functions for generating and saving galaxy mock catalogs
produced by AbacusSummit/AbacusHOD simulations.
"""

import time

import numpy as np
from astropy.io import fits


def save_mock_dict_to_fits(mock_dict, filename, overwrite=True):
    """
    Save a galaxy mock dictionary to a multi-extension FITS file.

    Each tracer (e.g. 'LRG', 'ELG') is stored in a separate BinTable
    extension. Scalar values are written to the extension header; 1D arrays
    become table columns. Multi-dimensional arrays are flattened to 1D.

    Parameters
    ----------
    mock_dict : dict
        Dictionary of the form {'LRG': {'x': array, ...}, 'ELG': {...}, ...}
        as returned by AbacusHOD.run_hod().
    filename : str
        Output path for the .fits file.
    overwrite : bool
        If True, overwrite an existing file (default: True).

    Returns
    -------
    str
        Path to the written file.
    """
    # Primary HDU with creation metadata
    prihdr = fits.Header()
    prihdr['CREATOR'] = 'save_mock_dict_to_fits'
    prihdr['DATE']    = time.strftime('%Y-%m-%d')
    prihdr['HISTORY'] = 'Mock catalog saved from mock_dict'
    hdul = [fits.PrimaryHDU(header=prihdr)]

    for tracer_name, tracer in mock_dict.items():
        cols = []
        header_meta = {}

        # Determine length N from the first 1D array found
        N = None
        for val in tracer.values():
            arr = np.asarray(val)
            if arr.ndim == 1:
                N = arr.shape[0]
                break
        if N is None:
            continue  # No 1D array found: skip this tracer

        for key, val in tracer.items():
            arr = np.asarray(val)

            # Scalars go into the FITS header
            if arr.shape == ():
                header_meta[key.upper()] = str(val)
                continue

            # Flatten multi-dimensional arrays
            if arr.ndim > 1:
                arr = arr.ravel()

            # Choose the FITS column format based on data type
            if arr.dtype.kind in ('U', 'S', 'O'):
                str_arr = arr.astype(str)
                maxlen = max(len(s) for s in str_arr) if str_arr.size > 0 else 1
                cols.append(fits.Column(name=key, format=f'{maxlen}A', array=str_arr))
            elif np.issubdtype(arr.dtype, np.integer):
                cols.append(fits.Column(name=key, format='K', array=arr.astype(np.int64)))
            else:
                cols.append(fits.Column(name=key, format='D', array=arr.astype(np.float64)))

        if not cols:
            continue

        # Create the BinTable extension for this tracer and attach metadata
        bintbl = fits.BinTableHDU.from_columns(cols, name=tracer_name)
        for meta_k, meta_v in header_meta.items():
            try:
                bintbl.header[meta_k] = meta_v
            except Exception:
                bintbl.header[meta_k[:8]] = meta_v  # Truncate key if too long
        hdul.append(bintbl)

    fits.HDUList(hdul).writeto(filename, overwrite=overwrite)
    return filename
