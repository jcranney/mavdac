#!/usr/bin/env python
"""
python utilities for MAVDAC
"""
import numpy as np
import mavdac


def coeffs_from_cogs(
    cogs, basis_function: mavdac.mavdac.BiVarPolyDistortions
) -> None:
    """Take a Vec<Vec<Centroid>> and calculate the coefficients of the
    basis functions.
    """
    fun = basis_function.sample_xy
    P = len(cogs)  # number of pinholes
    if P == 0:
        raise ValueError("no valid cogs, perhaps the threshold is too high?")
    N = len(cogs[0])  # number of shifts
    L = basis_function.ncoeffs  # number of polynomials
    R = 2  # number of dimensions
    d_ijkd = np.zeros([P, N, N, R])
    p_ijkd = np.zeros([P, N, N, R])
    f_ijkl = np.zeros([P, N, N, L])
    for i, cogs_pinhole in enumerate(cogs):
        for j, cogij in enumerate(cogs_pinhole):
            for k, cogik in enumerate(cogs_pinhole):
                d_ijkd[i, j, k, 0] = cogik.cogx - cogij.cogx
                d_ijkd[i, j, k, 1] = cogik.cogy - cogij.cogy
                p_ijkd[i, j, k, 0] = cogik.posx - cogij.posx
                p_ijkd[i, j, k, 1] = cogik.posy - cogij.posy
                for ell in range(L):
                    f_ijkl[i, j, k, ell] = (
                        fun(cogik.posx, cogik.posy, ell) -
                        fun(cogij.posx, cogij.posy, ell)
                    )
    b = d_ijkd - p_ijkd
    a = f_ijkl
    b = b.reshape([P*N*N, R]).T
    a = a.reshape([P*N*N, L]).T
    coeffs = np.linalg.solve(a @ a.T, a @ b.T)
    basis_function.load_coeffs(coeffs)


def run_mavdac(
    pattern: str, *, rad: int, flux_thresh: float, gridfile: str,
    poly_degree: int,
) -> mavdac.mavdac.BiVarPolyDistortions:
    """run the mavdac pipeline for a pattern of images, return sampleable
    distortion basis function"""
    imgs = mavdac.mavdac.load_images(pattern)
    im_shape = imgs[0].shape
    grid = mavdac.mavdac.Grid(gridfile)
    cogs = mavdac.mavdac.measure_cogs(imgs, grid, rad, flux_thresh)
    basis = mavdac.mavdac.BiVarPolyDistortions(poly_degree, im_shape)
    mavdac.coeffs_from_cogs(cogs, basis)
    imgs[0].draw_on_circles(grid, rad, 500)
    imgs[0].to_fits("/tmp/sample.fits")
    return basis
