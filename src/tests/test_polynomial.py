#!/usr/bin/env python
import numpy as np
import pytest
from typing import Callable, List, Tuple
import mavdac
from scipy.interpolate import LinearNDInterpolator  # type: ignore
psfs_file = "./src/tests/test_psf_gauss.fits"
gridfile = "./src/tests/grid.yaml"
mavis_dists = "./src/tests/mavis_dists.txt"
imshape = (4000, 4000)


def make_eval_points(n: int = 40, r: float = 15) -> List[Tuple[float, float]]:
    yy, xx = np.meshgrid(
        *[np.linspace(-r, r, n)]*2,
        indexing="ij",
    )
    xx = xx.flatten()
    yy = yy.flatten()
    rr = (xx**2 + yy**2)**0.5
    xx = xx[rr < r]
    yy = yy[rr < r]
    xx = xx*4000/30+2000
    yy = yy*4000/30+2000
    return list(zip(xx, yy))


def dists_none(p: np.ndarray) -> np.ndarray:
    return p*0


def dists_poly(p: np.ndarray) -> np.ndarray:
    np.random.seed(1234)
    poly = mavdac.mavdac.BiVarPolyDistortions(3, imshape)
    coeffs = np.array(poly.coeffs)
    poly.load_coeffs(1.0*np.random.randn(*coeffs.shape))
    return np.array([
        poly.eval_xy(x, y)
        for x, y in p
    ])


def dists_mavis(p: np.ndarray) -> np.ndarray:
    dist_func = load_dists_from_file(mavis_dists)
    return dist_func(p)


def load_dists_from_file(file: str) -> Callable:
    with open(file) as f:
        lines = f.readlines()
    platescale = 0.736  # "/mm
    pixelsize = 7.5e-3  # "/pixel
    scale = platescale/pixelsize
    pos = []
    dist = []
    for line in lines:
        s = line.split()
        try:
            float(s[0])
        except ValueError:
            continue
        pos.append((
            float(s[4]), float(s[5]),
        ))
        dist.append((
            float(s[6])-float(s[4]), float(s[7])-float(s[5]),
        ))
    dist_func = LinearNDInterpolator(
        np.array(pos)*scale+imshape[0]/2,
        np.array(dist)*scale,
        fill_value=0.0,
    )
    return dist_func


distortions = [
    # dists_none,
    # dists_poly,
    dists_mavis,
]

degrees = [5, 10, 15, 20, 25, 30]


@pytest.mark.parametrize("distortions", distortions)
@pytest.mark.parametrize("degree", degrees)
def test_basis_functions(
    distortions: Callable,
    degree: int,
    plot: bool = False,
) -> None:
    p_eval = np.array(make_eval_points(40, 15.0))
    d_true = distortions(p_eval)

    basis = mavdac.mavdac.BiVarPolyDistortions(degree, imshape)
    # we have a sampleable set of basis functions, now sample them at the eval
    # points:
    samples = np.array([
        [basis.sample_xy(x, y, i) for i in range(len(basis.coeffs))]
        for x, y in p_eval
    ])
    coeffs = np.linalg.solve(samples.T @ samples, samples.T @ d_true)
    basis.load_coeffs(coeffs=coeffs)
    d_est = np.array([basis.eval_xy(x, y) for x, y in p_eval])

    def rms(dists):
        return (dists.std(axis=0)**2).mean()**0.5

    # sample distortion function at evaluation points
    d_true -= d_true.mean(axis=0)[None, :]
    d_est -= d_est.mean(axis=0)[None, :]
    coord_mat = np.array([
        p_eval[:, 0],
        p_eval[:, 1],
        p_eval[:, 0]*0+1
    ]).T
    ps_filt = np.eye(coord_mat.shape[0]) - \
        coord_mat @ np.linalg.solve(coord_mat.T @ coord_mat, coord_mat.T)
    print(f"degree: {degree}")
    print("tt-ps-removed rms")
    print(f"input dist: {rms(ps_filt @ d_true)} pixels rms")
    print(f"  est dist: {rms(ps_filt @ d_est)} pixels rms")
    d_err = ps_filt @ (d_true - d_est)
    err = rms(d_err)
    print(f"     error: {err} pixels rms")
    rms(d_true - d_est)
    if plot:
        sf = 100.0
        for i in range(len(p_eval)):
            plt.subplot(1, 2, 1)
            plt.arrow(
                p_eval[i, 0], p_eval[i, 1],
                sf*d_true[i, 0], sf*d_true[i, 1],
                color="r", alpha=0.5, width=10,
            )
            plt.arrow(
                p_eval[i, 0], p_eval[i, 1],
                sf*d_est[i, 0], sf*d_est[i, 1],
                color="b", alpha=0.5, width=10,
            )
            plt.subplot(1, 2, 2)
            plt.arrow(
                p_eval[i, 0], p_eval[i, 1],
                sf*d_err[i, 0], sf*d_err[i, 1],
                color="g", alpha=0.5, width=10,
            )
    else:
        # assert err < 6.66e-3  # 50 uas rms
        pass


if __name__ == "__main__":
    plot = False
    # called as a script, not from pytest
    if plot:
        import matplotlib.pyplot as plt
        plt.ion()
    for distortion in distortions:
        for degree in degrees:
            if plot:
                plt.figure(figsize=(12, 8))
            test_basis_functions(
                distortions=distortion, plot=plot, degree=degree
            )
            if plot:
                plt.subplot(1, 2, 1)
                plt.axis("square")
                plt.legend(["d_true", "d_est"])
                plt.subplot(1, 2, 2)
                plt.axis("square")
                plt.legend(["d_err"])
