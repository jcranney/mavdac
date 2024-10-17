#!/usr/bin/env python
import numpy as np
import tempfile
import pytest
import os
import subprocess
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


def generate_image_mavisim_cli(
    shift_x: float, shift_y: float, *, dist_pixels: Callable,
    pinholes: List[Tuple[float, float]], filename: str
) -> None:
    """generate an image of the calibration source with some shift in x and y
    returning the name of the file"""
    xx_pixels, yy_pixels = np.array(pinholes).T
    xx_pixels += shift_x
    yy_pixels += shift_y
    d_pixels = dist_pixels(np.array([xx_pixels, yy_pixels]).T)
    xx_pixels += d_pixels[:, 0]
    yy_pixels += d_pixels[:, 1]
    xx_as = (xx_pixels * 30/imshape[0])-30/2
    yy_as = (yy_pixels * 30/imshape[0])-30/2
    header = f"--header={{\"xshift\":{shift_x},\"yshift\":{shift_y}}}"
    with tempfile.NamedTemporaryFile("w") as f:
        f.write("Star RA Dec X PM_X Y PM_Y Flux\n")
        for i, star in enumerate(zip(xx_as, yy_as)):
            f.write(f"{i} 0.0 0.0 {xx_as[i]} 0.0 {yy_as[i]} 0.0 100000.0\n")
        result = subprocess.run([
            "mavisim",
            f"-o={filename}",
            "--static=0",
            header,
            f.name,
            psfs_file,
        ])
        if result.returncode != 0:
            raise RuntimeError(
                f"mavisim failed, exitcode: {result.returncode}"
            )


def grid_perfect() -> List[Tuple[float, float]]:
    return grid_perturbed(0.0)


def grid_perturbed(std: float) -> List[Tuple[float, float]]:
    grid = mavdac.mavdac.Grid(gridfile)
    points = grid.all_points(*imshape)
    np.random.seed(1234)
    return [(
        point.x + std*np.random.randn(),
        point.y + std*np.random.randn(),
        ) for point in points]


def grid_scaled(std: float, scale: float) -> List[Tuple[float, float]]:
    return [(x*scale, y*scale) for (x, y) in grid_perturbed(std)]


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


@pytest.mark.parametrize("pinholes", [
    # grid_perfect(),
    # grid_perturbed(1.0),
    grid_scaled(1.0, 1.001)
])
@pytest.mark.parametrize("distortions", [
    # dists_none,
    dists_poly,
    dists_mavis
])
def test_mavdac_cli(
    pinholes: List[Tuple[float, float]], distortions: Callable
) -> None:
    """run the pipeline described in the docstring of this file"""
    p_eval = np.array(make_eval_points(40, 14.0))
    SHIFT_RAD: float = 50.0  # pixels
    NIMAGES: int = 5
    shifts = []
    for theta in np.linspace(0, 2*np.pi, NIMAGES+1)[:-1]:
        shifts.append([SHIFT_RAD*np.cos(theta), SHIFT_RAD*np.sin(theta)])

    with tempfile.TemporaryDirectory() as d:
        # create images
        i = 0
        for shift_x, shift_y in shifts:
            image_filename = os.path.join(d, f"img_{i:03d}.fits")
            generate_image_mavisim_cli(
                shift_x=shift_x, shift_y=shift_y,
                dist_pixels=lambda p: distortions(p)*0.1,
                filename=image_filename, pinholes=pinholes,
            )
            i += 1

        # create coords
        coord_path = os.path.join(d, "coords.txt")
        with open(coord_path, "w") as f:
            for x, y in p_eval:
                a = f"{x:},{y},\n"
                f.write(a)

        # run mavdac
        recov_path = os.path.join(d, "recov.txt")
        with open(recov_path, "w") as f:
            subprocess.run([
                "mavdac",
                os.path.join(d, "img_*.fits"),
                coord_path,
                f"--grid={gridfile}",
                "--radius=30",
                "--thresh=10000",
                "--degree=20",
            ], stdout=f)
        with open(recov_path, "r") as f:
            recov_lines = f.readlines()
        recov_lines = [line for line in recov_lines if len(line) > 0]
        d_est = np.zeros([len(recov_lines), 2])
        for ell, line in enumerate(recov_lines):
            d_est[ell] = [float(val) for val in line.split(",")[2:4]]
        d_est = np.array(d_est, dtype=np.float64)

    def rms(dists):
        return (dists.std(axis=0)**2).mean()**0.5

    # sample distortion function at evaluation points
    d_true = distortions(p_eval)
    d_true -= d_true.mean(axis=0)[None, :]
    d_est -= d_est.mean(axis=0)[None, :]
    # tip-tilt removed error:
    print("tt-removed rms")
    print(f"input dist: {rms(d_true)} pixels rms")
    print(f"  est dist: {rms(d_est)} pixels rms")
    err = rms(d_true - d_est)
    print(f"     error: {rms(d_true - d_est)} pixels rms")
    # assert err < 1e-3
    # plate scale removed error:
    coord_mat = np.array([
        p_eval[:, 0],
        p_eval[:, 1],
        p_eval[:, 0]*0+1
    ]).T
    ps_filt = np.eye(coord_mat.shape[0]) - \
        coord_mat @ np.linalg.solve(coord_mat.T @ coord_mat, coord_mat.T)
    print("tt-ps-removed rms")
    print(f"input dist: {rms(ps_filt @ d_true)} pixels rms")
    print(f"  est dist: {rms(ps_filt @ d_est)} pixels rms")
    d_err = ps_filt @ (d_true - d_est)
    err = rms(d_err)
    print(f"     error: {err} pixels rms")
    rms(d_true - d_est)
    assert err < 1e-5
