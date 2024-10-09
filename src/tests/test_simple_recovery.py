#!/usr/bin/env python
"""
The purpose of this test is to validate the distortion recovery accuracy of the
mavdac pipeline. The summary of this validation is as follows:
 1) generate a known distortion field, call it `dist_true`, and sample it at
    known positions across the field -> dist_true(p_eval)
 2) generate 3 noiseless images of a realistic calibration source with the true
    distortion present, each shifted by 10 pixels @ [0 deg, 120 deg, 240 deg]
    from the nominal pinhole position (note, no image at (0,0)-pixel shift),
 3) run these images through mavdac pipeline, and evaluate the recovered
    distortion function (dist_reco) at the same positions as in step 1,
    -> dist_reco(p_eval).
 4) compare the recovered and true distortion values, and confirm that they are
    within an acceptable range when ignoring global tip/tilt distortions.

Then, an additional test can be performed with an error in the shifting
quantity, and the comparison between truth and recovered distortions can be
made, showing an error in the plate scale modes but otherwise sufficient
accuracy.
"""
from astropy.io import fits
import mavisim
import numpy as np
import tempfile
import os
import subprocess
from mavdac import run_mavdac

psfs_file = "./src/tests/test_psf_gauss.fits"
gridfile = "./src/tests/grid.yaml"


def generate_image(
    shift_x: float, shift_y: float, dist_pixels: callable
) -> fits.hdu.PrimaryHDU:
    """generate an image of the calibration source with some shift in x and y,
    returning a formatted fits.FitsHDU"""

    from pydantic import BaseModel, ConfigDict

    class Source(BaseModel):
        model_config = ConfigDict(arbitrary_types_allowed=True)
        exp_time: float = None  # exposure time in seconds to simulate.
        star: np.ndarray = None  # unique ID of each star.
        flux: np.ndarray = None  # flux of each star.
        gauss_pos: np.ndarray = None  # X/Y-position of each star (in arcsec).
        gauss_cov: np.ndarray = None  # covariance of Gaussian kernel
        static_dist: np.ndarray = None  # static distortion to apply
        cov_mat: np.ndarray = None

    # create hex grid:
    pitch = 6.0
    width = 15.0
    overfill = 4.0
    xx, yy = np.meshgrid(
        pitch*np.arange(overfill*width/pitch),
        pitch*np.arange(overfill*width/pitch),
        indexing="xy",
    )
    xx = xx.flatten()
    yy = yy.flatten()
    xx -= overfill*width/2
    yy -= overfill*width/2
    p = np.concatenate([xx[None, :], yy[None, :]], axis=0)
    transform = np.array([
        [1.0, 0.5],
        [0.0, (3/4)**0.5]
    ])
    p = transform @ p
    p = p[:, ((p[0, :] > -width/2) & (p[0, :] < width/2))]
    p = p[:, ((p[1, :] > -width/2) & (p[1, :] < width/2))]
    p_pixels = (p + 30/2)*4000/30
    # add pinhole position error with 1 pixel of standard deviation:
    np.random.seed(1234)
    p_pixels = p_pixels + 1.0*np.random.randn(*p_pixels.shape)
    xx_pixels, yy_pixels = p_pixels
    xx_pixels += shift_x
    yy_pixels += shift_y
    # sample the true distortion at the actual positions of the pinholes after
    # shifting:
    xx_pixels -= 0.5
    yy_pixels -= 0.5
    d_pixels = dist_pixels(np.array([xx_pixels, yy_pixels]).T)
    xx_pixels += d_pixels[:, 0]
    yy_pixels += d_pixels[:, 1]
    xx_pixels += 0.5
    yy_pixels += 0.5
    xx = (xx_pixels * 30/4000)-30/2
    yy = (yy_pixels * 30/4000)-30/2

    nstars = xx.shape[0]
    source = Source(
        exp_time=1.0,
        star=np.arange(xx.shape[0]),
        flux=np.ones(nstars),
        gauss_pos=np.concatenate([xx[None, :], yy[None, :]], axis=0).T,
        gauss_cov=np.tile((np.eye(2)*1e-9)[None, :, :], reps=[nstars, 1, 1])
    )

    imgen = mavisim.ImageGenerator(
        12_000, source, psfs_file=psfs_file, which_psf=0,
    )
    imgen.main()
    im = imgen.get_rebinned_cropped(2, 30)
    im -= im.min()
    im /= im.max()
    im *= 30_000  # scale the image so that it's peak value is 30_000

    header = fits.Header()
    header["XSHIFT"] = shift_x
    header["YSHIFT"] = shift_y
    return (
        fits.hdu.PrimaryHDU(data=im, header=header),
        np.array([xx_pixels, yy_pixels]).T
    )


def test_mavdac():
    """run the pipeline described in the docstring of this file"""
    N_SAMPLES = 100  # number of unique points to evaluate the distortions
    RADIUS: int = 1000
    # define evaluation points
    p_eval = (np.random.rand(N_SAMPLES, 2)*2-1)*RADIUS+2000
    coeffs_true = np.array([
        [1.0, 0.0],
        [2.0, 0.0],
        [3.0, 4.0],
    ])

    SHIFT_RAD: float = 100.0  # pixels
    NIMAGES: int = 3
    shifts = []
    for theta in np.linspace(0, 2*np.pi, NIMAGES+1)[:-1]:
        shifts.append([SHIFT_RAD*np.cos(theta), SHIFT_RAD*np.sin(theta)])

    with tempfile.TemporaryDirectory() as d:
        # create images
        i = 0
        expected_cogs = []
        for shift_x, shift_y in shifts:
            hdu, cogs = generate_image(
                shift_x=shift_x, shift_y=shift_y,
                dist_pixels=lambda p: dist_eval(p, coeffs_true)
            )
            expected_cogs += list(cogs)
            hdu.writeto(os.path.join(d, f"img_{i:03d}.fits"))
            i += 1
        basis = run_mavdac(
            os.path.join(d, "img_*.fits"), rad=30, flux_thresh=10_000.0,
            gridfile=gridfile, poly_degree=3,
        )

    # sample distortion function at evaluation points
    d_true = dist_eval(p_eval, coeffs_true)
    d_est = np.array([
        basis.eval_xy(x, y)
        for x, y in p_eval
    ])
    err = ((d_true - d_est).std(axis=0)**2).mean()**0.5
    print(f"residual error: {err} pixels rms")
    assert err < 1e-4


def test_mavdac_cli():
    """run the pipeline described in the docstring of this file"""
    N_SAMPLES = 100  # number of unique points to evaluate the distortions
    RADIUS: int = 1000
    # define evaluation points
    p_eval = (np.random.rand(N_SAMPLES, 2)*2-1)*RADIUS+2000
    coeffs_true = np.array([
        [1.0, 0.0],
        [2.0, 0.0],
        [3.0, 4.0],
    ])

    SHIFT_RAD: float = 100.0  # pixels
    NIMAGES: int = 3
    shifts = []
    for theta in np.linspace(0, 2*np.pi, NIMAGES+1)[:-1]:
        shifts.append([SHIFT_RAD*np.cos(theta), SHIFT_RAD*np.sin(theta)])

    with tempfile.TemporaryDirectory() as d:
        # create images
        i = 0
        expected_cogs = []
        for shift_x, shift_y in shifts:
            hdu, cogs = generate_image(
                shift_x=shift_x, shift_y=shift_y,
                dist_pixels=lambda p: dist_eval(p, coeffs_true)
            )
            expected_cogs += list(cogs)
            hdu.writeto(os.path.join(d, f"img_{i:03d}.fits"))
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
                "--degree=3",
            ], stdout=f)
        with open(recov_path, "r") as f:
            recov_lines = f.readlines()
        d_est = []
        for line in recov_lines:
            if len(line) > 0:
                d_est.append(line.split(",")[2:4])
        d_est = np.array(d_est, dtype=np.float64)

    # sample distortion function at evaluation points
    d_true = dist_eval(p_eval, coeffs_true)
    d_true -= d_true.mean(axis=0)[None, :]
    d_est -= d_est.mean(axis=0)[None, :]
    err = ((d_true - d_est).std(axis=0)**2).mean()**0.5
    print(recov_lines)
    print(d_true)
    print(d_est)
    print(f"residual error: {err} pixels rms")
    assert err < 1e-4


def fun(x, y, ell):
    return ((x - 2000.0)/4000)**(ell+1)


def dist_eval(pos, coeffs):
    dist = np.array([
        np.array([[
            fun(x, y, ell) * coeffs[ell, 0], fun(x, y, ell) * coeffs[ell, 1]
        ] for ell in range(coeffs.shape[0])]).sum(axis=0)
        for x, y in pos
    ])
    return dist
