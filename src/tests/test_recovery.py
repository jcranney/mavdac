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


def dist_true(pos: np.ndarray) -> np.ndarray:
    """sample the true distortion function deterministically at points in
    the field"""
    dists = np.zeros_like(pos)
    dists[:, 0] = 1*((pos[:, 0] - 2000.0)/4000)**2
    return dists


def generate_image(shift_x: float, shift_y: float) -> fits.hdu.PrimaryHDU:
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
    pitch = 1.0
    width = 30.0
    overfill = 2.0
    xx, yy = np.meshgrid(
        np.arange(overfill*width/pitch),
        np.arange(overfill*width/pitch),
        indexing="xy",
    )
    xx = xx.flatten()
    yy = yy.flatten()
    xx -= overfill*width/pitch/2
    yy -= overfill*width/pitch/2
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
    p_pixels = p_pixels + np.random.randn(*p_pixels.shape)
    xx_pixels, yy_pixels = p_pixels
    xx_pixels += shift_x
    yy_pixels += shift_y
    # sample the true distortion at the actual positions of the pinholes after
    # shifting:
    dist_pixels = dist_true(np.array([xx_pixels, yy_pixels]).T)

    xx_pixels += dist_pixels[:, 0]
    yy_pixels += dist_pixels[:, 1]
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
        10_000, source, psfs_file="./test_psf_gauss.fits", which_psf=0,
    )
    imgen.main()
    im = imgen.get_rebinned_cropped(2, 30)
    im -= im.min()
    im /= im.max()
    im *= 30_000  # scale the image so that it's peak value is 30_000

    header = fits.Header()
    header["XSHIFT"] = shift_x
    header["YSHIFT"] = shift_y
    return fits.hdu.PrimaryHDU(data=im, header=header)


def run_mavdac(pattern: str, coordinate_file: str, distortion_file: str):
    """run the mavdac pipeline for a list of coordinates, and return the list
    of distortions evaluated at those coordinates"""
    task = subprocess.run([
        "mavdac",
        pattern,
        coordinate_file,
    ], stdout=subprocess.PIPE)
    with open(distortion_file, "w") as f:
        f.write(task.stdout.decode())


def test_mavdac():
    """run the pipeline described in the docstring of this file"""
    N_SAMPLES = 50  # number of unique points to evaluate the distortions
    PIXELS: int = 4000  # number of pixels across detector
    # define evaluation points
    p_eval = np.random.rand(N_SAMPLES, 2)*PIXELS
    # sample true distortion function at evaluation points
    d_true = dist_true(p_eval)

    SHIFT_RAD: float = 100.0  # pixels
    NIMAGES: int = 4
    shifts = []
    for theta in np.linspace(0, 2*np.pi, NIMAGES+1)[:-1]:
        shifts.append([SHIFT_RAD*np.cos(theta), SHIFT_RAD*np.sin(theta)])

    with tempfile.TemporaryDirectory() as d:
        # create images
        i = 0
        for shift_x, shift_y in shifts:
            hdu = generate_image(shift_x=shift_x, shift_y=shift_y)
            hdu.writeto(os.path.join(d, f"img_{i:03d}.fits"))
            i += 1

        # create coords
        coord_path = os.path.join(d, "coords.txt")
        with open(coord_path, "w") as f:
            for x, y in p_eval:
                f.write(f"{x},{y},\n")

        # run mavdac
        recov_path = os.path.join(d, "recov.txt")
        run_mavdac(os.path.join(d, "img_*.fits"), coord_path, recov_path)

        # read output
        with open(recov_path) as f:
            lines = f.readlines()
        lines = [
            [float(x) for x in line.split(",")[2:4]]
            for line in lines if len(line) > 0
        ]
        d_reco = np.array(lines)

    # compare d_true with d_reco
    print(d_true)
    print(d_reco)
    print(d_true.std(axis=0))
    print(d_reco.std(axis=0))
    print((d_true - d_reco).std(axis=0))
    return d_true, d_reco


d_true, d_reco = test_mavdac()
d_true = d_true - d_true.mean(axis=0)
d_reco = d_reco - d_reco.mean(axis=0)
print(f"quality: "
      f"{(((d_reco[:, 0] @ d_true[:, 0]) * (d_true.T @ d_true)[0, 0]**-1))}")
