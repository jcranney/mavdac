# mavdac
MAVIS Differential Astrometric Calibration. This image processing pipeline is designed to ingest images acquired during astrometric calibration of (e.g.) [MAVIS](https://mavis-ao.org), and output the estimated distortion field derived from these images. The bulk of the pipeline is written in Rust, but the most useful functions and the CLI tool are wrapped in Python.

# Install
```bash
pip install mavdac
```
To check that everything is set up correctly, test the CLI, e.g.,:
```bash
mavdac --help
```
You should see something like:
```
$ mavdac --help
usage: mavdac [-h] [--radius RADIUS] [--thresh THRESH] [--grid GRID] [--yes] [--degree DEGREE] pattern [coordinates]

mavis differential astrometric calibrator. For more info, see https://github.com/jcranney/mavdac

positional arguments:
  pattern          glob pattern for matching fits files to use for calibration
  coordinates      file containing coordinates to map through distortions

options:
  -h, --help       show this help message and exit
  --radius RADIUS  radius around nominal pinhole position to perform centroiding
  --thresh THRESH  threshold of summed flux under which a centroid is ignored
  --grid GRID      yaml file containing grid geometry definition, creates default if not present
  --yes, -y        answer "Yes" by default (e.g., when using non-interactive shells)
  --degree DEGREE  maximum order of bivariate polynomial used in distortion fitting
```


# Quickstart
The simplest usage of this software is by calling
```bash
mavdac "./some/pattern*.fits" ./coordinates.csv
```
which tries to parse all images on disk matching the glob pattern argument, and the coordinates in the comma separated human-readable coordinates file. If the calibrations are succesful, then the recovered distortions will be sampled at the specified coordinates and printed to stdout. So, to capture the output, you can do (e.g.):
```bash
mavdac "./some/pattern*.fits" ./coordinates.csv > output.csv
```



# Usage

There are a few ways to use `mavdac`, e.g.,:
 - Through the Python CLI,
 - As a Python module (part of another python-based application),
 - As a Rust crate (part of another rust-based application),

I'll assume you are a standard user, hoping to perform astrometric calibration of an optical system using the MAVIS Differential Astrometric Calibration technique. In that case, the workflow is as follows:

 1) **Collect images** of calibration source into `./some/dir/img_0*.fits`,
 2) **Define coordinates** in pixel space to evaluate the measured distortions, e.g., `./some/dir/coords`:
 ```
 0,0,
 0,8,
 0,16,
 0,32,
 ...
 4096,4088,
 4096,4096,
 ```
 3) **Write a `grid.yaml` file**, specifying the geometry of the pinhole mask/calibration source.
 4) **Run `mavdac`** pipeline, specifying image filename pattern and coordinates filename path, piping `stdout` to some `./results` file.

Below we provide the requirements for each of the above tasks.
### 1) Collecting images
 - Each image must be saved in FITS format, as any of: `{u8, f32, f64, i32, u32}` type. Note that internally, all images are converted to `f64` upon loading.
 - Each file must contain only a single 2D array in the primary HDU, and must contain the FITS header keywords:
   - `XSHIFT`: shift in pixels along x-axis of the calibration mask for that exposure (from some nominal "home" position),
   - `YSHIFT`: shift in pixels along y-axis of the calibration mask for that exposure (from some nominal "home" position),
### 2) Defining coordinates
 - The internal mavdac pipeline will fit a set of distortion basis functions (e.g., 2D bivariate polynomials) to the measured distortions. If you run the CLI without a second argument, e.g.:
 ```bash
mavdac "./some/dir/imgs_*.fits"
 ```
 then the software will print some status updates to stderr, and eventually the coefficients of the fitted basis functions to stdout. This is mostly unusable, unless you write your own software to parse these coefficients. The intended usage is to provide a set of coordinates in pixel space, in which case you will get the distortions at those coordinates printed to stdout. E.g., with a coordinates file, `./some/coords`:
 ```
0,0,
0,2000,
2000,0,
2000,2000,
 ```
 you can run:
 ```bash
 mavdac "./some/dir/imgs_*.fits" ./some/coords
 ```
 which will output:
 ```
0,0,0.0966850540113886,0.3300609989212976,
0,2000,0.07252484678180265,0.007773570343281411,
2000,0,0.004342546567941594,0.26372826456406706,
2000,2000,-0.0018909241959087534,-0.0036487693998900166,
 ```
 to stdout. The expected usage is something like:
 ```bash
mavdac "./imgs_*.fits" ./coords > ./out
 ```
 which cats only the stdout to a file, `./out`. I've spent a bit of time thinking about the simplest user interface with this software, and that is the best I can come up with, but if you have a better idea (perhaps inspired by experience using similar pipelines), please contact me or raise an issue in this repo.

### 3) Writing a `grid.yaml` file
 The coordinates of the pinhole grid are specified by a `grid.yaml` file, which currently only supports the "Hex" geometry. Here is an example:
  ```yaml
 !Hex
 pitch: 133.33333333333334
 rotation: 0.0
 offset:
   x: 0.0
   y: 0.0
  ```
  All of these parameters are required, they are described below:
   - `pitch`: distance between adjacent pinholes (in pixel units),
   - `rotation`: orientation of hexagonal pinhole geometry (in radians),
   - `offset`: shift of overall grid in pixel units, in `x` and `y` dimensions. Note that for `x=0` and `y=0`, there is a point in the grid at the exact middle of the image.
  
  Note that if you specify a grid file that does not exist, you will be asked if you want the default one. If you're feeling lazy, just say yes and then modify the default one that gets created.

### 4) Running `mavdac`
Finally, once you have:
 - the images in `./some/dir`,
 - a `grid.yaml` file,
 - and optionally, a `./coords` text file,
 
 you can run mavdac from the command line by:
```bash
mavdac "./some/dir/*.fits" ./coords > ./out
```


# Validating mavdac
All versions of mavdac < 1.0.0 are subject to breaking changes between any two releases. The tests and validations are being developed side-by-side with the software, so will become populated as the pipeline matures. It is envisioned that there will be 3 levels of validation for the software that should be re-satisfied for any release:
 - unit tests,
 - system tests, and
 - hardware tests.

## Unit tests
Unit tests for the rust package will be found within each rust module, e.g., within `./rust/src/lib.rs`, under `mod tests {...}`.

## System tests
System tests for the rust + python package (e.g., the CLI), can be found in `./src/tests`. These will include tests that generate images with known distortions that are piped through mavdac, and the measured distortion field can be checked for quality. These are already quite mature, and can be tested in any version `>=0.1.1` by running
```bash
pytest
```
in the root directory of this repo.

## Hardware tests
The final validation of this software package is through it's use on data acquired in a real optical system with introduced distortions. To demonstrate the challenges associated with this kind of validation - try to construct an experiment which can verify that the software works as intended. The following are the most direct:
- **measurement of known distortion**: By injecting a known distortion source into the optical path of a simple imaging system, one can estimate the level of distortion introduced by this optic, and compare it to the truth. The challenge here is finding some optical component that introduces a known distortion with some high level of absolute precision.
- **measurement of known calibration source**: By using a highly calibrated pinhole mask at the entrance of the optical system, one can perform a classical distortion calibration, measuring the observed spots and taking their offset as the "true distortion", then compare this with the recovered distortion using the differential calibration method. The challenge here is finding a calibration source which is accurate and stable down to 10s of nanometres. It is worth noting that the calibration object does not need to be the full pinhole grid, but could be a smaller set of sources with precise relative positioning that are moved around the field through different exposures.
- The third and most flexible method we considered allows the use of an **unknown distortion** with an **unknown calibration source** - that is, a calibration source with some unknown pinhole position errors, then to characterise the distortions with this mask. This will allow us to estimate the true positions of the pinholes (assuming that the mavdac pipeline works well). Then, we can insert the pinhole mask in a different orientation (e.g., rotated by 180 degrees), and measure the observed pinhole positions. These new pinhole positions should be the original estimated pinhole positions with distortion removed, rotated by 180 degress, plus the recovered distortion evaluated at the new pinhole positions. This should first be done without a strong distortion in the system, and then repeated with stronger distortions to assess the practical limitations of this method.


# todo:
 - provide an example of using mavdac + mavisim, simulating realistic images and distortions for MAVIS,
