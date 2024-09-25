# MAVDAC
MAVIS Differential Astrometric Calibration.

# Install
Eventually this will be pip installable, but for now, the easiest way to get going
is through `maturin` (which is `pip install`able). Clone this repo then build it:
```bash
git clone https://github.com/jcranney/mavdac
cd mavdac
maturin develop --release
```
If successful, you should now be able to run `mavdac` as per [Usage](#usage)

# Usage
There are a few ways to use `mavdac`, e.g.,:
 - Through the Python CLI,
 - As a Python module (part of another python-based application),
 - As a Rust crate (part of another rust-based application),

Here, I assume you are a standard user, hoping to perform astrometric calibration of an optical system using the MAVIS Differential Astrometric Calibration technique. In that case, the workflow is as follows:

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

### 4) Running `mavdac`
Finally, once you have:
 - the images in `./some/dir`,
 - a `grid.yaml` file,
 - and optionally, a `./coords` text file,
 
 you can run mavdac from the command line by:
```bash
mavdac "./some/dir/*.fits" ./coords > ./out
```
For a reminder, you can check `mavdac --help`:
mavis_acm/lab/inputs 
```
$ python -m mavdac --help
usage: process mavis differential astrometric calibrations [-h] pattern [coordinates]

positional arguments:
  pattern      glob pattern for matching fits files to use for calibration
  coordinates  file containing coordinates to map through distortions

options:
  -h, --help   show this help message and exit
```