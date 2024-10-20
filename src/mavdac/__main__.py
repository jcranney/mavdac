import argparse
import mavdac
import json
import os
import click
import sys
import glob
import numpy

parser = argparse.ArgumentParser(
    "mavdac",
    # usage="\n  mavdac <pattern> [coordinates]",
    description="""
    mavis differential astrometric calibrator.
    For more info, see https://github.com/jcranney/mavdac
    """
)
parser.add_argument(
    "pattern", type=str,
    help="glob pattern for matching fits files to use for calibration"
)
parser.add_argument(
    "coordinates", type=str, nargs="?",
    help="file containing coordinates to map through distortions",
)
parser.add_argument(
    "--radius", type=int, default=30,
    help="radius around nominal pinhole position to perform centroiding",
)
parser.add_argument(
    "--thresh", type=float, default=10e3,
    help="threshold of summed flux under which a centroid is ignored"
)
parser.add_argument(
    "--grid", type=str, default="grid.yaml",
    help=(
        "yaml file containing grid geometry definition, "
        "creates default if not present"
    ),
)
parser.add_argument(
    "--yes", "-y", action="count", default=0,
    help="answer \"Yes\" by default (e.g., when using non-interactive shells)"
)
parser.add_argument(
    "--degree", type=int, default=5,
    help="maximum order of bivariate polynomial used in distortion fitting"
)

args = parser.parse_args()

if len(glob.glob(args.pattern)) == 0:
    print(
        f"pattern does not match any files: {args.pattern}\naborting.",
        file=sys.stderr,
    )
    exit(1)

if not os.path.isfile(args.grid):
    if args.yes == 0 and not click.confirm(
        f"{args.grid} does not exist.\n"
        "would you like to use the default one and proceed?", default=True,
        err=True,
    ):
        print("aborting", file=sys.stderr)
        exit(2)
    mavdac.mavdac.Grid.Hex(
        pitch=100.0, rotation=0.0,
        offset=mavdac.mavdac.Vec2D(0.0, 0.0)
    ).to_yaml(args.grid)

try:
    basis = mavdac.run_mavdac(
        args.pattern, rad=args.radius, flux_thresh=args.thresh,
        gridfile=args.grid, poly_degree=args.degree
    )
except numpy.linalg.LinAlgError:
    print(
        """error: singular matrix, not enough independent observations.
try some of the following:
    - using a more diverse set of XSHIFT/YSHIFT,
    - using more pinholes (reduce the flux threshold),
    - using a smaller degree polynomial basis,
    - raising an issue at https://github.com/jcranney/mavdac/issues/new""",
        file=sys.stderr,
    )
    exit(3)

if args.coordinates:
    coordinates = mavdac.mavdac.get_coordinates(args.coordinates)
    for coordinate in coordinates:
        posx, posy = coordinate.pos
        distx, disty = basis.eval_xy(posx, posy)
        print(f"{posx},{posy},{distx},{disty}")
else:
    # no coordinates to eval, so lets just print the coefficients
    print(json.dumps(basis.coeffs, indent=2))
