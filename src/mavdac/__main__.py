import argparse
import mavdac
import json

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
    help="yaml file containing grid geometry definition",
)
parser.add_argument(
    "--degree", type=int, default=5,
    help="maximum order of bivariate polynomial used in distortion fitting"
)

args = parser.parse_args()

basis = mavdac.run_mavdac(
    args.pattern, rad=args.radius, flux_thresh=args.thresh,
    gridfile=args.grid, poly_degree=args.degree
)

if args.coordinates:
    coordinates = mavdac.get_coordinates(args.coordinates)
    for coordinate in coordinates:
        posx, posy = coordinate.pos
        distx, disty = basis.eval_xy(posx, posy)
        print(f"{posx},{posy},{distx},{disty}")
else:
    # no coordinates to eval, so lets just print the coefficients
    print(json.dumps(basis.coeffs, indent=2))
