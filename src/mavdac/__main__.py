import argparse
import mavdac

parser = argparse.ArgumentParser(
    "mavis differential astrometric calibrator",
    usage="\n  mavdac <pattern> [coordinates]",
    description="For more info, see https://github.com/jcranney/mavdac"
)
parser.add_argument(
    "pattern", type=str,
    help="glob pattern for matching fits files to use for calibration"
)
parser.add_argument(
    "coordinates", type=str, nargs="?",
    help="file containing coordinates to map through distortions",
)

args = parser.parse_args()
mavdac.cli_diff(args.pattern, args.coordinates)
