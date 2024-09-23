import argparse
import mavdac

parser = argparse.ArgumentParser(
    "process mavis differential astrometric calibrations"
)
parser.add_argument(
    "pattern", type=str,
    help="glob pattern for matching fits files to use for calibration"
)
parser.add_argument(
    "coordinates", type=str, nargs="?",
    help="file containing coordinates to map through distortions",
)

if __name__ == "__main__":
    args = parser.parse_args()
    mavdac.cli_diff(args.pattern, args.coordinates)
