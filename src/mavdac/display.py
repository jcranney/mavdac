import matplotlib.pyplot as plt


def plot_input(filename):
    with open(filename, "r") as f:
        lines = f.readlines()
    lines = [ell for ell in lines if len(ell) > 0]
    fig, ax = plt.subplots(1, 1, figsize=[8, 8])
    for line in lines:
        try:
            row = [float(a) for a in line.split(" ")]
        except ValueError:
            continue
        x = row[0]*3600*4000/30+2000
        y = row[1]*3600*4000/30+2000
        dx = (row[6] - row[4])*0.736*4000/30
        dy = (row[7] - row[5])*0.736*4000/30
        plt.arrow(x, y, dx, dy)


def plot_output(filename):
    with open(filename, "r") as f:
        lines = f.readlines()
    lines = [ell for ell in lines if len(ell) > 0]
    fig, ax = plt.subplots(1, 1, figsize=[8, 8])
    for line in lines:
        x, y, dx, dy = [float(a) for a in line.split(",")[:4]]
        plt.arrow(x, y, dx, dy)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("distortions")
    args = parser.parse_args()
    suffix = args.distortions.split(".")[-1]
    if suffix == "txt":
        plot_input(args.distortions)
    else:
        plot_output(args.distortions)
    plt.show()
