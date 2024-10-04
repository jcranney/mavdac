import matplotlib.pyplot as plt


def plot_output(filename):
    with open(filename, "r") as f:
        lines = f.readlines()
    lines = [ell for ell in lines if len(ell) > 0]
    fig, ax = plt.subplots(1, 1, figsize=[8, 8])
    for line in lines:
        x, y, dx, dy = [float(a) for a in line.split(",")[:4]]
        plt.arrow(x, y, dx, dy, width=20)


if __name__ == "__main__":
    plot_output("out")
    plt.show()
