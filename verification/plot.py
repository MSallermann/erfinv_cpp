import numpy as np
import matplotlib.pyplot as plt


def add_plot_digits(data_file_name):
    data = np.loadtxt(data_file_name)
    x = data[:, 0]
    n_digits = data[:, 3]
    plt.plot(x, n_digits, label=data_file_name)


def add_plot_erfinv(data_file_name):
    data = np.loadtxt(data_file_name)
    x = data[:, 0]
    erfinv_x = data[:, 1]
    plt.plot(x, erfinv_x, label=data_file_name)


data_files = [
    "series_5.txt",
    "series_10.txt",
    "newton.txt",
    "series_10_1halley.txt",
    "series_10_2halley.txt",
    # "halley.txt",
    # "series_15.txt",
    # "series_20.txt",
    # "winitzki.txt",
    # "full.txt",
    # "lut_100.txt",
]

for df in data_files:
    add_plot_digits(df)

plt.xlabel("x")
plt.ylabel("n_digits")
plt.legend()
plt.savefig("n_digits.png", dpi=300)
plt.show()

plt.close()
for df in data_files:
    add_plot_erfinv(df)

plt.xlabel("x")
plt.ylabel("erfinv")
plt.legend()
plt.show()
plt.savefig("erfinv.png", dpi=300)
