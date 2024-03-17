import numpy as np
import matplotlib.pyplot as plt

data_series_5 = np.loadtxt("series_5.txt")
data_series_10 = np.loadtxt("series_10.txt")
data_series_15 = np.loadtxt("series_15.txt")
data_series_20 = np.loadtxt("series_20.txt")
data_winitzki = np.loadtxt("winitzki.txt")
data_full = np.loadtxt("full.txt")
data_lut = np.loadtxt("lut_100.txt")

data_list = [
    data_series_5,
    data_series_10,
    data_series_15,
    data_series_20,
    data_winitzki,
    data_full,
    # data_lut,
]
label_list = [
    "erfinv_series<5>",
    "erfinv_series<10>",
    "erfinv_series<15>",
    "erfinv_series<20>",
    "erfinv_winitzki",
    "erfinv",
    # "erfinvLUT<100>",
]
for d, l in zip(data_list, label_list):
    plt.plot(d[:, 0], d[:, -1], label=l)

plt.legend()
plt.xlabel("x")
plt.ylabel("n_digits")
figname = "n_digits.png"
print(f"Saving figure to {figname}")
plt.savefig(figname, dpi=300)

plt.close()
for d, l in zip(data_list, label_list):
    plt.plot(d[:, 0], d[:, 1], label=l)
figname = "erfinv.png"
plt.legend()
plt.xlabel("x")
plt.ylabel("erfinv(x)")
print(f"Saving figure to {figname}")
plt.savefig(figname, dpi=300)
