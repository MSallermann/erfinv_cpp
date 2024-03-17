import numpy as np
import matplotlib.pyplot as plt


data_series = np.loadtxt("series.txt")
data_winitzki = np.loadtxt("winitzki.txt")
data_full = np.loadtxt("full.txt")

data_list = [data_series, data_winitzki, data_full]
label_list = ["detail::erfinv_series<double,5>", "detail::erfinv_winitzki", "erfinv"]
for d, l in zip(data_list, label_list):
    plt.plot(d[:, 0], d[:, -1], label=l)

plt.legend()

plt.xlabel("x")
plt.ylabel("n_digits")

plt.savefig("n_digits.png", dpi=300)
