import numpy as np
import matplotlib.pyplot as plt

dane = []
with open('gestosc_prawd.dat') as f:
    for line in f.readlines():
        dane.append(list(map(float, line.split())))

dane = np.array(dane)
L = np.linspace(0, 1, dane.shape[1])

plt.plot(L, dane[5])
plt.grid()
plt.show()
