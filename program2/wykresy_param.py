import numpy as np
import matplotlib.pyplot as plt

dane = []
with open('parametry.txt') as f:
    for line in f.readlines():
        dane.append(list(map(float, line.split())))

dane = np.array(dane)
t = dane[:, 0]
dane = dane[:, 1:]

plt.plot(t, dane)
plt.grid()
plt.legend(('N', 'x', 'E'))
plt.show()
