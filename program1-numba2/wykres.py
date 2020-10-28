import matplotlib.pyplot as plt
import numpy as np

data = [[], [], [], [], []]
with open('properties.txt') as f:
    for line in f.readlines():
        values = list(map(float, line.split()))

        data[0].append(values[0])
        data[1].append(values[1])
        data[2].append(values[2])
        data[3].append(values[3])
        data[4].append(values[4])

    data = np.array(data)

plt.plot(data[0], data[1:].T)
plt.legend(('H', 'V', 'T', 'P'))
plt.grid()
plt.show()
