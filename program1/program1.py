import numpy as np
import matplotlib.pyplot as plt

save_txt_file = False
plot_results = False


def save_to_file(plik='wyniki.xyz'):  # Funkcja od zapisu kryształu do pliku
    with open(plik, 'w') as f:
        f.write(f'{N}\n\n')
        for indx, (x, y, z) in enumerate(r):
            f.write(f'atom{indx}\t{x}\t{y}\t{z}\n')


def plot_momentum():  # Funkcja do histogramów pędu
    fig, axs = plt.subplots(1, 3)
    for (indx, title), ax in zip(enumerate('x y z'.split()), axs):
        ax.hist(p[:, indx], bins=30)
        ax.set_title(f'Pęd w osi {title}')
        ax.grid()
    plt.show()


def load_initial_cond(plik='dane.txt'):  # Wczytanie warunków początkowych z pliku
    global n, m, a, T0, N, k
    '''
    n - liczba atomów w jednej osi,
    m - masa atomu,
    a - stała atomowa,
    T0 - temperatura początkowa,
    N - całkowita liczba atomów,
    k - stała Boltzmanna
    '''
    data = []
    with open(plik) as f:
        for line in f.readlines():
            for val in line.split():
                try:
                    data.append(float(val))
                except:
                    pass  # val nie jest liczbą - odrzucenie komentarzy

    data[0] = int(data[0])
    k = 8.31e-3
    n, m, a, T0 = data
    N = n**3


# Wczytanie war. początkowych
load_initial_cond()

# Wektory położenia i pędu
r = np.array([[0, 0, 0] for i in range(N)]).astype(np.float)  # położenia
p = np.array([[0, 0, 0] for i in range(N)]).astype(np.float)  # pędy

# wektory potrzebne do początkowych pozycji atomów
b0 = np.array([a, 0, 0])
b1 = np.array([a/2, a*np.sqrt(3)/2, 0])
b2 = np.array([a/2, a*np.sqrt(3)/6, a*np.sqrt(2/3)])

# Ustawienie atomów w siatce rombu
for i0 in range(n):
    for i1 in range(n):
        for i2 in range(n):
            i = i0 + i1*n + i2*n**2
            r[i, :] = (i0-(n-1)/2)*b0 + (i1-(n-1)/2)*b1 + (i2-(n-1)/2)*b2

if save_txt_file:
    save_to_file()

# Losowanie pędów
random_01 = np.random.random(size=(N, 3))
random_01[random_01 == 0] = 1
E_kin = np.array(-1/2 * k * T0 * np.log(random_01))
for indx, (x, y, z) in enumerate(np.random.randint(2, size=(N, 3)) * 2 - 1):
    p[indx, :] = x*np.sqrt(2*m*E_kin[indx, 0]), y*np.sqrt(2*m*E_kin[indx, 1]), z*np.sqrt(2*m*E_kin[indx, 2])

p = p - sum(p)/N  # wyeliminowanie ruchu środka masy

if plot_results:
    plot_momentum()
