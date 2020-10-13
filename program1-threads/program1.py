import numpy as np
import multiprocessing
import time
np.set_printoptions(linewidth=np.inf, precision=3)

save_txt_file = False
plot_results = False

position_file = open('avs.txt', 'w')
properties_file = open('properties.txt', 'w')


# Funkcja od zapisu kryształu do pliku
def save_to_file(plik='wyniki.xyz'):
    with open(plik, 'w') as f:
        f.write(f'{N}\n\n')
        for indx, (x, y, z) in enumerate(r):
            f.write(f'atom{indx}\t{x}\t{y}\t{z}\n')


# Funkcja do histogramów pędu
def plot_momentum():
    import matplotlib.pyplot as plt

    fig, axs = plt.subplots(1, 3)
    for (indx, title), ax in zip(enumerate('x y z'.split()), axs):
        ax.hist(p[:, indx], bins=30)
        ax.set_title(f'Pęd w osi {title}')
        ax.grid()
    plt.show()


# Wczytanie warunków początkowych z pliku
def load_initial_cond(plik='dane.txt'):
    global n, m, a, T0, N, k, epsilon, L, f, R, tau, S_o, S_d, S_out, S_xyz
    '''
    n - liczba atomów w jednej osi,
    m - masa atomu,
    a - stała atomowa,
    T0 - temperatura początkowa,
    N - całkowita liczba atomów,
    k - stała Boltzmanna,
    epsilon - minimum potencjału,
    L - promień kuli w której znajduje się gaz,
    f - stała sprężystości,
    R - odległość międzyatomowa dla której minimum potencjału,
    '''
    data = []
    with open(plik) as f:
        for line in f.readlines():
            for val in line.split():
                try:
                    data.append(float(val))
                except:
                    pass  # val nie jest liczbą - odrzucenie komentarzy

    k = 8.31e-3
    n, m, a, T0, epsilon, L, f, R, tau, S_o, S_d, S_out, S_xyz = data
    n = int(n)
    S_o = int(S_o)
    S_d = int(S_d)
    S_out = int(S_out)
    S_xyz = int(S_xyz)
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
for indx, (x, y, z) in enumerate(np.random.choice((-1, 1), size=(N, 3))):
    p[indx, :] = x*np.sqrt(2*m*E_kin[indx, 0]), y*np.sqrt(2*m*E_kin[indx, 1]), z*np.sqrt(2*m*E_kin[indx, 2])

p = p - sum(p)/N  # wyeliminowanie ruchu środka masy

if plot_results:
    plot_momentum()

# Obliczenie potencjału i sił
Vs = np.zeros(N)
Vp = np.zeros((N, N))
Fs = np.zeros((N, 3))
Fp = np.zeros((N, N, 3))
shm1 = multiprocessing.RawArray('d', N)
shm2 = multiprocessing.RawArray('d', N*N)
shm3 = multiprocessing.RawArray('d', N*3)
shm4 = multiprocessing.RawArray('d', N*N*3)
Vs = np.frombuffer(shm1, dtype=np.float).reshape(Vs.shape)
Vp = np.frombuffer(shm2, dtype=np.float).reshape(Vp.shape)
Fs = np.frombuffer(shm3, dtype=np.float).reshape(Fs.shape)
Fp = np.frombuffer(shm4, dtype=np.float).reshape(Fp.shape)
V = None
F = None
P = None
T = None
H = None
T_avg = 0
P_avg = 0
H_avg = 0
t = 0
i = [int(i*N/multiprocessing.cpu_count()) for i in range(multiprocessing.cpu_count())]
j = [int((i+1)*N/multiprocessing.cpu_count()) for i in range(multiprocessing.cpu_count())]


def task(i_min, i_max, N, r, f, L, R, Vs, Vp, Fs, Fp):  # funkcja wątku
    for indx1 in range(i_min, min(i_max, N)):
        particle1 = r[indx1]
        r_abs = np.linalg.norm(particle1)
        if r_abs >= L:
            Vs[indx1] = 1/2 * f * (r_abs-L)**2
            Fs[indx1] = f*(L-r_abs)*particle1/r_abs

        for indx2, particle2 in enumerate(r[:indx1]):
            r_abs = np.linalg.norm(particle1 - particle2)
            Vp[indx1, indx2] = epsilon * ((R / r_abs)**12 - 2*(R/r_abs)**6)
            Fp[indx1, indx2] = 12*epsilon * ((R/r_abs)**12 - (R/r_abs)**6) * ((particle1-particle2)/r_abs**2)
            Fp[indx2, indx1] = -Fp[indx1, indx2]


# Obliczenie sił i potencjałów
def count_forces():
    global Vs, Vp, Fs, Fp, V, F
    Fp[:] = 0

    processes = [multiprocessing.Process(target=task, args=(i[ii], j[ii], N, r, f, L, R, Vs, Vp, Fs, Fp)) for ii in range(multiprocessing.cpu_count())]
    for process in processes:
        process.start()
    for process in processes:
        process.join()

    V = np.sum(Vs) + np.sum(Vp)
    F = Fs + np.sum(Fp, axis=1)


# Obliczanie ciśnienia
def count_pressure():
    global Fs, L, P
    P = 1/(4*np.pi*L**2) * np.sum(np.linalg.norm(Fs, axis=1))


# Obliczenie energii
def count_energy():
    global E_kin
    E_kin = np.linalg.norm(p, axis=1)**2 / (2*m)


# Obliczanie temperatury
def count_temperature():
    global T, p, N, k, m
    count_energy()
    T = 2/(3*N*k) * np.sum(E_kin)


# Obliczanie Hamiltonianu
def count_hamiltonian():
    global p, V, H
    H = np.sum(np.linalg.norm(p, axis=1)**2 / (2*m)) + V


# Całkowanie równania ruchu
def integrate():
    global p, F, r, tau, t
    p = p + F*tau/2
    r = r + p*tau/m
    count_forces()
    count_pressure()
    p = p + F*tau/2
    t += tau


# Zapisanie parametrów układu
def save_properties():
    global H, V, T, P, t
    properties_file.write(f'{t}\t{H}\t{V}\t{T}\t{P}\n')


# Zapisanie położeń atomów
def save_positions():
    global r, E_kin, N
    position_file.write(f'{N}\n\n')
    for (x, y, z), E in zip(r, E_kin):
        position_file.write(f'Ar\t{x}\t{y}\t{z}\t{E}\n')


# Symulacja całego układu
def symulacja():
    global S_o, S_d, S_out, S_xyz, T_avg, H_avg, P_avg

    for i in range(S_o + S_d):
        integrate()
        count_temperature()
        count_hamiltonian()

        if not i % S_out:
            save_properties()

        if not i % S_xyz:
            save_positions()

        if i >= S_o:
            T_avg += T
            P_avg += P
            H_avg += H

    T_avg /= S_d
    P_avg /= S_d
    H_avg /= S_d


# Wstępne obliczenie sił i uruchomienie symulacji
count_forces()
start = time.time()
symulacja()
end = time.time()

print(f'Elapse time: {end-start}s.')

position_file.close()
properties_file.close()
