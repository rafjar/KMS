import numpy as np
from numba import jit
from tqdm import tqdm


# Wczytanie warunków początkowych
def load_initial_cond(plik='dane.txt'):
    data = []

    with open(plik) as f:
        for line in f.readlines():
            for val in line.split():
                try:
                    data.append(float(val))
                except:
                    pass
    k = 8.31e-3
    n, m, a, T0, epsilon, L, f, R, tau, S_o, S_d, S_out, S_xyz = data
    n = int(n)
    S_o = int(S_o)
    S_d = int(S_d)
    S_out = int(S_out)
    S_xyz = int(S_xyz)
    N = n**3

    return k, n, N, m, a, T0, epsilon, L, f, R, tau, S_o, S_d, S_out, S_xyz


# Ustawienie atomów w siatce rombu
def set_positions(r_org, n, b0, b1, b2):
    r = np.zeros(r_org.shape)
    for i0 in range(n):
        for i1 in range(n):
            for i2 in range(n):
                i = i0 + i1*n + i2*n**2
                r[i, :] = (i0-(n-1)/2)*b0 + (i1-(n-1)/2)*b1 + (i2-(n-1)/2)*b2
    return r


# Losowanie pędów
def randomize_momentum(N, k, T0, m):
    p = np.array([[0, 0, 0] for i in range(N)]).astype(np.float)
    random_01 = np.random.random(size=(N, 3))
    random_01[random_01 == 0] = 1
    E_kin = np.array(-1/2 * k * T0 * np.log(random_01))
    for indx, (x, y, z) in enumerate(np.random.choice((-1, 1), size=(N, 3))):
        p[indx, :] = x*np.sqrt(2*m*E_kin[indx, 0]), y*np.sqrt(2*m*E_kin[indx, 1]), z*np.sqrt(2*m*E_kin[indx, 2])

    return p - sum(p)/N  # wyeliminowanie ruchu środka masy


# Liczenie sił i potencjałów
def count_FV(N, r, L, f, epsilon, R):
    Vs = np.zeros(N)
    Vp = np.zeros(N)
    Fs = np.zeros((N, 3))
    Fp = np.zeros((N, 3))

    for indx1, particle1 in enumerate(r):
        r_abs = np.linalg.norm(particle1)

        if r_abs > L:
            Vs[indx1] = f/2 * (r_abs - L)**2
            Fs[indx1] = f*(L-r_abs) * particle1/r_abs

        for indx2, particle2 in enumerate(r[:indx1]):
            r_abs = np.linalg.norm(particle1 - particle2)
            Vp[indx1] += epsilon * ((R/r_abs)**12 - 2*(R/r_abs)**6)
            F = 12*epsilon * ((R/r_abs)**12 - (R/r_abs)**6) * (particle1 - particle2)/r_abs**2
            Fp[indx1] += F
            Fp[indx2] -= F

    return np.sum(Vs) + np.sum(Vp), Fs + Fp, Fs


# Obliczenie ciśnienia
def count_pressure(Fs, L):
    P = np.linalg.norm(Fs, axis=-1) / (4 * np.pi * L**2)
    return P


# Obliczenie energii
def count_energy(p, m):
    E_kin = np.linalg.norm(p, axis=-1)**2 / (2*m)
    return E_kin


# Obliczenie temperatury
def count_temperature(p, m, N, k, E_kin):
    T = 2/(3*N*k) * np.sum(E_kin)
    return T


# Obliczenie Hamiltonianiu
def count_hamiltonian(p, m, V):
    H = np.sum(np.linalg.norm(p, axis=-1)**2 / (2*m)) + V
    return H


# Całkowanie równiania ruchu
def integrate(p, F, tau, m, N, r, L, f, epsilon, R, t):
    p += F*tau/2
    r += p*tau/m
    V, F, Fs = count_FV(N, r, L, f, epsilon, R)
    P = count_pressure(Fs, L)
    p += F*tau/2
    t += tau

    return t, V, F, P, r, p


# Zapisanie parametrów układu
def save_properties(properties_file, t, H, V, T, P):
    properties_file.write(f'{t}\t{H}\t{V}\t{T}\t{P}\n')


# Zapisanie położeń atomów
def save_positions(position_file, N, r, E_kin):
    position_file.write(f'{N}\n\n')
    for (x, y, z), E in zip(r, E_kin):
        position_file.write(f'Ar\t{x}\t{y}\t{z}\t{E}\n')


# Symulacja całego układu
def symulacja(S_o, S_d, p, F, tau, m, N, r, L, f, epsilon, R, t, k, S_out, S_xyz, properties_file, position_file):
    T_avg = 0.
    P_avg = 0.
    H_avg = 0.

    for i in tqdm(range(S_o + S_d)):
        t, V, F, P, r, p = integrate(p, F, tau, m, N, r, L, f, epsilon, R, t)
        E_kin = count_energy(p, m)
        T = count_temperature(p, m, N, k, E_kin)
        H = count_hamiltonian(p, m, V)

        if not i % S_out:
            save_properties(properties_file, t, H, V, T, P)

        if not i % S_xyz:
            save_positions(position_file, N, r, E_kin)

        if i >= S_o:
            T_avg += T
            P_avg += P
            H_avg += H

    T_avg /= S_d
    P_avg /= S_d
    H_avg /= S_d

    return T_avg, P_avg, H_avg


def main():
    position_file = open('avs.txt', 'w')
    properties_file = open('properties.txt', 'w')

    # Wczytanie warunków początkowych
    k, n, N, m, a, T0, epsilon, L, f, R, tau, S_o, S_d, S_out, S_xyz = load_initial_cond()

    # Wektory położenia i pędu
    r = np.array([[0, 0, 0] for i in range(N)]).astype(np.float)  # położenia
    p = np.array([[0, 0, 0] for i in range(N)]).astype(np.float)  # pędy

    # wektory potrzebne do początkowych pozycji atomów
    b0 = np.array([a, 0, 0])
    b1 = np.array([a/2, a*np.sqrt(3)/2, 0])
    b2 = np.array([a/2, a*np.sqrt(3)/6, a*np.sqrt(2/3)])

    # Ustawienie położeń początkowych
    r = set_positions(r, n, b0, b1, b2)
    t = 0.

    # Losowanie pędów
    p = randomize_momentum(N, k, T0, m)

    V, F, Fs = count_FV(N, r, L, f, epsilon, R)

    T_avg, P_avg, H_avg = symulacja(S_o, S_d, p, F, tau, m, N, r, L, f, epsilon, R, t, k, S_out, S_xyz, properties_file, position_file)

    position_file.close()
    properties_file.close()


if __name__ == '__main__':
    main()
