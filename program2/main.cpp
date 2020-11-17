#include <iostream>
#include <vector>
#include <cmath>

int main() {
    double L = 1; // Długość odcinka
    int N = 1001; // Liczba punktów na siatce
    double delta_x = 1./N;
    double tau = 0; // Czas początkowy
    double delta_tau = 1e-4;
    int kroki = 4000;

    // Stałe fizyczne
    const double e = -1.6e-19; // Ładunek elektronu
    const double m = 9.1e-31; // Masa elektronu
    const double h_bar = 6.58e-16; // h kreślone w eV*s

    int n = 1;
    double K0 = 0;
    double Omega = 0;

    double kappa = K0*e*m*L*L*L / (h_bar*h_bar);
    double omega = Omega*m*L*L / h_bar;

    std::vector<double> psi_R(N), psi_I(N); // Tablice funkcji falowych
    std::vector<double> H_R(N), H_I(N); // Talice hamiltonianów

    // Inicjalizacja funkcji falowych
    for(int k=0; k<N; k++)
        psi_R[k] = std::sqrt(2) * std::sin(n * M_PI * k*delta_x);

    // Inicjalizacja Hamiltonianów
    for(int k=1; k<N-1; k++) {
        H_R[k] = -1./2 * ((psi_R[k-1] - 2*psi_R[k] + psi_R[k+1]) / (delta_x*delta_x))
                + kappa * (k*delta_x - 1./2) * psi_R[k] * std::sin(omega*tau);

        H_I[k] = -1./2 * ((psi_I[k-1] - 2*psi_I[k] + psi_I[k+1]) / (delta_x*delta_x))
                + kappa * (k*delta_x - 1./2) * psi_I[k] * std::sin(omega*tau);
    }

    // Pętla symulacji
    for(int i=0; i<kroki; ++i) {
        for(int k1=0; k1<N; ++k1) {
            psi_R[k1] = psi_R[k1] + H_I[k1] * delta_tau/2; // Obliczenie psi_R w czasie tau + d_tau/2

            tau += delta_tau/2;
            // Obliczenie H_R w czasie tau + d_tau/2
            for(int k2=1; k2<N-1; k2++)
                H_R[k2] = -1./2 * ((psi_R[k2-1] - 2*psi_R[k2] + psi_R[k2+1]) / (delta_x*delta_x))
                        + kappa * (k2*delta_x - 1./2) * psi_R[k2] * std::sin(omega*tau);

            psi_I[k1] = psi_I[k1] - H_R[k1] * delta_tau; // Obliczenie psi_I w czasie tau + d_tau

            tau += delta_tau/2;
            // Obliczenie H_I w czasie tau + d_tau
            for(int k2=1; k2<N-1; k2++)
                H_I[k2] = -1./2 * ((psi_I[k2-1] - 2*psi_I[k2] + psi_I[k2+1]) / (delta_x*delta_x))
                        + kappa * (k2*delta_x - 1./2) * psi_I[k2] * std::sin(omega*tau);

            psi_R[k1] = psi_R[k1] + H_I[k1] * delta_tau/2; // Obliczenie psi_R w czasue tau + d_tau
        }
    }
    
}