#include <vector>
#include <cmath>
#include <fstream>
#include <string>
#include <iostream>

inline double count_H (const std::vector<double> &psi, const int &k, const double &delta_x, const double &kappa, const double &omega, const double &tau) {
    if(k == 0 || k == (int)psi.size()-1)
        return 0;
    return -0.5 * ((psi[k-1] - 2*psi[k] + psi[k+1]) / (delta_x*delta_x)) 
           + kappa * (k*delta_x - 0.5) * psi[k] * std::sin(omega * tau);
}

int main() {
    std::ifstream dane("dane.txt"); // Plik z danymi początkowymi
    std::ofstream gestosc_prawd("gestosc_prawd.dat"); // Otworzenie pliku do zapisu gęstości prawdopodobieństwa
    std::ofstream parametry("parametry.txt"); // Zapisanie parametrów układu

    int N; // Liczba punktów na siatce
    dane >> N; dane.ignore(100, '\n'); // Wczytanie N z pliku
    double delta_x = 1./N;
    double tau = 0; // Czas początkowy
    double delta_tau;
    dane >> delta_tau; dane.ignore(100, '\n'); // Wczytanie delta_tau z pliku
    int kroki;
    dane >> kroki; dane.ignore(100, '\n'); // Wczytanie liczby kroków z pliku

    // Inne śmieszne parametry
    int n;
    double kappa;
    double omega;
    dane >> n; dane.ignore(100, '\n'); // Wczytanie n z pliku
    dane >> kappa; dane.ignore(100, '\n'); // Wczytanie K0 z pliku
    dane >> omega; dane.ignore(100, '\n'); // Wczytanie Omega z pliku

    // Wczytanie co ile zapisywać dane
    int zapis_parametrow, zapis_gestosci_prawd;
    dane >> zapis_parametrow; dane.ignore(100, '\n');
    dane >> zapis_gestosci_prawd;

    std::vector<double> psi_R(N+1), psi_I(N+1); // Tablice funkcji falowych
    std::vector<double> H_R(N+1), H_I(N+1); // Talice hamiltonianów

    // Inicjalizacja funkcji falowych
    for(int k=0; k<=N; ++k)
        psi_R[k] = std::sqrt(2) * std::sin(n * M_PI * k*delta_x);

    // Inicjalizacja Hamiltonianów
    for(int k=1; k<=N-1; ++k) {
        H_R[k] = count_H(psi_R, k, delta_x, kappa, omega, tau);

        H_I[k] = count_H(psi_I, k, delta_x, kappa, omega, tau);
    }

    // Pętla symulacji
    for(int i=0; i<kroki; ++i) {
        for(int k=0; k<=N; ++k) {
            // Obliczenie psi_R w czasie tau + d_tau/2
            psi_R[k] += H_I[k] * delta_tau/2;

            tau += delta_tau/2; // Czas tau + d_tau/2
            // Obliczenie H_R w czasie tau + d_tau/2
            H_R[k] = count_H(psi_R, k, delta_x, kappa, omega, tau);

            // Obliczenie psi_I w czasie tau + d_tau
            psi_I[k] -= H_R[k] * delta_tau;

            tau += delta_tau/2; // Czas tau + d_tau
            // Obliczenie H_I w czasie tau + d_tau
            H_I[k] = count_H(psi_I, k, delta_x, kappa, omega, tau);

            // Obliczenie psi_R w czasie tau + d_tau
            psi_R[k] += H_I[k] * delta_tau/2;
        }

        // Zapis gęstości prawdopodobieństwa położenia
        if(!(i%zapis_gestosci_prawd)) {
            for(int j=0; j<N; ++j)
                gestosc_prawd << psi_R[j]*psi_R[j] + psi_I[j]*psi_I[j] << '\t';
            gestosc_prawd << '\n';
        }

        // Zapis parametrów układu do pliku
        if(!(i%zapis_parametrow)) {
            // Obliczenie H_R w czasie tau
            for(int k=1; k<=N-1; ++k) 
                H_R[k] = count_H(psi_R, k, delta_x, kappa, omega, tau);

            double N_zapis = 0, x_zapis = 0, E_zapis = 0;
            for(int k=0; k<=N; ++k) {
                N_zapis += psi_R[k]*psi_R[k] + psi_I[k]*psi_I[k];
                x_zapis += k*delta_x*(psi_R[k]*psi_R[k] + psi_I[k]*psi_I[k]);
                E_zapis += psi_R[k]*H_R[k] + psi_I[k]*H_I[k];
            }
            N_zapis *= delta_x;
            x_zapis *= delta_x;
            E_zapis *= delta_x;
            parametry << tau << '\t' << N_zapis << '\t' << x_zapis << '\t' << E_zapis << '\n';
        }
    }

    parametry.close();
    gestosc_prawd.close();
}