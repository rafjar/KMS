#include <vector>
#include <cmath>
#include <fstream>
#include <string>
#include <iostream>

int main() {
    std::ifstream dane("dane.txt"); // Plik z danymi początkowymi
    std::ofstream gestosc_prawd("gestosc_prawd.dat"); // Otworzenie pliku do zapisu gęstości prawdopodobieństwa
    std::ofstream parametry("parametry.txt"); // Zapisanie parametrów układu

    double L = 1; // Długość odcinka
    int N; // Liczba punktów na siatce
    dane >> N; dane.ignore(100, '\n'); // Wczytanie N z pliku
    double delta_x = 1./N;
    double tau = 0; // Czas początkowy
    double delta_tau;
    dane >> delta_tau; dane.ignore(100, '\n'); // Wczytanie delta_tau z pliku
    int kroki;
    dane >> kroki; dane.ignore(100, '\n'); // Wczytanie liczby kroków z pliku

    // Stałe fizyczne
    const double e = -1.6e-19; // Ładunek elektronu
    const double m = 9.1e-31; // Masa elektronu
    const double h_bar = 6.58e-16; // h kreślone w eV*s

    // Inne śmieszne parametry
    int n;
    double K0;
    double Omega;
    dane >> n; dane.ignore(100, '\n'); // Wczytanie n z pliku
    dane >> K0; dane.ignore(100, '\n'); // Wczytanie K0 z pliku
    dane >> Omega; dane.ignore(100, '\n'); // Wczytanie Omega z pliku

    double kappa = K0*e*m*L*L*L / (h_bar*h_bar);
    double omega = Omega*m*L*L / h_bar;

    // Wczytanie co ile zapisywać dane
    int zapis_parametrow, zapis_gestosci_prawd;
    dane >> zapis_parametrow; dane.ignore(100, '\n');
    dane >> zapis_gestosci_prawd;

    std::vector<double> psi_R(N), psi_I(N); // Tablice funkcji falowych
    std::vector<double> H_R(N), H_I(N); // Talice hamiltonianów

    // Inicjalizacja funkcji falowych
    for(int k=0; k<N; ++k)
        psi_R[k] = std::sqrt(2) * std::sin(n * M_PI * k*delta_x);

    // Inicjalizacja Hamiltonianów
    for(int k=1; k<N-1; ++k) {
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
            for(int k2=1; k2<N-1; ++k2)
                H_R[k2] = -1./2 * ((psi_R[k2-1] - 2*psi_R[k2] + psi_R[k2+1]) / (delta_x*delta_x))
                        + kappa * (k2*delta_x - 1./2) * psi_R[k2] * std::sin(omega*tau);

            psi_I[k1] = psi_I[k1] - H_R[k1] * delta_tau; // Obliczenie psi_I w czasie tau + d_tau

            tau += delta_tau/2;
            // Obliczenie H_I w czasie tau + d_tau
            for(int k2=1; k2<N-1; ++k2)
                H_I[k2] = -1./2 * ((psi_I[k2-1] - 2*psi_I[k2] + psi_I[k2+1]) / (delta_x*delta_x))
                        + kappa * (k2*delta_x - 1./2) * psi_I[k2] * std::sin(omega*tau);

            psi_R[k1] = psi_R[k1] + H_I[k1] * delta_tau/2; // Obliczenie psi_R w czasue tau + d_tau
        }

        // Zapis gęstości prawdopodobieństwa położenia
        if(!(i%zapis_gestosci_prawd)) {
            for(int j=0; j<N; ++j)
                gestosc_prawd << psi_R[j]*psi_R[j] + psi_I[j]*psi_I[j] << '\t';
            gestosc_prawd << '\n';
        }

        // Zapis parametrów układu do pliku
        if(!(i%zapis_parametrow)) {
            double N_zapis = 0, x_zapis = 0, E_zapis = 0;
            for(int j=0; j<N; j++) {
                N_zapis += psi_R[j]*psi_R[j] + psi_I[j]*psi_I[j];
                x_zapis += j*delta_x*(psi_R[j]*psi_R[j] + psi_I[j]*psi_I[j]);
                E_zapis += psi_R[j]*H_R[j] + psi_I[j]*H_I[j];
            }
            N_zapis *= delta_x;
            x_zapis *= delta_x;
            E_zapis *= delta_x;
            parametry << N_zapis << '\t' << x_zapis << '\t' << E_zapis << '\n';
        }
    }

    parametry.close();
    gestosc_prawd.close();
}