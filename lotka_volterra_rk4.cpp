#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>

/**
 * @file lotka_volterra_rk4.cpp
 * @brief Implementasi metode Runge-Kutta orde-4 untuk model Lotka-Volterra
 * @author 2306250554_R.AishaSyauqiRamadhani
 * 
 * Program ini menyelesaikan sistem persamaan diferensial model Lotka-Volterra
 * menggunakan metode Runge-Kutta orde-4. Model ini menggambarkan interaksi
 * antara populasi predator dan prey dalam suatu ekosistem.
 */

// Struktur untuk menyimpan parameter model Lotka-Volterra
struct ModelParameters {
    double a; // Laju pertumbuhan prey
    double b; // Efisiensi predasi
    double c; // Laju kematian predator
    double d; // Efisiensi konversi prey menjadi predator
};

// Struktur untuk menyimpan kondisi simulasi
struct SimulationParameters {
    double x0;    // Populasi awal prey
    double y0;    // Populasi awal predator
    double h;     // Ukuran langkah
    double tMax;  // Waktu simulasi maksimum
    std::string outputFile; // Nama file output
};

/**
 * @brief Fungsi laju perubahan populasi prey
 * 
 * @param t Waktu
 * @param x Populasi prey
 * @param y Populasi predator
 * @param params Parameter model
 * @return double Laju perubahan populasi prey (dx/dt)
 */
double dxdt(double t, double x, double y, const ModelParameters& params) {
    return params.a * x - params.b * x * y;
}

/**
 * @brief Fungsi laju perubahan populasi predator
 * 
 * @param t Waktu
 * @param x Populasi prey
 * @param y Populasi predator
 * @param params Parameter model
 * @return double Laju perubahan populasi predator (dy/dt)
 */
double dydt(double t, double x, double y, const ModelParameters& params) {
    return -params.c * y + params.d * x * y;
}

/**
 * @brief Implementasi metode Runge-Kutta orde-4 untuk model Lotka-Volterra
 * 
 * @param modelParams Parameter model Lotka-Volterra
 * @param simParams Parameter simulasi
 */
void rungeKutta4(const ModelParameters& modelParams, const SimulationParameters& simParams) {
    // Inisialisasi variabel
    double t = 0.0;
    double x = simParams.x0;
    double y = simParams.y0;
    
    // Membuka file output
    std::ofstream outFile(simParams.outputFile);
    if (!outFile.is_open()) {
        std::cerr << "Error: Tidak dapat membuka file output!" << std::endl;
        return;
    }
    
    // Menulis header file
    outFile << "t,x,y" << std::endl;
    
    // Menulis kondisi awal
    outFile << std::fixed << std::setprecision(6);
    outFile << t << "," << x << "," << y << std::endl;
    
    // Iterasi metode Runge-Kutta
    while (t < simParams.tMax) {
        // Menghitung k1
        double k1x = simParams.h * dxdt(t, x, y, modelParams);
        double k1y = simParams.h * dydt(t, x, y, modelParams);
        
        // Menghitung k2
        double k2x = simParams.h * dxdt(t + simParams.h/2, x + k1x/2, y + k1y/2, modelParams);
        double k2y = simParams.h * dydt(t + simParams.h/2, x + k1x/2, y + k1y/2, modelParams);
        
        // Menghitung k3
        double k3x = simParams.h * dxdt(t + simParams.h/2, x + k2x/2, y + k2y/2, modelParams);
        double k3y = simParams.h * dydt(t + simParams.h/2, x + k2x/2, y + k2y/2, modelParams);
        
        // Menghitung k4
        double k4x = simParams.h * dxdt(t + simParams.h, x + k3x, y + k3y, modelParams);
        double k4y = simParams.h * dydt(t + simParams.h, x + k3x, y + k3y, modelParams);
        
        // Update nilai x dan y
        x = x + (k1x + 2*k2x + 2*k3x + k4x) / 6.0;
        y = y + (k1y + 2*k2y + 2*k3y + k4y) / 6.0;
        
        // Update waktu
        t += simParams.h;
        
        // Menulis hasil ke file
        outFile << t << "," << x << "," << y << std::endl;
    }
    
    // Menutup file
    outFile.close();
    
    std::cout << "Simulasi selesai. Hasil disimpan di " << simParams.outputFile << std::endl;
}

/**
 * @brief Fungsi untuk menghitung karakteristik dinamika populasi
 * 
 * @param filename Nama file hasil simulasi
 */
void analyzeResults(const std::string& filename) {
    std::ifstream inFile(filename);
    if (!inFile.is_open()) {
        std::cerr << "Error: Tidak dapat membuka file hasil!" << std::endl;
        return;
    }
    
    std::string header;
    std::getline(inFile, header); // Membaca header
    
    std::vector<double> times, preyPop, predatorPop;
    double t, x, y;
    char comma;
    
    // Membaca data dari file
    while (inFile >> t >> comma >> x >> comma >> y) {
        times.push_back(t);
        preyPop.push_back(x);
        predatorPop.push_back(y);
    }
    
    inFile.close();
    
    // Mencari nilai maksimum dan minimum
    double maxPrey = preyPop[0], minPrey = preyPop[0];
    double maxPred = predatorPop[0], minPred = predatorPop[0];
    
    for (size_t i = 1; i < preyPop.size(); ++i) {
        if (preyPop[i] > maxPrey) maxPrey = preyPop[i];
        if (preyPop[i] < minPrey) minPrey = preyPop[i];
        
        if (predatorPop[i] > maxPred) maxPred = predatorPop[i];
        if (predatorPop[i] < minPred) minPred = predatorPop[i];
    }
    
    // Menghitung rata-rata
    double avgPrey = 0.0, avgPred = 0.0;
    for (size_t i = 0; i < preyPop.size(); ++i) {
        avgPrey += preyPop[i];
        avgPred += predatorPop[i];
    }
    avgPrey /= preyPop.size();
    avgPred /= predatorPop.size();
    
    // Menampilkan hasil analisis
    std::cout << "\nAnalisis Hasil Simulasi:" << std::endl;
    std::cout << "------------------------" << std::endl;
    std::cout << "Populasi Prey:" << std::endl;
    std::cout << "  Minimum: " << minPrey << std::endl;
    std::cout << "  Maksimum: " << maxPrey << std::endl;
    std::cout << "  Rata-rata: " << avgPrey << std::endl;
    std::cout << "  Amplitudo: " << (maxPrey - minPrey) / 2.0 << std::endl;
    
    std::cout << "\nPopulasi Predator:" << std::endl;
    std::cout << "  Minimum: " << minPred << std::endl;
    std::cout << "  Maksimum: " << maxPred << std::endl;
    std::cout << "  Rata-rata: " << avgPred << std::endl;
    std::cout << "  Amplitudo: " << (maxPred - minPred) / 2.0 << std::endl;
    
    // Estimasi periode osilasi (metode sederhana)
    // Mencari dua puncak berturut-turut untuk prey
    double period = 0.0;
    bool foundPeak = false;
    size_t peakIndex = 0;
    
    for (size_t i = 1; i < preyPop.size() - 1; ++i) {
        if (preyPop[i] > preyPop[i-1] && preyPop[i] > preyPop[i+1]) {
            if (!foundPeak) {
                foundPeak = true;
                peakIndex = i;
            } else {
                period = times[i] - times[peakIndex];
                break;
            }
        }
    }
    
    if (period > 0.0) {
        std::cout << "\nPeriode Osilasi (estimasi): " << period << std::endl;
    } else {
        std::cout << "\nPeriode Osilasi: Tidak dapat diestimasi dari data" << std::endl;
    }
}

int main() {
    // Mendefinisikan parameter model
    ModelParameters modelParams = {
        1.0,  // a: laju pertumbuhan prey
        0.5,  // b: efisiensi predasi
        0.75, // c: laju kematian predator
        0.25  // d: efisiensi konversi
    };
    
    // Mendefinisikan parameter simulasi
    SimulationParameters simParams = {
        10.0,           // x0: populasi awal prey
        5.0,            // y0: populasi awal predator
        0.01,           // h: ukuran langkah
        20.0,           // tMax: waktu simulasi maksimum
        "lotka_volterra_results.csv" // outputFile: nama file output
    };
    
    // Menampilkan informasi simulasi
    std::cout << "Simulasi Model Lotka-Volterra dengan Metode Runge-Kutta Orde-4" << std::endl;
    std::cout << "=============================================================" << std::endl;
    std::cout << "Parameter Model:" << std::endl;
    std::cout << "  a (laju pertumbuhan prey): " << modelParams.a << std::endl;
    std::cout << "  b (efisiensi predasi): " << modelParams.b << std::endl;
    std::cout << "  c (laju kematian predator): " << modelParams.c << std::endl;
    std::cout << "  d (efisiensi konversi): " << modelParams.d << std::endl;
    
    std::cout << "\nParameter Simulasi:" << std::endl;
    std::cout << "  Populasi awal prey (x0): " << simParams.x0 << std::endl;
    std::cout << "  Populasi awal predator (y0): " << simParams.y0 << std::endl;
    std::cout << "  Ukuran langkah (h): " << simParams.h << std::endl;
    std::cout << "  Waktu simulasi maksimum: " << simParams.tMax << std::endl;
    std::cout << "  File output: " << simParams.outputFile << std::endl;
    
    // Jalankan simulasi
    std::cout << "\nMenjalankan simulasi..." << std::endl;
    rungeKutta4(modelParams, simParams);
    
    // Analisis hasil
    analyzeResults(simParams.outputFile);
    
    return 0;
}