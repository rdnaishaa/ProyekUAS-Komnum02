# PROYEK UAS KOMNUM-02

# Implementasi Metode Runge-Kutta Orde-4 untuk Model Lotka-Volterra

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![C++](https://img.shields.io/badge/C++-11-blue.svg)](https://isocpp.org/)
[![Python](https://img.shields.io/badge/Python-3.x-green.svg)](https://python.org/)

*Simulasi dinamika populasi predator-prey menggunakan metode numerik canggih*

---

## Deskripsi

Model Lotka-Volterra adalah sistem persamaan diferensial biasa nonlinear yang menggambarkan interaksi antara dua populasi dalam ekosistem: **predator** dan **prey** (mangsa). Model ini diformulasikan sebagai berikut:

$$\frac{dx}{dt} = ax - bxy$$

$$\frac{dy}{dt} = -cy + dxy$$

| Variabel | Deskripsi |
|----------|-----------|
| **x(t)** | Ukuran populasi prey pada waktu t |
| **y(t)** | Ukuran populasi predator pada waktu t |
| **a** | Laju pertumbuhan intrinsik prey (tanpa adanya predator) |
| **b** | Laju predasi per predator |
| **c** | Laju kematian intrinsik predator (tanpa adanya prey) |
| **d** | Efisiensi konversi prey menjadi predator |

Proyek ini mengimplementasikan metode Runge-Kutta orde-4 untuk menyelesaikan sistem persamaan diferensial ini dan menganalisis dinamika populasi yang dihasilkan.

---

## Fitur

| **Implementasi** | **Analisis** | **Output** |
|:---:|:---:|:---:|
| Metode Runge-Kutta orde-4 | Parameter model yang dapat disesuaikan | Output hasil simulasi dalam format CSV |
| Kondisi awal yang dapat disesuaikan | Analisis statistik hasil simulasi | Visualisasi dengan Python (opsional) |
| Sistem persamaan diferensial | Estimasi periode osilasi | Analisis nilai min/max/rata-rata/amplitudo |

---

## Cara Menggunakan

### Persyaratan

- Compiler C++ yang mendukung C++11 atau lebih tinggi (g++, clang, dll.)
- Python 3.x dengan matplotlib dan pandas (opsional, untuk visualisasi)

### Kompilasi dan Menjalankan

**1. Clone repositori ini**

```bash
git clone https://github.com/username/lotka-volterra-rk4.git
cd lotka-volterra-rk4
```

**2. Kompilasi program C++**

```bash
g++ -o lotka_volterra_rk4 lotka_volterra_rk4.cpp -std=c++11
```

**3. Jalankan program**

```bash
./lotka_volterra_rk4
```

**4. (Opsional) Visualisasi hasil dengan Python**

```bash
python visualize_results.py
```

---

### Menyesuaikan Parameter

Anda dapat menyesuaikan parameter model dan simulasi dengan mengedit nilai-nilai dalam fungsi `main()`:

**Parameter Model:**
```cpp
ModelParameters modelParams = {
    1.0,  // a: laju pertumbuhan prey
    0.5,  // b: efisiensi predasi
    0.75, // c: laju kematian predator
    0.25  // d: efisiensi konversi
};
```

**Parameter Simulasi:**
```cpp
SimulationParameters simParams = {
    10.0,           // x0: populasi awal prey
    5.0,            // y0: populasi awal predator
    0.01,           // h: ukuran langkah
    20.0,           // tMax: waktu simulasi maksimum
    "lotka_volterra_results.csv" // outputFile: nama file output
};
```

---

## 📊 Hasil

Simulasi menghasilkan dinamika osilasi periodik yang khas dari model Lotka-Volterra, dengan karakteristik sebagai berikut:

### Populasi Prey (x)

| Metrik | Nilai |
|--------|-------|
| **Nilai minimum** | ~0.174 |
| **Nilai maksimum** | ~13.158 |
| **Rata-rata** | ~3.024 |
| **Amplitudo** | ~6.492 |

### Populasi Predator (y)

| Metrik | Nilai |
|--------|-------|
| **Nilai minimum** | ~0.194 |
| **Nilai maksimum** | ~7.507 |
| **Rata-rata** | ~2.291 |
| **Amplitudo** | ~3.656 |

### Periode Osilasi
**Estimasi periode: ~9.39 unit waktu**

### 📈 Visualisasi

| Deskripsi | Gambar |
|-----------|---------|
| **Dinamika populasi prey-predator dan phase portrait** | ![Dinamika Populasi](https://i.imgur.com/GY4lfxi.png) |

---

## Metode Runge-Kutta Orde-4

Metode Runge-Kutta orde-4 (RK4) adalah metode numerik untuk menyelesaikan persamaan diferensial biasa. Untuk persamaan diferensial orde pertama:

$$\frac{dy}{dt} = f(t, y)$$

Metode RK4 menggunakan formula:

$$k_1 = h \cdot f(t_n, y_n)$$

$$k_2 = h \cdot f\left(t_n + \frac{h}{2}, y_n + \frac{k_1}{2}\right)$$

$$k_3 = h \cdot f\left(t_n + \frac{h}{2}, y_n + \frac{k_2}{2}\right)$$

$$k_4 = h \cdot f(t_n + h, y_n + k_3)$$

$$y_{n+1} = y_n + \frac{k_1 + 2k_2 + 2k_3 + k_4}{6}$$

**Komponen Formula:**

| Simbol | Deskripsi |
|--------|-----------|
| **h** | Ukuran langkah |
| **t_n** | Waktu pada langkah ke-n |
| **y_n** | Nilai variabel pada langkah ke-n |

Untuk sistem persamaan diferensial seperti model Lotka-Volterra, metode ini diterapkan untuk setiap variabel (x dan y) secara simultan.

---

## 📁 Struktur Kode

```
Struktur Proyek
├── lotka_volterra_rk4.cpp       # Program utama yang mengimplementasikan RK4
├── visualize_results.py         # Script Python untuk visualisasi hasil
├── lotka_volterra_results.csv   # File output yang berisi hasil simulasi
└── README.md                    # Dokumentasi proyek
```

### Struktur Program C++

Program C++ terdiri dari beberapa bagian utama:

**1. Struktur Data**
- **`ModelParameters`**: Menyimpan parameter model Lotka-Volterra (a, b, c, d)
- **`SimulationParameters`**: Menyimpan parameter simulasi (x₀, y₀, h, tMax, outputFile)

**2. Fungsi Laju Perubahan**
- **`dxdt()`**: Menghitung laju perubahan populasi prey
- **`dydt()`**: Menghitung laju perubahan populasi predator

**3. Implementasi Runge-Kutta**
- **`rungeKutta4()`**: Mengimplementasikan metode Runge-Kutta orde-4 untuk model Lotka-Volterra

**4. Analisis Hasil**
- **`analyzeResults()`**: Menganalisis hasil simulasi (nilai min/max, rata-rata, amplitudo, periode)

**5. Fungsi Main**
- Menginisialisasi parameter model dan simulasi
- Menjalankan simulasi dan menganalisis hasil

---

## 🔍 Analisis Teoretis

### Titik Kesetimbangan

Model Lotka-Volterra memiliki dua titik kesetimbangan:

1. **(0, 0)**: Titik kesetimbangan trivial (tidak stabil)
2. **(c/d, a/b)**: Titik kesetimbangan non-trivial

Untuk parameter default (a=1.0, b=0.5, c=0.75, d=0.25), titik kesetimbangan non-trivial adalah **(3.0, 2.0)**.

### Konservasi Energi

Model Lotka-Volterra klasik memiliki sifat konservatif, dengan kuantitas yang dikonservasi:

$$E = dx - c\ln(x) + by - a\ln(y)$$

**Nilai E tetap konstan** sepanjang orbit dalam ruang fase, yang menjelaskan mengapa solusi membentuk kurva tertutup.

---

## Tentang Proyek

**Dibuat oleh**: R. Aisha Syauqi Ramadhani (2306250554)  
**Mata Kuliah**: Komputasi Numerik  
**Universitas**: Universitas Indonesia  

---

## Lisensi

Proyek ini dilisensikan di bawah [Lisensi MIT](LICENSE).

---

## Referensi

1. Chapra, S. C., & Canale, R. P. (2015). *Numerical Methods for Engineers* (7th ed.). McGraw-Hill Education.

---
