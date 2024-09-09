#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <fstream>

#define N_SECONDS 10

double triangle_impulse(double x, double T, double t, double x_m, double phi) {
    x = x - phi;
    if (x < 0) {
        while (x < 0) x = x + T;
    }
    while (x > T) x = x - T;
    if (x <= t / 2) return 2 * x_m * x / t;
    else if (t / 2 < x && x < t) return 2 * x_m - 2 * x_m * x / t;
    else return 0;
}

void fill_csv_triangle(double F, double quantization_max, double T, double t, double phi) {
    int N_SAMPLES = N_SECONDS * F;
    double quantization_min = 0.;
    uint8_t quantization_levels_num = 64;
    uint8_t* digital_signal = new uint8_t[N_SAMPLES];
    for (int i = 0; i < N_SAMPLES; i++) digital_signal[i] = (round(triangle_impulse(i / F, T, t, quantization_max, phi) * (quantization_levels_num - 1))) /
        (quantization_max - quantization_min);
    const char csv_file_name[64] = "data.csv";
    std::ofstream csv_file;
    csv_file.open(csv_file_name);
    csv_file << "time,signal,triangle impulse\n";
    for (size_t i = 0; i < N_SAMPLES; ++i)
    {
        double signal_val = double(digital_signal[i]) * (quantization_max - quantization_min) /
            (quantization_levels_num - 1) + quantization_min;

        csv_file << (i / F) << "," << signal_val << "," << triangle_impulse(i / F, T, t, quantization_max, phi) << "\n";
    }
    csv_file.close();
}

int main() {

    double F = 10.;
    double quantization_max = 2.5;
    double T = 6.;
    double t = 4.;
    double phi = 0.;
    fill_csv_triangle(F, quantization_max, T, t, phi);

    return 0;
}