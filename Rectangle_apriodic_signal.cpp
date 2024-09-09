#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <ctime>

#define N_SECONDS 100

double rectangle_impulse(double x, double t, double x_m) {
    srand(time(0));
    double T1 = (double)rand() / RAND_MAX * (N_SECONDS - 2*t);
    double T2 = T1 + t + 0.01 + (double)rand() / RAND_MAX * (N_SECONDS - T1 - 2*t - 0.01);
    if ((x > T1 && x < T1 + t) || (x > T2 && x < T2 + t)) return x_m;
    else return 0;
}

void fill_csv_rectangle(double F, double quantization_max, double t) {
    int N_SAMPLES = N_SECONDS * F + 1;
    double quantization_min = 0.;
    uint8_t quantization_levels_num = 64;
    uint8_t* digital_signal = new uint8_t[N_SAMPLES];
    for (int i = 0; i < N_SAMPLES; i++) digital_signal[i] = (round(rectangle_impulse(i / F, t, quantization_max) * (quantization_levels_num - 1))) /
        (quantization_max - quantization_min);
    const char csv_file_name[64] = "data.csv";
    std::ofstream csv_file;
    csv_file.open(csv_file_name);
    csv_file << "time,signal,rectangle impulse\n";
    for (size_t i = 0; i < N_SAMPLES; ++i)
    {
        double signal_val = double(digital_signal[i]) * (quantization_max - quantization_min) /
            (quantization_levels_num - 1) + quantization_min;

        csv_file << (i / F) << "," << signal_val << "," << rectangle_impulse(i / F, t, quantization_max) << "\n";
    }
    csv_file.close();
}

int main() {

    double F = 40.; // частота дискретизации
    double quantization_max = 1.;
    double t = 4.;
    fill_csv_rectangle(F, quantization_max, t);

    return 0;
}