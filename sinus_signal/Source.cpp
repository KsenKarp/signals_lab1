#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <fstream>

#define N_SECONDS 10


void fill_csv_sin(double F, double quantization_min, double quantization_max, double phi, double T) {
    uint8_t quantization_levels_num = 64; // количество уровней равномерного квантования
    int N_SAMPLES = N_SECONDS * F + 1;
    uint8_t* digital_signal = new uint8_t[N_SAMPLES]; // коды квантования цифрового сигнала
    double k = M_PI * 2 / T; //(M_PI*2/T) - коэффициент k в синусе
    for (int i = 0; i < N_SAMPLES; i++) digital_signal[i] = (round((sin(i * k/ F + phi) + 1) * (quantization_levels_num - 1))) /
        (quantization_max - quantization_min);
    const char csv_file_name[64] = "data.csv";
    std::ofstream csv_file;
    csv_file.open(csv_file_name);
    csv_file << "time,signal,sine\n";
    for (size_t i = 0; i < N_SAMPLES; ++i)
    {
        double signal_val = double(digital_signal[i]) * (quantization_max - quantization_min) / 
            (quantization_levels_num - 1) + quantization_min;

        csv_file << (i / F) << "," << signal_val << /*"," << sin(k * (i / F) + phi) << */"\n";
    }
    csv_file.close();
    delete digital_signal;
}

int main()  {

    double F = 10.; // частота дискретизации
    double quantization_min = -1., quantization_max = 1.;
    double phi = 2*M_PI;
    double T = 2 * M_PI;
    fill_csv_sin(F, quantization_min, quantization_max, phi, T);

    return 0;
}