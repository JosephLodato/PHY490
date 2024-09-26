#include <iostream>
#include <complex>
#include <cmath>
#include <vector>
#include <valarray> // For FFT
#include <chrono>
#include <thread>

using namespace std;
using namespace std::this_thread;     // sleep_for, sleep_until
using namespace std::chrono_literals; // ns, us, ms, s, h, etc.
using std::chrono::system_clock;

#define M_PI 3.14159265358979323846

using Complex = std::complex<double>;
using CArray = std::valarray<Complex>;

// FFT function
void fft(CArray &x)
{
    const size_t N = x.size();
    if (N <= 1)
        return;

    // Divide
    CArray even = x[std::slice(0, N / 2, 2)];
    CArray odd = x[std::slice(1, N / 2, 2)];

    // Conquer
    fft(even);
    fft(odd);

    // Combine
    for (size_t k = 0; k < N / 2; ++k)
    {
        Complex t = std::polar(1.0, -2 * M_PI * k / N) * odd[k];
        x[k] = even[k] + t;
        x[k + N / 2] = even[k] - t;
    }
}

void applyWindowFunction(CArray &signal, double M)
{
    for (size_t i = 0; i < signal.size(); i++)
    {
        double t = i * (M / signal.size()); // Scale t across the window
        double windowFunction = 0.5 * (1 - cos(2 * M_PI * t / M));
        signal[i] *= windowFunction;
    }
}

int main(int argc, char const *argv[])
{
    const int frequencyViewStop = 20; // Just for adjusting the range of frequencies plotted
    // Set up input signal
    const double tStart = -2.5;
    const double tStop = 120.0; // NOTE: Increasing this increases the resolution of the graph
    const double dT = 0.001;

    const int N = static_cast<int>((tStop - tStart) / dT); // Total number of samples
    CArray signal(N);

    // Populate signal with the sampled input signal, apply cosine input
    for (int i = 0; i < N; i++)
    {
        double t = tStart + i * dT;
        signal[i] = cos(1 * (2 * M_PI) * t);
        signal[i] += cos(4 * (2 * M_PI) * t);
        signal[i] += sin(8 * (2 * M_PI) * t);
        signal[i] += cos(18 * (2 * M_PI) * t);

    }

    double M = tStop - tStart;
    applyWindowFunction(signal, M);

    // Perform FFT
    fft(signal);

    // Output results
    for (size_t i = 0; i < N / 2; i++) // Only positive frequencies
    {
        double freq = i / ((tStop - tStart)); // Calculate the frequency corresponding to this bin
        if (freq > frequencyViewStop)
        {
            break;
        }
        printf("%f %f\n", freq, std::abs(signal[i]));
    }

    return 0;
}