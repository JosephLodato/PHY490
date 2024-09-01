// #include <cmath.h>
#include <iostream>
#include <stdio.h>   /* Standard Library of Input and Output */
#include <complex.h> /* Standard Library of Complex Numbers */
#include <complex>
#include <cmath>
#include <vector>
#include <array>
#include <chrono>
#include <thread>

using namespace std;
using namespace std::this_thread;     // sleep_for, sleep_until
using namespace std::chrono_literals; // ns, us, ms, s, h, etc.
using std::chrono::system_clock;

#define M_PI 3.14159265358979323846

/* Inverse Fourier transform
 *  f(t) = 1/(2π) ∫F(ω)e^(iωt) dω
 *
 * f(t): original function in the time domain.
 * F(ω): fouier transform of f(t) in frequency domain
 * w: Angular Frequency
 * 1/(2π) norm factor
 *
 * Forward Fourier transform
 *  F(ω) = ∫f(t)e^(-iωt) dt
 *
 */

std::complex<double>
complexExponential(double f, double t)
{
    std::complex<double> exponent(0, -2 * M_PI * f * t); // -2πi * f * t
    return std::exp(exponent);
}

array<double, 2> finalCalculation(double f, double t)
{
    std::complex<double> result = (cos(6 * M_PI * t) + 1) * complexExponential(f, t);

    array<double, 2> temp;
    temp[0] = result.real();
    temp[1] = result.imag();

    return temp;
}

void DFT(double f, double t, double tStop, double dT)
{
    vector<array<double, 2>> points;

    while (t < tStop)
    {
        t += dT;
        points.push_back(finalCalculation(f, t));
    }

    // display all points
    for (int i = 0; i < points.size(); i++)
    {
        printf("c %f %f 0.01\n", points[i][0], points[i][1]);
    }
    printf("l -2 0 2 0\n");
    printf("l 0 -2 0 2\n");
    printf("t -2 2 \n Wrapping Frequency: %f \n", f);

    printf("F\n");
    points.clear();
}
int main()
{
    vector<array<double, 2>> points;

    double f = 1;      // Frequency
    double t = 1.0;    // Time
    double dT = 0.001; // Time Step
    double dF = 0.001; // Frequency Step
    
    double tStop = 5; // Max Time simulated for each f
    

    // calculate points across time
    while (f <= 10)
    {
        f += dF;
        DFT(f, t, tStop, dT);
        sleep_for(10ms);

    }

    system("pause");
    return 0;
}