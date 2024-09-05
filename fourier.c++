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

std::complex<double> complexExponential(double f, double t)
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

array<double, 2> calculateCenter(vector<array<double, 2>> points, double dT)
{
    // this must do the integral section of the equation
    // for visualizing it we must multiply my 1 / (t2-t1)
    // this gets dropped later

    int n = points.size();
    complex<double> integral(0.0, 0.0);

    // First and last points
    complex<double> first(points[0][0], points[0][1]);
    complex<double> last(points[n - 1][0], points[n - 1][1]);

    // Apply the trapezoidal rule
    integral += first + last;

    for (int i = 1; i < n - 1; i++)
    {
        complex<double> current(points[i][0], points[i][1]);
        integral += 2.0 * current;
    }

    integral *= (dT / 2.0);

    // Calculate the x/y coordinate in polar form
    double r = abs(integral);
    double theta = arg(integral);

    double x_polar = r * cos(theta);
    double y_polar = r * sin(theta);

    array<double, 2> temp;
    temp[0] = x_polar;
    temp[1] = y_polar;

    //printf("!Center Point Location: %f %f \n", temp[0], temp[1]);

    return temp;
}

void DFT(double f, double t, double tStop, double dT, bool usingAnim)
{
    vector<array<double, 2>> points;

    while (t < tStop)
    {
        t += dT;
        points.push_back(finalCalculation(f, t));
    }

    if (usingAnim)
    {

        // display all points
        for (int i = 0; i < points.size(); i++)
        {
            printf("c %f %f 0.01\n", points[i][0], points[i][1]);
            // printf("!c %f %f 0.01\n", points[i][0], points[i][1]);
        }
        printf("l -2 0 2 0\n");
        printf("l 0 -2 0 2\n");
        printf("t -2 2 \n Wrapping Frequency: %f \n", f);

        array<double, 2> temp = calculateCenter(points, dT);
        for (int i = 0; i < 50; i++)
        {
            printf("c %f %f 0.1\n", temp[0], temp[1]);
        }

        printf("F\n");
        points.clear();
    }
    else if(!usingAnim){
        array<double, 2> temp = calculateCenter(points, dT);
        printf("%f %f\n", f, temp[1]);    
    }
}

int main(int argc, char const *argv[])
{
    bool usingAnim;
    if(argc != 2){
        printf("INVALID USE: [a/g]");
        return 0;
    }
    else{
        if(*argv[1] == 'a')
        {
            usingAnim = true; 
        }
        else
        {
            usingAnim = false;
        }
    }
    
    
    vector<array<double, 2>> points;

    double f = 0;    // Frequency
    double t = 1.0;    // Time
    double dT = 0.001; // Time Step
    double dF = 0.001; // Frequency Step

    double tStop = 5; // Max Time simulated for each f
    double fStop = 5;
    // calculate points across time
    while (f <= 10)
    {
        f += dF;
        DFT(f, t, tStop, dT, usingAnim);
        if(usingAnim){sleep_for(10ms);}
    }

    return 0;
}