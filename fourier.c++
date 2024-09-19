/* To Do:
 * Fix timestep/tstop to not cut off mid frequency
 * Refactor/clean up code
 *  

*/


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


array<double, 2> calculationHelper(double f, double t, double M)
{
    // M represents the total time the program for each frequency runs through
    double inputSignal = (cos(2 * M_PI * t));   
    double windowFunction = 0.5 * (1 - cos(2 * M_PI * t / M)); 

    // double windowFunction = 0.5*(1+cos( (2*M_PI*(t + 2*M)) / M )); // 5 is the delta time so 2M
    
    inputSignal *= windowFunction;


    std::complex<double> result = inputSignal * (cexp(-2 * M_PI * I * f * t));

    

    array<double, 2> temp;
    temp[0] = result.real();
    temp[1] = result.imag();

    return temp;
}

array<double, 2> calculateIntegral(vector<array<double, 2>> points, double dT)
{
    // Calculates the integral using a vector of all the points and the trapezoidal method

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

    return temp;
}

void DFT(double f, double t, double tStop, double dT, bool usingAnim, double M)
{
    vector<array<double, 2>> points;

    // Calcualtes all the that make up the integral for the given wrapping frequency
    while (t < tStop)
    {
        t += dT;
        points.push_back(calculationHelper(f, t, M)); // Keeps a list of the results for all t values
    }

    // Printing results in either anim format or plot format
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

        array<double, 2> temp = calculateIntegral(points, dT);
        for (int i = 0; i < 50; i++)
        {
            printf("c %f %f 0.1\n", temp[0], temp[1]);
        }

        printf("F\n");
        points.clear();
    }
    else if (!usingAnim)
    {
        array<double, 2> temp = calculateIntegral(points, dT);
        printf("%f %f\n", f, hypot(temp[0], temp[1]));
    }
}

int main(int argc, char const *argv[])
{
    bool usingAnim = false;
    if (argc != 2)
    {
        printf("INVALID USE - Please Use one of the formats: \n ./a.out a | anim or ./a.out p | plot\n");
        return 0;
    }
    else if (*argv[1] == 'a')
    {
        usingAnim = true;
    }

    double f = 0;       // Starting Wrapping Frequency
    double t = -2.5;    // Stating Time
    double dT = 0.001;  // Time Step

    double dF = 0.001; // Frequency Step
    // double dF = 2; // Frequency Step


    double tStop = 10; // Max Time simulated for each f
    double fStop = 5;   // Max frequency the simulation runs to

    // calculate points across time
    while (f <= fStop)
    {
        DFT(f, t, tStop, dT, usingAnim, abs(t) + abs(tStop));
        f += dF;
    }

    return 0;
}