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

// Function to apply FFT on all rows
void fft2D_rows(vector<CArray> &matrix)
{
    for (auto &row : matrix)
    {
        fft(row);
    }
}

// Transpose the matrix
void transpose(vector<CArray> &matrix)
{
    const size_t rows = matrix.size();
    const size_t cols = matrix[0].size();
    vector<CArray> transposed(cols, CArray(rows));

    for (size_t i = 0; i < rows; ++i)
    {
        for (size_t j = 0; j < cols; ++j)
        {
            transposed[j][i] = matrix[i][j];
        }
    }
    matrix = std::move(transposed);
}

// Function to perform 2D FFT
void fft2D(vector<CArray> &matrix)
{
    // Perform FFT on all rows
    fft2D_rows(matrix);

    // Transpose to turn rows into columns
    transpose(matrix);

    // Perform FFT on all "rows" (which were originally columns)
    fft2D_rows(matrix);

    // Transpose back to the original layout
    transpose(matrix);
}

// Helper function to apply a window function to each row
void applyWindowFunction(CArray &signal, double M)
{
    for (size_t i = 0; i < signal.size(); i++)
    {
        double t = i * (M / signal.size()); // Scale t across the window
        double windowFunction = 0.5 * (1 - cos(2 * M_PI * t / M));
        signal[i] *= windowFunction;
    }
}

int main()
{
    const int N = 8; // Small 8x8 checkerboard for demonstration

    // Create a 2D checkerboard pattern
    vector<CArray> checkerboard(N, CArray(N));

    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            // Alternate between 0 and 1 to create a checkerboard pattern
            checkerboard[i][j] = ((i / 2) % 2 == (j / 2) % 2) ? -1.0 : 1.0;
        }
    }

    // More complex signal
    // for (int i = 0; i < N; i++)
    // {
    //     for (int j = 0; j < N; j++)
    //     {
    //         checkerboard[i][j] = sin(2 * M_PI * i / N) + cos(2 * M_PI * j / N);
    //     }
    // }

    // Print the original checkerboard pattern
    cout << "Original Checkerboard Pattern:" << endl;
    for (const auto &row : checkerboard)
    {
        for (const auto &value : row)
        {
            cout << value.real() << " ";
        }
        cout << endl;
    }

    // Perform 2D FFT
    fft2D(checkerboard);

    // Output the 2D FFT result (real and imaginary parts)
    cout << "\n2D FFT Result (Real and Imaginary Parts):" << endl;
    for (const auto &row : checkerboard)
    {
        for (const auto &value : row)
        {
            cout << "(" << abs(value) << ") ";
        }
        cout << endl;
    }
}
