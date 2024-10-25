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

// Reads in images that have been passed through readImage.py and updates the red, green, and blue arrays, as well as height and width.
void readfile(vector<CArray> &red, vector<CArray> &green, vector<CArray> &blue, int &width, int &height, const char *filename)
{
    FILE *fp = fopen(filename, "r");

    if (!fp)
    {
        printf("You told me to read file %s but I can't so I am sad and going to go away now.\n", filename);
        exit(-2);
    }

    fscanf(fp, "%d %d", &width, &height);

    printf("Width = %d, height = %d\n", width, height);

    red.resize(width, CArray(height));
    green.resize(width, CArray(height));
    blue.resize(width, CArray(height));

    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < width; j++)
        {
            int ii, jj;
            double rhere, ghere, bhere;
            fscanf(fp, "%d %d : %lf %lf %lf", &ii, &jj, &rhere, &ghere, &bhere);

            if (i != ii || j != jj)
            {
                printf("Warning -- expected pixel %d, %d, but got pixel %d, %d\n", i, j, ii, jj);
                exit(-1);
            }

            red[j][i] = rhere;
            green[j][i] = ghere;
            blue[j][i] = bhere;
        }
    }
    fclose(fp);
}

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

int main()
{
    vector<CArray> red, green, blue;
    int width = 0, height = 0;

    const char filename[] = "img.txt";

    // Call the readFile function
    readfile(red, green, blue, width, height, filename);

    // Test to see if the image is imported correctly
    // for (int i = 0; i < width; i++)
    // {
    //     for (int j = 0; j < height; j++)
    //     {
    //         printf("Pixel (%d, %d): Red = %f, Green = %f, Blue = %f\n", i, j, red[i][j].real(), green[i][j].real(), blue[i][j].real());
    //     }
    // }


    // Perform 2D FFT on all color channels
    fft2D(red);
    fft2D(green);
    fft2D(blue);

    // print the output in a format the converter likes
    for (int i = 0; i < width; i++)
    {
        for (int j = 0; j < height; j++)
        {
            printf("%d %d %f %f %f\n", i, j, red[i][j].real(), green[i][j].real(), blue[i][j].real());
        }
    }

    return 0;
}
