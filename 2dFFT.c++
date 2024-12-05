#include <iostream>
#include <complex>
#include <cmath>
#include <vector>
#include <valarray> // For FFT
#include <chrono>
#include <thread>
#include <map>

using namespace std;
using namespace std::this_thread;     // sleep_for, sleep_until
using namespace std::chrono_literals; // ns, us, ms, s, h, etc.
using std::chrono::system_clock;

#define M_PI 3.14159265358979323846

using Complex = std::complex<double>;
using CArray = std::valarray<Complex>;

int cutoffs = 0;
/*
 * HOW TO USE FILES:
 * 1. Run readImage.py with this command
 *      python3 readImage.py [image path] >> img.txt
 * 2. Run this file with this command
 *      ./2dFFT.o >> FFTIMG.txt
 * 3. Run createImg.py to turn it back into an image
 *      python3 createImg.py FFTImg.txt transformed.jpg 8

*/

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

    printf("%d %d\n", width, height);

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

// Inverse FFT function
void ifft(CArray &x)
{
    // Conjugate the complex numbers
    for (auto &elem : x)
        elem = std::conj(elem);

    // Apply FFT
    fft(x);

    // Conjugate the complex numbers again
    for (auto &elem : x)
        elem = std::conj(elem);

    // Scale the result
    // x /= x.size();
    for (auto &value : x)
    {
        value /= static_cast<double>(x.size());
    }
}

// Function to apply iFFT on all rows
void ifft2D_rows(vector<CArray> &matrix)
{
    for (auto &row : matrix)
    {
        ifft(row);
    }
}

// 2D Inverse FFT function
void ifft2D(vector<CArray> &matrix, int width, int height)
{
    // Conjugate before applying FFT for inverse calculation
    for (auto &row : matrix)
        for (auto &val : row)
            val = std::conj(val);

    // Apply FFT on all rows
    fft2D_rows(matrix);

    // Transpose the matrix
    transpose(matrix);

    // Apply FFT on all rows (which were originally columns)
    fft2D_rows(matrix);

    // Transpose back to the original layout
    transpose(matrix);

    // Conjugate and normalize the result by total size to complete iFFT
    const size_t totalSize = matrix.size() * matrix[0].size();
    for (auto &row : matrix)
        for (auto &val : row)
            val = std::conj(val) / static_cast<double>(width * height); // Conjugate again and scale down
}

// Function to shift the FFT output
void fftShift(vector<CArray> &matrix)
{
    const size_t rows = matrix.size();
    const size_t cols = matrix[0].size();
    const size_t halfRows = rows / 2;
    const size_t halfCols = cols / 2;

    // Swap quadrants diagonally
    for (size_t i = 0; i < halfRows; ++i)
    {
        for (size_t j = 0; j < halfCols; ++j)
        {
            // Top-left <--> Bottom-right
            std::swap(matrix[i][j], matrix[i + halfRows][j + halfCols]);

            // Top-right <--> Bottom-left
            std::swap(matrix[i][j + halfCols], matrix[i + halfRows][j]);
        }
    }
}

/*
 * Methods that opperate on trtransformed images:
 */

// Function to calculate the Euclidean distance
double distanceFromCenter(int u, int v, int centerU, int centerV)
{
    return sqrt((u - centerU) * (u - centerU) + (v - centerV) * (v - centerV));
}

// High-pass filter: keeps frequencies above the cutoff
void highPassFilter(vector<CArray> &matrix, int width, int height, double cutoff)
{
    int centerU = width / 2;
    int centerV = height / 2;

    for (int u = 0; u < width; ++u)
    {
        for (int v = 0; v < height; ++v)
        {
            double dist = distanceFromCenter(u, v, centerU, centerV);
            if (dist < cutoff)
            {
                matrix[u][v] = 0; // Zero out low frequencies
                cutoffs++;
                // std::cerr << "!Pixel Cut off\n";
            }
        }
    }
}

// Low-pass filter: keeps frequencies below the cutoff
void lowPassFilter(vector<CArray> &matrix, int width, int height, double cutoff)
{
    int centerU = width / 2;
    int centerV = height / 2;

    for (int u = 0; u < width; ++u)
    {
        for (int v = 0; v < height; ++v)
        {
            double dist = distanceFromCenter(u, v, centerU, centerV);
            if (dist > cutoff)
            {
                matrix[u][v] = 0; // Zero out high frequencies
            }
        }
    }
}

// Band-pass filter: keeps frequencies within a specific range
void bandPassFilter(vector<CArray> &matrix, int width, int height, double lowCutoff, double highCutoff)
{
    int centerU = width / 2;
    int centerV = height / 2;

    for (int u = 0; u < width; ++u)
    {
        for (int v = 0; v < height; ++v)
        {
            double dist = distanceFromCenter(u, v, centerU, centerV);
            if (dist < lowCutoff || dist > highCutoff)
            {
                matrix[u][v] = 0; // Zero out frequencies outside the band
            }
        }
    }
}

// ****************************************************************************
// ********************** JWST Miror Extraction Methods ***********************
// ****************************************************************************

// Helper function to calculate angle from center
double calculateAngle(int x, int y, int centerX, int centerY)
{
    return atan2(y - centerY, x - centerX) * (180.0 / M_PI); // Convert to degrees
}

// Analyze diffraction spikes to extract mirror geometry
void analyzeDiffractionSpikes(const vector<CArray> &matrix, int width, int height, const char *channelName)
{
    int centerX = width / 2;
    int centerY = height / 2;

    // // Find the maximum intensity (excluding DC component)
    // double maxIntensity = 0.0;
    // for (int x = 0; x < width; ++x)
    // {
    //     for (int y = 0; y < height; ++y)
    //     {
    //         if (x == centerX && y == centerY)
    //             continue; // Skip DC component
    //         maxIntensity = max(maxIntensity, abs(matrix[x][y]));
    //     }
    // }
    // double threshold = 0.9 * maxIntensity; // Set threshold as 50% of the max intensity

    // Threshold to identify high-intensity areas (spikes)
        double threshold = 0.18 * abs(matrix[centerX][centerY]); // 10% of DC component THis got 6/8 Spikes might be good enough

    // Map to store angles and intensities of spikes
    std::map<double, double> spikeAngles;

    for (int x = 0; x < width; ++x)
    {
        for (int y = 0; y < height; ++y)
        {
            if (abs(matrix[x][y]) > threshold)
            {
                double angle = calculateAngle(x, y, centerX, centerY);
                double distance = distanceFromCenter(x, y, centerX, centerY);

                // Normalize angle to 0-360 degrees
                if (angle < 0)
                    angle += 360;

                // Save the strongest intensity for each angle
                if (spikeAngles.find(angle) == spikeAngles.end() || abs(matrix[x][y]) > spikeAngles[angle])
                {
                    spikeAngles[angle] = abs(matrix[x][y]);
                }
            }
        }
    }

    // Print out spike angles for debugging
    std::cerr << "Spike analysis for " << channelName << " channel:\n";
    for (const auto &entry : spikeAngles)
    {
        std::cerr << "Angle: " << entry.first << " degrees, Intensity: " << entry.second << "\n";
    }
}

// Use this command to execute. It gets rid of using two commands and an middle txt file
//     ./a.o | python3 createImg.py transformed.jpg 0
int main()
{
    vector<CArray> red, green, blue;
    int width = 0, height = 0;

    const char filename[] = "JWST.txt";
    // const char filename[] = "img.txt";

    // Call the readFile function
    readfile(red, green, blue, width, height, filename);

    // Perform 2D FFT on all color channels
    fft2D(red);
    fft2D(green);
    fft2D(blue);

    // A fucntion to shift the DC componant of the transform back into the center of the image
    // Rather then on the 4 cornors like it was before.
    fftShift(red);
    fftShift(green);
    fftShift(blue);

    // Analyze diffraction spikes for mirror geometry
    analyzeDiffractionSpikes(red, width, height, "Red");
    analyzeDiffractionSpikes(green, width, height, "Green");
    analyzeDiffractionSpikes(blue, width, height, "Blue");

    // std::cerr << "Num pixels cut off" << cutoffs << endl; // Goes to stderr

    for (int y = 0; y < height; y++)
    {
        for (int x = 0; x < width; x++)
        {
            printf("%d %d %f %f %f\n", x, y, abs(red[x][y]), abs(green[x][y]), abs(blue[x][y]));
        }
    }
    return 0;
}
