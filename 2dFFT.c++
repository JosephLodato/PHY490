#include <iostream>
#include <complex>
#include <valarray>

using namespace std;

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

// ****************************************************************************
// ********************** JWST Miror Extraction Methods ***********************
// ****************************************************************************

// Helper function to calculate angle from center of the image
double calculateAngle(int x, int y, int centerX, int centerY)
{
    return atan2(y - centerY, x - centerX) * (180.0 / M_PI); // Convert to degrees
}

// Analyze diffraction spikes to extract mirror geometry
void analyzeDiffractionSpikes(const vector<CArray> &matrix, int width, int height, const char *channelName)
{
    int centerX = width / 2;
    int centerY = height / 2;

    // Number of angle bins (360 degrees / bin width)
    const int numBins = 360;
    std::vector<double> angleBins(numBins, 0.0); // Accumulate intensity per angle

    // Loop through the matrix and calculate angle for each point
    for (int x = 0; x < width; ++x)
    {
        for (int y = 0; y < height; ++y)
        {
            // Skip the center (DC component)
            if (x == centerX && y == centerY)
                continue;

            double angle = calculateAngle(x, y, centerX, centerY);
            // Normalize the angle to the range [0, 360)
            if (angle < 0)
                angle += 360;

            // Convert angle to bin index
            int bin = static_cast<int>(angle);
            if (bin >= 0 && bin < numBins)
            {
                angleBins[bin] += abs(matrix[x][y]); // Accumulate intensity
            }
        }
    }

    // Now analyze the angleBins to detect spikes
    std::vector<double> spikeAngles;

    // the 0.56 corresponds to 56% of max intensity, play with this value as this will not work for all images
    double peakThreshold = 0.56 * *std::max_element(angleBins.begin(), angleBins.end()); 
    std::vector<double> finalSpikeAngles;

    for (int i = 0; i < numBins; ++i)
    {
        if (angleBins[i] > peakThreshold)
        {
            spikeAngles.push_back(i); // Record the angle of the spike
        }
    }
    
    std::sort(spikeAngles.begin(), spikeAngles.end());

    // Use std::unique with a lambda to remove adjacent angles within 1 degree
    auto end = unique(spikeAngles.begin(), spikeAngles.end(), [](int a, int b) {
        return abs(a - b) <= 1; // Remove if the difference is 1 or less
    });
    spikeAngles.erase(end, spikeAngles.end());


    // Print out the detected spike angles and their corresponding intensities
    std::cerr << "Detected spike angles for " << channelName << " channel:\n";
    for (auto angle : spikeAngles)
    {
        std::cerr << "Angle: " << angle << " degrees, Intensity: " << angleBins[angle] << "\n"; // used for debuging 
        std::cout << angle << " " << angleBins[angle] << "\n"; // passed on to python script
    }
}

int main()
{
    vector<CArray> red, green, blue;
    int width = 0, height = 0;

    // This string is the path to the txt file representation of the image you want to read in
    // See README for more info if needed
    const char filename[] = "JWST.txt";

    // Call the readFile function
    readfile(red, green, blue, width, height, filename);

    // Perform 2D FFT on all color channels
    fft2D(red);
    fft2D(green);
    fft2D(blue);

    // Shift the DC componant of the transform back into the center of the image
    // Rather then on the 4 corner like it was before.
    fftShift(red);
    fftShift(green);
    fftShift(blue);


    // Analyze diffraction spikes for mirror geometry
    // We are only analyzing one of the color channels, I chose red but you can pick whichever channel works best
    analyzeDiffractionSpikes(red, width, height, "Red");
    // analyzeDiffractionSpikes(green, width, height, "Green");
    // analyzeDiffractionSpikes(blue, width, height, "Blue");
    return 0;
}
