#include <iostream>

using namespace std;



/*
 * Called on the command line with 
 *    ./imageReadTest.c++ 
 * no need to pass in a file path
*/

void readfile(double **&red, double **&blue, double **&green, int &width, int &height, char *filename)
{
    FILE *fp = fopen(filename, "r");

    if (!fp)
    {
        printf("You told me to read file %s but I can't so I am sad and going to go away now.\n", filename);
        exit(-2);
    }

    fscanf(fp, "%d %d", &width, &height);

    printf("Width = %d, height = %d\n", width, height);

    red = (double **)malloc(sizeof(double *) * width);
    green = (double **)malloc(sizeof(double *) * width);
    blue = (double **)malloc(sizeof(double *) * width);

    for (int i = 0; i < width; i++)
    {
        red[i] = (double *)malloc(sizeof(double) * height);
        green[i] = (double *)malloc(sizeof(double) * height);
        blue[i] = (double *)malloc(sizeof(double) * height);
    }

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
}

int main(int argc, char const *argv[])
{
    double **red = nullptr;
    double **green = nullptr;
    double **blue = nullptr;

    int width = 0, height = 0;

    char filename[] = "img.txt";

    // Call the readfile function
    readfile(red, blue, green, width, height, filename);

    std::cout << "Image Data (Width x Height): " << width << " x " << height << std::endl;
    for (int i = 0; i < width; i++)
    {
        for (int j = 0; j < height; j++)
        {
            printf("Pixel (%d, %d): Red = %f, Green = %f, Blue = %f\n", i, j, red[i][j], green[i][j], blue[i][j]);
        }
    }
}
