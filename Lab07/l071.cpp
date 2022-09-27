#include <stdio.h>
#include <iostream>
#include <opencv2/opencv.hpp>
#include <opencv2/highgui/highgui.hpp>

using namespace std;
using namespace cv;

ofstream output("results.txt");

int main(int argc, char** argv )
{
    if ( argc != 2 )
    {
        printf("usage: DisplayImage.out <Image_Path>\n");
        return -1;
    }

    Mat image, gray, blurred, hough;
    image = imread(argv[1], 1);
    hough = imread(argv[1], 1);

    if ( !image.data )
    {
        printf("No image data \n");
        return -1;
    }
    
    vector<Vec3f> vect;
    gray = imread( argv[1], IMREAD_GRAYSCALE );
    medianBlur(gray, blurred, 3);
    HoughCircles(blurred, vect, HOUGH_GRADIENT, 1, 160, 100, 70, 40, 180);
    
    int p = 0, n = 0, d = 0, q = 0, h = 0, sum;
    for(size_t i = 0; i < vect.size(); i++)
    {
        Vec3i v = vect[i];
        Point c = Point(v[0], v[1]);
        int radius = v[2];
        if(radius > 70 && radius < 81) {
            circle(hough, c, radius, Scalar(0, 255, 0), 3, LINE_AA);
            d += 1;
        } else if(radius < 90) {
            circle(hough, c, radius, Scalar(255, 255, 0), 3, LINE_AA);
            p += 1;
        } else if(radius < 105) {
            circle(hough, c, radius, Scalar(0, 255, 255), 3, LINE_AA);
            n += 1;
        } else if(radius < 120) {
            circle(hough, c, radius, Scalar(255, 0, 0), 3, LINE_AA);
            q += 1;
        } else if(radius > 160 && radius < 172) {
            circle(hough, c, radius, Scalar(0, 0, 255), 3, LINE_AA);
            h += 1;
        }
    }

    output << p << " Pennies" << endl;
    output << n << " Nickels" << endl;
    output << d << " Dimes" << endl;
    output << q << " Quarters" << endl;
    output << h << " Half Dollars" << endl;
    sum = p + (5*n) + (10*d) + (25*q) + (50*h);
    output << "Total sum: $" << sum/100 << "." << sum%100;

    imwrite("./imagec.jpg", hough);
    waitKey(0);

    return 0;
}