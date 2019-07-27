/*
 * File:   ImgP.hh
 * Author: Mohamed TarekIbnZiad
 *
 * Created on September 16, 2015, 7:55 AM
 */
#pragma once
#include <vector>
#include <math/mpz_class.hh>
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/highgui/highgui.hpp"

#include<iostream>

using namespace std;
using namespace cv;

class ImgP {
 public:
    ImgP();
    Mat AverageFilterOpenCV(Mat src);
    Mat AverageFilter(Mat src, vector < vector<double> > &filter);
    Mat SobelFilterOpenCV(Mat src);
    Mat SobelFilter(Mat src);
    Mat convolution(Mat src, vector < vector<double> > &filter);
    Mat AdjustBrightness(Mat src, int value);
    Mat NegativeImage(Mat src);
    Mat AddSaltPepperNoise(Mat src);
    Mat equalHist(Mat src);
    Mat highBoostFilter(Mat src, double &A);
    Mat EdgeDetectionFilter(Mat src, char &type, char &direction);

};
