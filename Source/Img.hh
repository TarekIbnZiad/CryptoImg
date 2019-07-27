/*
 * File:   Img.hh
 * Author: Mohamed TarekIbnZiad
 *
 * Created on September 16, 2015, 4:11 AM
 */

#ifndef IMG_HH
#define	IMG_HH

#endif	/* IMG_HH */

#pragma once
#include <vector>
#include <crypto/paillier.hh>

#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/highgui/highgui.hpp"

#include<iostream>

using namespace std;
using namespace cv;

class Img {
 public:
    Img(vector < vector<encnum> > &_A);
    void setPbKey( Paillier &_p);

    vector < vector<encnum> > NegativeImageH(const vector < vector<encnum> > &src_enc);
    vector < vector<encnum> > AdjustBrightnessH(const vector < vector<encnum> > &src_enc, const encnum &value_enc);
    vector < vector<encnum> > convolutionH(const vector < vector<encnum> > &src_enc, const vector < vector<double> > &filter);
    vector < vector<encnum> > morphH(const vector < vector<encnum> > &src_enc, const vector < vector<double> > &element);
    vector < vector<encnum> > equalHistH(const vector < vector<encnum> > &src_enc, int &rows, int &cols);
    vector < vector<encnum> > highBoostFilterH(const vector < vector<encnum> > &src_enc, double &A);
    vector < vector<encnum> > highBoostFilterH2(const vector < vector<encnum> > &src_enc, const vector < vector<encnum> > &Avg_enc, double &A);
    vector < vector<encnum> > EdgeDetectionFilterH(const vector < vector<encnum> > &src_enc, char &type, char &direction);

    /* Paillier key */
    Paillier *p;

    /* Source image */
    vector < vector<encnum> > A;

};
