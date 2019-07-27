/*
 * File:   test_Paillier_Image.cc
 * Author: Mohamed TarekIbnZiad
 *
 * Created on September 16, 2015, 7:48 AM
 */

#include <assert.h>
#include <vector>
#include <crypto/paillier.hh>
#include <Img.hh>
#include <ImgP.hh>
#include <crypto/gm.hh>
#include <NTL/ZZ.h>
#include <gmpxx.h>
#include <math/util_gmp_rand.h>

#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/highgui/highgui.hpp"

#include <ctime>

#include<iostream>

#define PI 3.14159265359

using namespace cv;
using namespace std;
using namespace NTL;


/*
 * convert an image from Mat type to double type used as input in encryption
 */
vector < vector<double> >
Mat2Double(Mat src)
{
    vector < vector<double> > A(src.rows, vector<double>(src.cols));
    Scalar s;
    for (int i = 0; i < src.rows; i++) {
        for (int j = 0; j < src.cols; j++){
            s = src.at<uchar>(i,j);
            A[i][j] = (int)s[0];
        }
    }
    return A;
}

/*
 * convert an image from double type to Mat type (opencv)
 */
Mat
Double2Mat(vector < vector<double> > &A)
{
    int rows = A.size();
    int cols = A.data()->size();
    Mat dst(rows,cols,CV_8UC1);
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++){
            dst.at<uchar>(i,j) = (int)A[i][j];  //Through away the float part or round it?
            //dst.at<uchar>(i,j) = round(A[i][j]);
        }
    }
    return dst;
}


/*
 * convert an image from encnum type used in encryption to Mat type (opencv)
 */
Mat
encnum2Mat(vector < vector<encnum> > &A)
{
    int rows = A.size();
    int cols = A.data()->size();
    Mat dst(rows,cols,CV_8UC1);
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++){
            dst.at<uchar>(i,j) = (int)A[i][j].mantissa.get_d(); //double to int
        }
    }
    return dst;
}


double
getRMSE(const Mat& I1, const Mat& I2)
{
    Mat s1;
    absdiff(I1, I2, s1);       // |I1 - I2|
    s1.convertTo(s1, CV_32F);  // cannot make a square on 8 bits
    s1 = s1.mul(s1);           // |I1 - I2|^2

    Scalar s = sum(s1);        // sum elements per channel
    double sse = s.val[0] + s.val[1] + s.val[2]; // sum channels
    double mse  = sse / (double)(I1.channels() * I1.total());
    double rmse = sqrt(mse);

    return rmse;
 }

/*
 * Get the absolute values for a resultant image
 */
vector < vector <double> >
getAbsolute(vector < vector<double> > &A)
{
    int rows = A.size();
    int cols = A.data()->size();
    vector < vector<double> > B(rows, vector<double>(cols));
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++){
            if (A[i][j] < 0 )
                B[i][j] = A[i][j] * -1;
            else
                B[i][j] = A[i][j];
        }
    }
    return B;
}

/*
 * Add the absolute values of 2 resultant images (H and V)
 * IF accepted and implemented on Mob --> combine it with get absolute
 */
vector < vector <double> >
AddAbsolute(vector < vector<double> > &A, vector < vector<double> > &B)
{
    int rows = A.size();
    int cols = A.data()->size();
    vector < vector<double> > C(rows, vector<double>(cols));
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++){
                C[i][j] = A[i][j] + B[i][j];
        }
    }
    return C;
}

/*
 * Creates an average filter with a specified kernel size
 */
vector < vector <double> >
CreateBoxFilter(int kSize){

    vector < vector <double> > filter (kSize, vector<double>(kSize));
    for (int i = 0; i < kSize; i++)
        for(int j = 0;j< kSize; j++)
            filter[i][j] = 1.0/(kSize*kSize);
    return filter;
}

/*
 * Creates a Gaussian filter with a specified kernel size
 */
vector < vector <double> >
CreateGaussianFilter(int kSize, double sigma){

    vector < vector <double> > filter (kSize, vector<double>(kSize));
    if (sigma < 0)
        sigma = 0.3*((kSize-1)*0.5 - 1) + 0.8;

    double alpha = 1/(2*PI*sigma*sigma);
    //(0,0) is in the center
    for (int x = -kSize/2; x <= kSize/2; x++){
        for(int y = -kSize/2;y<= kSize/2; y++){
            //Needs To be checked! http://prntscr.com/8h6q9b
            filter[x+(kSize/2)][y+(kSize/2)] = alpha*exp(-(x*x + y*y)/(2*sigma*sigma));
        }
    }
    return filter;
}

/*
 * Post processing for erosion
 */
vector < vector <double> >
erodePost(vector < vector<double> > &A)
{
    int rows = A.size();
    int cols = A.data()->size();
    vector < vector<double> > B(rows, vector<double>(cols));
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++){
            if (A[i][j] == (9*255) ) //All the values under the mask are 255s
                B[i][j] = 1*255;
            else
                B[i][j] = 0*255;
        }
    }
    return B;
}

/*
 * Post processing for dilation
 */
vector < vector <double> >
dilatePost(vector < vector<double> > &A)
{
    int rows = A.size();
    int cols = A.data()->size();
    vector < vector<double> > B(rows, vector<double>(cols));
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++){
            if (A[i][j] >= 1*255 && A[i][j] <= 9*255 )
                B[i][j] = 1*255;
            else //All the values under the mask are 0s
                B[i][j] = 0*255;
        }
    }
    return B;
}

/*
 * Takes an image as input
 * returns a 2d vector with the 1st columns represents colors and
 * the 2nd represents their freq.
 * The vector length is the total intensity levels G
 */
vector < vector <double> >
countHist(vector < vector<double> > src)
{
    int rows = src.size();
    int cols = src.data()->size();
    vector <int> Frequencies (256);
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++){
            Frequencies[src[i][j]]++;
        }
    }

    vector < vector<double> > ColorFreq(256, vector<double>(2));
    int G = 0;
    for (int i=0; i<256; i++){
        if(Frequencies[i] != 0){
            ColorFreq[G][0] = i;
            ColorFreq[G][1] = Frequencies[i];
            G++;
        }
    }
    return ColorFreq;
}

/*
 * Post processing Histogram function that constructs the new image based on the new histogram counts
 */
Mat
updateHist(vector < vector<double> > &newFreqHist, Mat src)
{
    int rows = src.rows;
    int cols = src.cols;
    int G = newFreqHist.size();

    vector <int> newColors (256);
    for (int i = 0; i<256; i++)
        newColors[i] = i;

    for (int i = 0; i<G; i++){
        newColors[newFreqHist[i][0]] = round(newFreqHist[i][1]);
    }

    Mat dst(rows,cols,CV_8UC1);
    Scalar s;
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++){
            s = src.at<uchar>(i,j);
            dst.at<uchar>(i,j) = newColors[(int)s[0]];
        }
    }
    return dst;
}


void
testImagePaillier()
{
    //Measuring time
    struct timespec t0,t1;
    uint64_t t;

    //Generate Paillier key
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t0);
    gmp_randstate_t randstate;
    gmp_randinit_default(randstate);
    gmp_randseed_ui(randstate,time(NULL));

    auto sk = Paillier_priv::keygen(randstate,256,2); //600 ,256 , 128
    Paillier_priv pp(sk,randstate);

    auto pk = pp.pubkey();
    mpz_class n = pk[0];
    Paillier p(pk,randstate);
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t1);
    t = (((uint64_t)t1.tv_sec) - ((uint64_t)t0.tv_sec) )* 1000000000 + (t1.tv_nsec - t0.tv_nsec);
    cerr << "Key Generation : "<<  ((double)t/1000000) <<"ms" << endl;
    cout << "key ready" << endl;

    //read image
    Mat src = imread( "cameraman.JPG", 0); //cameraman.JPG
    imshow( "Original Image", src );

    int rows = src.rows;
    int cols = src.cols;
    vector < vector<double> > A(rows, vector<double>(cols));
    A = Mat2Double(src);

    //Encrypt image
    vector < vector<encnum> > A_enc;
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t0);
    A_enc = p.encryptMatrix(A);
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t1);
    t = (((uint64_t)t1.tv_sec) - ((uint64_t)t0.tv_sec) )* 1000000000 + (t1.tv_nsec - t0.tv_nsec);
    cerr << "public encryption: "<<  ((double)t/1000000) <<"ms per image" << endl;
    cout << "Image encrypted " << endl;
    //imshow( "Encrypted Image", encnum2Mat(A_enc) );

    //The client only can decrypt a matrix
    vector < vector<double> > A_dec;
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t0);
    A_dec = pp.decryptMatrix(A_enc);
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t1);
    t = (((uint64_t)t1.tv_sec) - ((uint64_t)t0.tv_sec) )* 1000000000 + (t1.tv_nsec - t0.tv_nsec);
    cerr << "private decryption: "<<  ((double)t/1000000) <<"ms per image" << endl;
    cout << "Image decrypted " << endl;
    Mat dst;
    dst = Double2Mat(A_dec);
    imshow( "Decrypted Image", dst );
    cout << "Root Mean Square Error in decryption: "<<getRMSE(src,dst)<<endl;

    //Create image to work in plain domain
    ImgP img1P;

    //Create image to work in encrypted domain
    Img img1(A_enc);
    img1.setPbKey(p);

    //Debug code
    cout << "Test sub Paillier numbers...\n" << flush;
    double a = 255 , b = 5;
    encnum a_enc = p.encrypt_f(a);
    encnum b_enc = p.encrypt_f(b);
    b_enc = p.constMult_f(1.0 , b_enc);
    //b_enc = p.encrypt_f(pp.decrypt_f(b_enc)); //Re-encrypt
    encnum sum = p.sub_f(a_enc, b_enc);
    cout << a << " - " << b << " = " << pp.decrypt_f(sum)<< endl<< endl;
    cout << "---------------------------------------\n" << flush;

    //Negate image in plain domain
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t0);
    dst = img1P.NegativeImage(src);
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t1);
    t = (((uint64_t)t1.tv_sec) - ((uint64_t)t0.tv_sec) )* 1000000000 + (t1.tv_nsec - t0.tv_nsec);
    cerr << "Negative Image (Plain): "<<  ((double)t/1000000) <<"ms" << endl;
    imshow( "Negative Image", dst );

    //Negate image in encrypted domain
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t0);
    vector < vector<encnum> > Neg_enc = img1.NegativeImageH(A_enc);
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t1);
    t = (((uint64_t)t1.tv_sec) - ((uint64_t)t0.tv_sec) )* 1000000000 + (t1.tv_nsec - t0.tv_nsec);
    cerr << "Negative Image (Enc): "<<  ((double)t/1000000) <<"ms per image" << endl;
    //imshow( "Encrypted Negative Image", encnum2Mat(Neg_enc) );

    vector < vector<double> > Neg_dec = pp.decryptMatrix(Neg_enc);
    imshow( "Decrypted Negative Image", Double2Mat(Neg_dec) );
    cout << "Root Mean Square Error in Negative image: "<<getRMSE(dst,Double2Mat(Neg_dec))<<endl;

    //Test Brightness
    int value = 10;
    //Increase brightness in plain domain
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t0);
    dst = img1P.AdjustBrightness(src, value);
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t1);
    t = (((uint64_t)t1.tv_sec) - ((uint64_t)t0.tv_sec) )* 1000000000 + (t1.tv_nsec - t0.tv_nsec);
    cerr << "Brightness Image (Plain): "<<  ((double)t/1000000) <<"ms per image" << endl;
    imshow( "Brightness Image", dst );

    //Increase brightness in encrypted domain
    encnum value_enc = p.encrypt_f(value);
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t0);
    vector < vector<encnum> > Bright_enc = img1.AdjustBrightnessH(A_enc, value_enc);
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t1);
    t = (((uint64_t)t1.tv_sec) - ((uint64_t)t0.tv_sec) )* 1000000000 + (t1.tv_nsec - t0.tv_nsec);
    cerr << "Brightness Image (Enc): "<<  ((double)t/1000000) <<"ms per image" << endl;
    //imshow( "Encrypted Brightness Image", encnum2Mat(Bright_enc) );

    vector < vector<double> > Bright_dec = pp.decryptMatrix(Bright_enc);
    imshow( "Decrypted Brightness Image", Double2Mat(Bright_dec) );
    cout << "Root Mean Square Error in Brightness image: "<<getRMSE(dst,Double2Mat(Bright_dec))<<endl;
    //post-processing to fix the image
    //Bright_dec = FixRange(Bright_dec);
    //imshow( "Decrypted Brightness Image", Double2Mat(Bright_dec) );

    //Test convolution
    vector < vector<double> > filter;
    //filter = {{1,0,-1},{2,0,-2},{1,0,-1}};//Vertical edges
    filter = {{1,2,1},{0,0,0},{-1,-2,-1}};  //Horizontal edges
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t0);
    vector < vector<encnum> > conv_enc = img1.convolutionH(A_enc, filter);
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t1);
    t = (((uint64_t)t1.tv_sec) - ((uint64_t)t0.tv_sec) )* 1000000000 + (t1.tv_nsec - t0.tv_nsec);
    cerr << "Convolution (Enc): "<<  ((double)t/1000000) <<"ms per image" << endl;
    //imshow( "Encrypted convoluted Image", encnum2Mat(conv_enc) );
    vector < vector<double> > conv_dec = pp.decryptMatrix(conv_enc);
    //Temp solution to the -ve value problem in edge detection
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t0);
    vector < vector<double> > conv_dec_fixed = getAbsolute(conv_dec);
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t1);
    t = (((uint64_t)t1.tv_sec) - ((uint64_t)t0.tv_sec) )* 1000000000 + (t1.tv_nsec - t0.tv_nsec);
    cerr << "Convolution (Post Processing): "<<  ((double)t/1000000) <<"ms per image" << endl;
    imshow( "Decrypted convoluted Image", Double2Mat(conv_dec_fixed) );

    //Test convolution in Plain domain
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t0);
    dst = img1P.convolution(src, filter);
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t1);
    t = (((uint64_t)t1.tv_sec) - ((uint64_t)t0.tv_sec) )* 1000000000 + (t1.tv_nsec - t0.tv_nsec);
    cerr << "Convolution (Plain): "<<  ((double)t/1000000) <<"ms per image" << endl;
    //Take care of rounding and borders to reduce error
    imshow( "convoluted Image", dst );
    cout << "Root Mean Square Error in convoluted image: "<<getRMSE(dst,Double2Mat(conv_dec))<<endl;
    cout << "Root Mean Square Error in convoluted image (Absolute): "<<getRMSE(dst,Double2Mat(conv_dec_fixed))<<endl;

    //Test Average filter
    //filter = {{1.0/9,1.0/9,1.0/9},{1.0/9,1.0/9,1.0/9},{1.0/9,1.0/9,1.0/9}};
    filter = CreateBoxFilter(3);
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t0);
    vector < vector<encnum> > Avg_enc = img1.convolutionH(A_enc, filter);
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t1);
    t = (((uint64_t)t1.tv_sec) - ((uint64_t)t0.tv_sec) )* 1000000000 + (t1.tv_nsec - t0.tv_nsec);
    cerr << "Average (Enc): "<<  ((double)t/1000000) <<"ms per image" << endl;
    //imshow( "Encrypted averaged Image", encnum2Mat(Avg_enc) );
    vector < vector<double> > Avg_dec = pp.decryptMatrix(Avg_enc);
    //vector < vector<double> > Avg_dec_d = Mat
    imshow( "Decrypted averaged Image", Double2Mat(Avg_dec) );

    //Test average filter in plain domain
    //dst = img1P.AverageFilter(src,filter);
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t0);
    dst = img1P.convolution(src, filter);
    vector < vector<double> >  aa = Mat2Double(dst);//Debug Line
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t1);
    t = (((uint64_t)t1.tv_sec) - ((uint64_t)t0.tv_sec) )* 1000000000 + (t1.tv_nsec - t0.tv_nsec);
    cerr << "Average (Plain): "<<  ((double)t/1000000) <<"ms per image" << endl;
    imshow( "Averaged Image", dst );
    cout << "Root Mean Square Error in Averaged image: "<<getRMSE(dst,Double2Mat(Avg_dec))<<endl; //Error due to rounding

    //Test Gaussian Filtering
    //filter = CreateGaussianFilter(5,1);

    //Laplacian:Calculates the Laplacian of an image.
    //filter = {{0,1,0},{1,-4,1},{0,1,0}};

    //Test Histogram equalization
    //Pre-processing
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t0);
    vector < vector<double> > ColorFreq = countHist(Mat2Double(src));
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t1);
    t = (((uint64_t)t1.tv_sec) - ((uint64_t)t0.tv_sec) )* 1000000000 + (t1.tv_nsec - t0.tv_nsec);
    cerr << "Histogram Equalization (Enc_Pre): "<<  ((double)t/1000000) <<"ms" << endl;

    vector < vector<encnum> > ColorFreq_enc = p.encryptMatrix(ColorFreq);

    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t0);
    vector < vector<encnum> > newColorFreq_enc = img1.equalHistH(ColorFreq_enc, src.rows, src.cols);
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t1);
    t = (((uint64_t)t1.tv_sec) - ((uint64_t)t0.tv_sec) )* 1000000000 + (t1.tv_nsec - t0.tv_nsec);
    cerr << "Histogram Equalization (Enc): "<<  ((double)t/1000000) <<"ms" << endl;

    vector < vector<double> > newColorFreq_dec = pp.decryptMatrix(newColorFreq_enc);

    //Post-processing
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t0);
    Mat equalizedSrc = updateHist(newColorFreq_dec, src);
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t1);
    t = (((uint64_t)t1.tv_sec) - ((uint64_t)t0.tv_sec) )* 1000000000 + (t1.tv_nsec - t0.tv_nsec);
    cerr << "Histogram Equalization (Enc_post): "<<  ((double)t/1000000) <<"ms" << endl;
    imshow( "Decrypted Equalized Image", equalizedSrc );
    //OpenCV --> TODO
    //imshow( "Equalized Image Open CV", srcB );
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t0);
    dst = img1P.equalHist(src);
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t1);
    t = (((uint64_t)t1.tv_sec) - ((uint64_t)t0.tv_sec) )* 1000000000 + (t1.tv_nsec - t0.tv_nsec);
    cerr << "Histogram Equalization (Plain): "<<  ((double)t/1000000) <<"ms" << endl;
    imshow( "Equalized Image ", dst );
    cout << "Root Mean Square Error in Equalized image: "<<getRMSE(dst,equalizedSrc)<<endl;
    //To eliminate error here --> use precision = 0.00000001 with 128-bit key

    //Test Edge filters
    //In plain domain
    char direction = 'H'; //'H' for horizontal and 'V' for vertical
    char type = 'S'; //'R' Roberts, 'P' Prewitt, 'S' Sobel, 'R' Robinson, 'K' Kirsch
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t0);
    dst = img1P.EdgeDetectionFilter(src, type, direction);
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t1);
    t = (((uint64_t)t1.tv_sec) - ((uint64_t)t0.tv_sec) )* 1000000000 + (t1.tv_nsec - t0.tv_nsec);
    cerr << "Edge filter in H (Plain): "<<  ((double)t/1000000) <<"ms" << endl;
    imshow( "Edge Image H", dst );
    //In encrypted domain
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t0);
    vector < vector<encnum> > Edge_enc = img1.EdgeDetectionFilterH(A_enc, type, direction);
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t1);
    t = (((uint64_t)t1.tv_sec) - ((uint64_t)t0.tv_sec) )* 1000000000 + (t1.tv_nsec - t0.tv_nsec);
    cerr << "Edge filter (Enc): "<<  ((double)t/1000000) <<"ms" << endl;
    imshow( "Encrypted Edge Image", encnum2Mat(Edge_enc) );
    vector < vector<double> > Edge_dec = pp.decryptMatrix(Edge_enc);
    imshow( "Decrypted Edge Image H", Double2Mat(Edge_dec) );
    cout << "Root Mean Square Error in Edge image H: "<<getRMSE(dst,Double2Mat(Edge_dec))<<endl;
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t0);
    vector < vector<double> > Edge_dec_fixed = getAbsolute(Edge_dec);
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t1);
    t = (((uint64_t)t1.tv_sec) - ((uint64_t)t0.tv_sec) )* 1000000000 + (t1.tv_nsec - t0.tv_nsec);
    cerr << "Edge filter in H (Enc_post): "<<  ((double)t/1000000) <<"ms" << endl;
    imshow( "Decrypted Edge Image in H (Absolute)", Double2Mat(Edge_dec_fixed) );
    cout << "Root Mean Square Error in Edge image H (Absolute): "<<getRMSE(dst,Double2Mat(Edge_dec_fixed))<<endl;

    //In vertical direction
    direction = 'V'; //'H' for horizontal and 'V' for vertical
    type = 'S'; //'R' Roberts, 'P' Prewitt, 'S' Sobel, 'R' Robinson, 'K' Kirsch
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t0);
    dst = img1P.EdgeDetectionFilter(src, type, direction);
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t1);
    t = (((uint64_t)t1.tv_sec) - ((uint64_t)t0.tv_sec) )* 1000000000 + (t1.tv_nsec - t0.tv_nsec);
    cerr << "Edge filter in V (Plain): "<<  ((double)t/1000000) <<"ms" << endl;
    imshow( "Edge Image V", dst );
    //In encrypted domain
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t0);
    Edge_enc = img1.EdgeDetectionFilterH(A_enc, type, direction);
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t1);
    t = (((uint64_t)t1.tv_sec) - ((uint64_t)t0.tv_sec) )* 1000000000 + (t1.tv_nsec - t0.tv_nsec);
    cerr << "Edge filter (Enc): "<<  ((double)t/1000000) <<"ms" << endl;
    //imshow( "Encrypted Edge Image V", encnum2Mat(Edge_enc) );
    Edge_dec = pp.decryptMatrix(Edge_enc);
    imshow( "Decrypted Edge Image V", Double2Mat(Edge_dec) );
    cout << "Root Mean Square Error in Edge image V: "<<getRMSE(dst,Double2Mat(Edge_dec))<<endl;
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t0);
    vector < vector<double> > Edge_dec_fixedV = getAbsolute(Edge_dec);
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t1);
    t = (((uint64_t)t1.tv_sec) - ((uint64_t)t0.tv_sec) )* 1000000000 + (t1.tv_nsec - t0.tv_nsec);
    cerr << "Edge filter in V (Enc_post): "<<  ((double)t/1000000) <<"ms" << endl;
    imshow( "Decrypted Edge Image in V (Absolute)", Double2Mat(Edge_dec_fixedV) );
    cout << "Root Mean Square Error in Edge image V (Absolute): "<<getRMSE(dst,Double2Mat(Edge_dec_fixedV))<<endl;

    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t0);
    vector < vector<double> > EdgeTotal = AddAbsolute(Edge_dec_fixed,Edge_dec_fixedV);
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t1);
    t = (((uint64_t)t1.tv_sec) - ((uint64_t)t0.tv_sec) )* 1000000000 + (t1.tv_nsec - t0.tv_nsec);
    cerr << "Edge filter (Add abs_post): "<<  ((double)t/1000000) <<"ms" << endl;
    imshow( "Decrypted Edge Image Total ", Double2Mat(EdgeTotal) );


    //Test morphological operations
    // Load an image from file
    Mat srcB = imread( "butterfly-11.jpg", 0); //0 --> gray, 1 --> RGB image
    srcB = srcB > 128; // Convert gray to binary
    //show the loaded image
    imshow( "Original Binary Image", srcB );

    Mat dstB;
    dstB = srcB.clone(); //Copy of src as initial value

    int element_shape = MORPH_RECT;
    Mat element = getStructuringElement(element_shape, Size(3, 3), Point(-1, -1) );
    vector < vector<double> > element_d = Mat2Double(element);

    int rowsB = srcB.rows;
    int colsB = srcB.cols;
    vector < vector<double> > B(rowsB, vector<double>(colsB));
    B = Mat2Double(srcB);
    //Encrypt image
    vector < vector<encnum> > B_enc = p.encryptMatrix(B);
    //imshow( "Encrypted Binary Image", encnum2Mat(B_enc) );

    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t0);
    vector < vector<encnum> > morph_enc = img1.morphH(B_enc, element_d);
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t1);
    t = (((uint64_t)t1.tv_sec) - ((uint64_t)t0.tv_sec) )* 1000000000 + (t1.tv_nsec - t0.tv_nsec);
    cerr << "Morph (Enc): "<<  ((double)t/1000000) <<"ms" << endl;
    //imshow( "Encrypted Morphological Image", encnum2Mat(morph_enc) );
    vector < vector<double> > morph_dec = pp.decryptMatrix(morph_enc);
    //Post processing for erosion/dilation
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t0);
    vector < vector<double> > morph_dec_erode = erodePost(morph_dec);
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t1);
    t = (((uint64_t)t1.tv_sec) - ((uint64_t)t0.tv_sec) )* 1000000000 + (t1.tv_nsec - t0.tv_nsec);
    cerr << "Morph (Enc_post_Erode): "<<  ((double)t/1000000) <<"ms" << endl;
    imshow( "Decrypted eroded Image", Double2Mat(morph_dec_erode) );

    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t0);
    erode(srcB, dstB, element);
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t1);
    t = (((uint64_t)t1.tv_sec) - ((uint64_t)t0.tv_sec) )* 1000000000 + (t1.tv_nsec - t0.tv_nsec);
    cerr << "Morph (Plain_Erode): "<<  ((double)t/1000000) <<"ms" << endl;
    imshow("Erode OpenCV",dstB);
    cout << "Root Mean Square Error in eroded image: "<<getRMSE(dstB,Double2Mat(morph_dec_erode))<<endl;

    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t0);
    vector < vector<double> > morph_dec_dilate = dilatePost(morph_dec);
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t1);
    t = (((uint64_t)t1.tv_sec) - ((uint64_t)t0.tv_sec) )* 1000000000 + (t1.tv_nsec - t0.tv_nsec);
    cerr << "Morph (Enc_post_Dilate): "<<  ((double)t/1000000) <<"ms" << endl;
    imshow( "Decrypted dilated Image", Double2Mat(morph_dec_dilate) );
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t0);
    dilate(srcB, dstB, element);
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t1);
    t = (((uint64_t)t1.tv_sec) - ((uint64_t)t0.tv_sec) )* 1000000000 + (t1.tv_nsec - t0.tv_nsec);
    cerr << "Morph (Plain_Dilate): "<<  ((double)t/1000000) <<"ms" << endl;
    imshow("Dilate OpenCV",dstB);
    cout << "Root Mean Square Error in dilated image: "<<getRMSE(dstB,Double2Mat(morph_dec_dilate))<<endl;

    //wait for a key press infinitely
    waitKey(0);
}

void
testImagePaillier_N(int n_iteration)
{
    //Measuring time
    struct timespec t0,t1;
    uint64_t t;

    //Generate Paillier key
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t0);
    gmp_randstate_t randstate;
    gmp_randinit_default(randstate);
    gmp_randseed_ui(randstate,time(NULL));

    auto sk = Paillier_priv::keygen(randstate,128,2); //600 ,256 , 128
    Paillier_priv pp(sk,randstate);

    auto pk = pp.pubkey();
    mpz_class n = pk[0];
    Paillier p(pk,randstate);
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t1);
    t = (((uint64_t)t1.tv_sec) - ((uint64_t)t0.tv_sec) )* 1000000000 + (t1.tv_nsec - t0.tv_nsec);
    cerr << "Key Generation : "<<  ((double)t/1000000) <<"ms" << endl;
    cout << "key ready" << endl;

    //read image
    Mat src = imread( "cameraman.JPG", 0);
    imshow( "Original Image", src );

    int rows = src.rows;
    int cols = src.cols;
    vector < vector<double> > A(rows, vector<double>(cols));
    A = Mat2Double(src);

    //Encrypt image
    vector < vector<encnum> > A_enc;
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t0);
    for (size_t i = 0; i < n_iteration; i++) {
        A_enc = p.encryptMatrix(A);
    }
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t1);
    t = (((uint64_t)t1.tv_sec) - ((uint64_t)t0.tv_sec) )* 1000000000 + (t1.tv_nsec - t0.tv_nsec);
    cerr << "public encryption: "<<  ((double)t/1000000)/n_iteration <<"ms per image" << endl;
    cout << "Image encrypted " << endl;
    //imshow( "Encrypted Image", encnum2Mat(A_enc) );

    //The client only can decrypt a matrix
    vector < vector<double> > A_dec;
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t0);
    for (size_t i = 0; i < n_iteration; i++) {
        A_dec = pp.decryptMatrix(A_enc);
    }
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t1);
    t = (((uint64_t)t1.tv_sec) - ((uint64_t)t0.tv_sec) )* 1000000000 + (t1.tv_nsec - t0.tv_nsec);
    cerr << "private decryption: "<<  ((double)t/1000000)/n_iteration <<"ms per image" << endl;
    cout << "Image decrypted " << endl;
    Mat dst;
    dst = Double2Mat(A_dec);
    imshow( "Decrypted Image", dst );
    cout << "Root Mean Square Error in decryption: "<<getRMSE(src,dst)<<endl;

    //Create image to work in plain domain
    ImgP img1P;

    //Create image to work in encrypted domain
    Img img1(A_enc);
    img1.setPbKey(p);

    //Negate image in plain domain
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t0);
    for (size_t i = 0; i < n_iteration; i++) {
        dst = img1P.NegativeImage(src);
    }
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t1);
    t = (((uint64_t)t1.tv_sec) - ((uint64_t)t0.tv_sec) )* 1000000000 + (t1.tv_nsec - t0.tv_nsec);
    cerr << "Negative Image (Plain): "<<  ((double)t/1000000)/n_iteration <<"ms" << endl;
    imshow( "Negative Image", dst );

    //Negate image in encrypted domain
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t0);
    vector < vector<encnum> > Neg_enc;
    for (size_t i = 0; i < n_iteration; i++) {
        Neg_enc = img1.NegativeImageH(A_enc);
    }
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t1);
    t = (((uint64_t)t1.tv_sec) - ((uint64_t)t0.tv_sec) )* 1000000000 + (t1.tv_nsec - t0.tv_nsec);
    cerr << "Negative Image (Enc): "<<  ((double)t/1000000)/n_iteration <<"ms per image" << endl;
    //imshow( "Encrypted Negative Image", encnum2Mat(Neg_enc) );

    vector < vector<double> > Neg_dec = pp.decryptMatrix(Neg_enc);
    imshow( "Decrypted Negative Image", Double2Mat(Neg_dec) );
    cout << "Root Mean Square Error in Negative image: "<<getRMSE(dst,Double2Mat(Neg_dec))<<endl;

    //Test Brightness
    int value = 10;
    //Increase brightness in plain domain
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t0);
    for (size_t i = 0; i < n_iteration; i++) {
        dst = img1P.AdjustBrightness(src, value);
    }
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t1);
    t = (((uint64_t)t1.tv_sec) - ((uint64_t)t0.tv_sec) )* 1000000000 + (t1.tv_nsec - t0.tv_nsec);
    cerr << "Brightness Image (Plain): "<<  ((double)t/1000000)/n_iteration <<"ms per image" << endl;
    imshow( "Brightness Image", dst );

    //Increase brightness in encrypted domain
    encnum value_enc = p.encrypt_f(value);
    vector < vector<encnum> > Bright_enc;
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t0);
    for (size_t i = 0; i < n_iteration; i++) {
        Bright_enc = img1.AdjustBrightnessH(A_enc, value_enc);
    }
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t1);
    t = (((uint64_t)t1.tv_sec) - ((uint64_t)t0.tv_sec) )* 1000000000 + (t1.tv_nsec - t0.tv_nsec);
    cerr << "Brightness Image (Enc): "<<  ((double)t/1000000)/n_iteration <<"ms per image" << endl;
    //imshow( "Encrypted Brightness Image", encnum2Mat(Bright_enc) );

    vector < vector<double> > Bright_dec = pp.decryptMatrix(Bright_enc);
    imshow( "Decrypted Brightness Image", Double2Mat(Bright_dec) );
    cout << "Root Mean Square Error in Brightness image: "<<getRMSE(dst,Double2Mat(Bright_dec))<<endl;
    //post-processing to fix the image
    //Bright_dec = FixRange(Bright_dec);
    //imshow( "Decrypted Brightness Image", Double2Mat(Bright_dec) );

    //Test convolution
    vector < vector<double> > filter;
    //filter = {{1,0,-1},{2,0,-2},{1,0,-1}};//Vertical edges
    filter = {{1,2,1},{0,0,0},{-1,-2,-1}};  //Horizontal edges
    vector < vector<encnum> > conv_enc;
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t0);
    for (size_t i = 0; i < n_iteration; i++) {
        conv_enc = img1.convolutionH(A_enc, filter);
    }
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t1);
    t = (((uint64_t)t1.tv_sec) - ((uint64_t)t0.tv_sec) )* 1000000000 + (t1.tv_nsec - t0.tv_nsec);
    cerr << "Convolution (Enc): "<<  ((double)t/1000000)/n_iteration <<"ms per image" << endl;
    //imshow( "Encrypted convoluted Image", encnum2Mat(conv_enc) );
    vector < vector<double> > conv_dec = pp.decryptMatrix(conv_enc);
    //Temp solution to the -ve value problem in edge detection
    vector < vector<double> > conv_dec_fixed;
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t0);
    for (size_t i = 0; i < n_iteration; i++) {
        conv_dec_fixed = getAbsolute(conv_dec);
    }
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t1);
    t = (((uint64_t)t1.tv_sec) - ((uint64_t)t0.tv_sec) )* 1000000000 + (t1.tv_nsec - t0.tv_nsec);
    cerr << "Convolution (Post Processing): "<<  ((double)t/1000000)/n_iteration <<"ms per image" << endl;
    imshow( "Decrypted convoluted Image", Double2Mat(conv_dec_fixed) );

    //Test convolution in Plain domain
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t0);
    for (size_t i = 0; i < n_iteration; i++) {
        dst = img1P.convolution(src, filter);
    }
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t1);
    t = (((uint64_t)t1.tv_sec) - ((uint64_t)t0.tv_sec) )* 1000000000 + (t1.tv_nsec - t0.tv_nsec);
    cerr << "Convolution (Plain): "<<  ((double)t/1000000)/n_iteration <<"ms per image" << endl;
    //Take care of rounding and borders to reduce error
    imshow( "convoluted Image", dst );
    cout << "Root Mean Square Error in convoluted image: "<<getRMSE(dst,Double2Mat(conv_dec))<<endl;
    cout << "Root Mean Square Error in convoluted image (Absolute): "<<getRMSE(dst,Double2Mat(conv_dec_fixed))<<endl;

    //Test Average filter
    //filter = {{1.0/9,1.0/9,1.0/9},{1.0/9,1.0/9,1.0/9},{1.0/9,1.0/9,1.0/9}};
    filter = CreateBoxFilter(3);
    vector < vector<encnum> > Avg_enc;
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t0);
    for (size_t i = 0; i < n_iteration; i++) {
        Avg_enc = img1.convolutionH(A_enc, filter);
    }
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t1);
    t = (((uint64_t)t1.tv_sec) - ((uint64_t)t0.tv_sec) )* 1000000000 + (t1.tv_nsec - t0.tv_nsec);
    cerr << "Average (Enc): "<<  ((double)t/1000000)/n_iteration <<"ms per image" << endl;
    //imshow( "Encrypted averaged Image", encnum2Mat(Avg_enc) );
    vector < vector<double> > Avg_dec = pp.decryptMatrix(Avg_enc);
    //vector < vector<double> > Avg_dec_d = Mat
    imshow( "Decrypted averaged Image", Double2Mat(Avg_dec) );

    //Test average filter in plain domain
    //dst = img1P.AverageFilter(src,filter);
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t0);
    for (size_t i = 0; i < n_iteration; i++) {
        dst = img1P.convolution(src, filter);
    }
    //vector < vector<double> >  aa = Mat2Double(dst);//Debug Line
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t1);
    t = (((uint64_t)t1.tv_sec) - ((uint64_t)t0.tv_sec) )* 1000000000 + (t1.tv_nsec - t0.tv_nsec);
    cerr << "Average (Plain): "<<  ((double)t/1000000)/n_iteration <<"ms per image" << endl;
    imshow( "Averaged Image", dst );
    cout << "Root Mean Square Error in Averaged image: "<<getRMSE(dst,Double2Mat(Avg_dec))<<endl; //Error due to rounding

    //Test Gaussian Filtering
    //filter = CreateGaussianFilter(5,1);

    //Laplacian:Calculates the Laplacian of an image.
    //filter = {{0,1,0},{1,-4,1},{0,1,0}};

    //Test Histogram equalization
    //Pre-processing
    vector < vector<double> > ColorFreq;
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t0);
    for (size_t i = 0; i < n_iteration; i++) {
        ColorFreq = countHist(Mat2Double(src));
    }
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t1);
    t = (((uint64_t)t1.tv_sec) - ((uint64_t)t0.tv_sec) )* 1000000000 + (t1.tv_nsec - t0.tv_nsec);
    cerr << "Histogram Equalization (Enc_Pre): "<<  ((double)t/1000000)/n_iteration <<"ms" << endl;

    vector < vector<encnum> > ColorFreq_enc = p.encryptMatrix(ColorFreq);

    vector < vector<encnum> > newColorFreq_enc;
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t0);
    for (size_t i = 0; i < n_iteration; i++) {
        newColorFreq_enc = img1.equalHistH(ColorFreq_enc, src.rows, src.cols);
    }
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t1);
    t = (((uint64_t)t1.tv_sec) - ((uint64_t)t0.tv_sec) )* 1000000000 + (t1.tv_nsec - t0.tv_nsec);
    cerr << "Histogram Equalization (Enc): "<<  ((double)t/1000000)/n_iteration <<"ms" << endl;

    vector < vector<double> > newColorFreq_dec = pp.decryptMatrix(newColorFreq_enc);

    //Post-processing
    Mat equalizedSrc;
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t0);
    for (size_t i = 0; i < n_iteration; i++) {
        equalizedSrc = updateHist(newColorFreq_dec, src);
    }
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t1);
    t = (((uint64_t)t1.tv_sec) - ((uint64_t)t0.tv_sec) )* 1000000000 + (t1.tv_nsec - t0.tv_nsec);
    cerr << "Histogram Equalization (Enc_post): "<<  ((double)t/1000000)/n_iteration <<"ms" << endl;
    imshow( "Decrypted Equalized Image", equalizedSrc );
    //OpenCV --> TODO
    //imshow( "Equalized Image Open CV", srcB );
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t0);
    for (size_t i = 0; i < n_iteration; i++) {
        dst = img1P.equalHist(src);
    }
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t1);
    t = (((uint64_t)t1.tv_sec) - ((uint64_t)t0.tv_sec) )* 1000000000 + (t1.tv_nsec - t0.tv_nsec);
    cerr << "Histogram Equalization (Plain): "<<  ((double)t/1000000)/n_iteration <<"ms" << endl;
    imshow( "Equalized Image ", dst );
    cout << "Root Mean Square Error in Equalized image: "<<getRMSE(dst,equalizedSrc)<<endl;
    //To eliminate error here --> use precision = 0.00000001 with 128-bit key

    //Test Edge filters
    //In plain domain
    char direction = 'H'; //'H' for horizontal and 'V' for vertical
    char type = 'S'; //'R' Roberts, 'P' Prewitt, 'S' Sobel, 'R' Robinson, 'K' Kirsch
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t0);
    for (size_t i = 0; i < n_iteration; i++) {
        dst = img1P.EdgeDetectionFilter(src, type, direction);
    }
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t1);
    t = (((uint64_t)t1.tv_sec) - ((uint64_t)t0.tv_sec) )* 1000000000 + (t1.tv_nsec - t0.tv_nsec);
    cerr << "Edge filter in H (Plain): "<<  ((double)t/1000000)/n_iteration <<"ms" << endl;
    imshow( "Edge Image H", dst );
    //In encrypted domain
    vector < vector<encnum> > Edge_enc;
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t0);
    for (size_t i = 0; i < n_iteration; i++) {
        Edge_enc = img1.EdgeDetectionFilterH(A_enc, type, direction);
    }
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t1);
    t = (((uint64_t)t1.tv_sec) - ((uint64_t)t0.tv_sec) )* 1000000000 + (t1.tv_nsec - t0.tv_nsec);
    cerr << "Edge filter H(Enc): "<<  ((double)t/1000000)/n_iteration <<"ms" << endl;
    imshow( "Encrypted Edge Image H", encnum2Mat(Edge_enc) );
    vector < vector<double> > Edge_dec = pp.decryptMatrix(Edge_enc);
    imshow( "Decrypted Edge Image H", Double2Mat(Edge_dec) );
    cout << "Root Mean Square Error in Edge image H: "<<getRMSE(dst,Double2Mat(Edge_dec))<<endl;
    vector < vector<double> > Edge_dec_fixed;
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t0);
    for (size_t i = 0; i < n_iteration; i++) {
        Edge_dec_fixed = getAbsolute(Edge_dec);
    }
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t1);
    t = (((uint64_t)t1.tv_sec) - ((uint64_t)t0.tv_sec) )* 1000000000 + (t1.tv_nsec - t0.tv_nsec);
    cerr << "Edge filter in H(Enc_post): "<<  ((double)t/1000000)/n_iteration <<"ms" << endl;
    imshow( "Decrypted Edge Image H(Absolute)", Double2Mat(Edge_dec_fixed) );
    cout << "Root Mean Square Error in Edge image H(Absolute): "<<getRMSE(dst,Double2Mat(Edge_dec_fixed))<<endl;

    //In vertical direction
    direction = 'V'; //'H' for horizontal and 'V' for vertical
    type = 'S'; //'R' Roberts, 'P' Prewitt, 'S' Sobel, 'R' Robinson, 'K' Kirsch
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t0);
    for (size_t i = 0; i < n_iteration; i++) {
        dst = img1P.EdgeDetectionFilter(src, type, direction);
    }
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t1);
    t = (((uint64_t)t1.tv_sec) - ((uint64_t)t0.tv_sec) )* 1000000000 + (t1.tv_nsec - t0.tv_nsec);
    cerr << "Edge filter in V (Plain): "<<  ((double)t/1000000)/n_iteration <<"ms" << endl;
    imshow( "Edge Image V", dst );
    //In encrypted domain
    //vector < vector<encnum> > Edge_enc;
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t0);
    for (size_t i = 0; i < n_iteration; i++) {
        Edge_enc = img1.EdgeDetectionFilterH(A_enc, type, direction);
    }
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t1);
    t = (((uint64_t)t1.tv_sec) - ((uint64_t)t0.tv_sec) )* 1000000000 + (t1.tv_nsec - t0.tv_nsec);
    cerr << "Edge filter V(Enc): "<<  ((double)t/1000000)/n_iteration <<"ms" << endl;
    imshow( "Encrypted Edge Image V", encnum2Mat(Edge_enc) );
    Edge_dec = pp.decryptMatrix(Edge_enc);
    imshow( "Decrypted Edge Image V", Double2Mat(Edge_dec) );
    cout << "Root Mean Square Error in Edge image V: "<<getRMSE(dst,Double2Mat(Edge_dec))<<endl;
    vector < vector<double> > Edge_dec_fixedV;
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t0);
    for (size_t i = 0; i < n_iteration; i++) {
        Edge_dec_fixedV = getAbsolute(Edge_dec);
    }
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t1);
    t = (((uint64_t)t1.tv_sec) - ((uint64_t)t0.tv_sec) )* 1000000000 + (t1.tv_nsec - t0.tv_nsec);
    cerr << "Edge filter in V(Enc_post): "<<  ((double)t/1000000)/n_iteration <<"ms" << endl;
    imshow( "Decrypted Edge Image V(Absolute)", Double2Mat(Edge_dec_fixedV) );
    cout << "Root Mean Square Error in Edge image V(Absolute): "<<getRMSE(dst,Double2Mat(Edge_dec_fixedV))<<endl;

	vector < vector<double> > EdgeTotal;
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t0);
	  for (size_t i = 0; i < n_iteration; i++) {
        EdgeTotal = AddAbsolute(Edge_dec_fixed,Edge_dec_fixedV);
    }
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t1);
    t = (((uint64_t)t1.tv_sec) - ((uint64_t)t0.tv_sec) )* 1000000000 + (t1.tv_nsec - t0.tv_nsec);
    cerr << "Edge filter (Add abs_post): "<<  ((double)t/1000000)/n_iteration <<"ms" << endl;
    imshow( "Decrypted Edge Image Total ", Double2Mat(EdgeTotal) );

    //Test morphological operations
    // Load an image from file
    Mat srcB = imread( "butterfly-11.jpg", 0); //0 --> gray, 1 --> RGB image
    srcB = srcB > 128; // Convert gray to binary
    //show the loaded image
    imshow( "Original Binary Image", srcB );

    Mat dstB;
    dstB = srcB.clone(); //Copy of src as initial value

    int element_shape = MORPH_RECT;
    Mat element = getStructuringElement(element_shape, Size(3, 3), Point(-1, -1) );
    vector < vector<double> > element_d = Mat2Double(element);

    int rowsB = srcB.rows;
    int colsB = srcB.cols;
    vector < vector<double> > B(rowsB, vector<double>(colsB));
    B = Mat2Double(srcB);
    //Encrypt image
    vector < vector<encnum> > B_enc = p.encryptMatrix(B);
    //imshow( "Encrypted Binary Image", encnum2Mat(B_enc) );

    vector < vector<encnum> > morph_enc;
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t0);
    for (size_t i = 0; i < n_iteration; i++) {
        morph_enc = img1.morphH(B_enc, element_d);
    }
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t1);
    t = (((uint64_t)t1.tv_sec) - ((uint64_t)t0.tv_sec) )* 1000000000 + (t1.tv_nsec - t0.tv_nsec);
    cerr << "Morph (Enc): "<<  ((double)t/1000000)/n_iteration <<"ms" << endl;
    //imshow( "Encrypted Morphological Image", encnum2Mat(morph_enc) );
    vector < vector<double> > morph_dec = pp.decryptMatrix(morph_enc);
    //Post processing for erosion/dilation
    vector < vector<double> > morph_dec_erode;
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t0);
    for (size_t i = 0; i < n_iteration; i++) {
        morph_dec_erode = erodePost(morph_dec);
    }
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t1);
    t = (((uint64_t)t1.tv_sec) - ((uint64_t)t0.tv_sec) )* 1000000000 + (t1.tv_nsec - t0.tv_nsec);
    cerr << "Morph (Enc_post_Erode): "<<  ((double)t/1000000)/n_iteration <<"ms" << endl;
    imshow( "Decrypted eroded Image", Double2Mat(morph_dec_erode) );

    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t0);
    for (size_t i = 0; i < n_iteration; i++) {
        erode(srcB, dstB, element);
    }
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t1);
    t = (((uint64_t)t1.tv_sec) - ((uint64_t)t0.tv_sec) )* 1000000000 + (t1.tv_nsec - t0.tv_nsec);
    cerr << "Morph (Plain_Erode): "<<  ((double)t/1000000)/n_iteration <<"ms" << endl;
    imshow("Erode OpenCV",dstB);
    cout << "Root Mean Square Error in eroded image: "<<getRMSE(dstB,Double2Mat(morph_dec_erode))<<endl;

    vector < vector<double> > morph_dec_dilate;
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t0);
    for (size_t i = 0; i < n_iteration; i++) {
        morph_dec_dilate = dilatePost(morph_dec);
    }
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t1);
    t = (((uint64_t)t1.tv_sec) - ((uint64_t)t0.tv_sec) )* 1000000000 + (t1.tv_nsec - t0.tv_nsec);
    cerr << "Morph (Enc_post_Dilate): "<<  ((double)t/1000000)/n_iteration <<"ms" << endl;
    imshow( "Decrypted dilated Image", Double2Mat(morph_dec_dilate) );
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t0);
    for (size_t i = 0; i < n_iteration; i++) {
        dilate(srcB, dstB, element);
    }
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t1);
    t = (((uint64_t)t1.tv_sec) - ((uint64_t)t0.tv_sec) )* 1000000000 + (t1.tv_nsec - t0.tv_nsec);
    cerr << "Morph (Plain_Dilate): "<<  ((double)t/1000000)/n_iteration <<"ms" << endl;
    imshow("Dilate OpenCV",dstB);
    cout << "Root Mean Square Error in dilated image: "<<getRMSE(dstB,Double2Mat(morph_dec_dilate))<<endl;

    //wait for a key press infinitely
    waitKey(0);
}

/*
 * Testing New Floating point Paillier with Images
 */
int main(int argc, char** argv) {

    testImagePaillier(); /*For generating images*/
    //testImagePaillier_N(10); /*N = 10 to measure time*/
    return 0;
}
