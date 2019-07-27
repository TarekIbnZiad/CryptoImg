/*
 * Main File For testing
 * Added by Mohamed TarekIbnZiad
 */

#include <assert.h>
#include <vector>
#include <crypto/paillier.hh>
#include <Img.hh>
#include <crypto/gm.hh>
#include <NTL/ZZ.h>
#include <gmpxx.h>
#include <math/util_gmp_rand.h>

#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/highgui/highgui.hpp"

#include <ctime>

#include<iostream>

using namespace cv;
using namespace std;
using namespace NTL;

/*
 * Apply average filter in Plain domain using open cv ready made function (blur)
 */
Mat
AverageFilterOpenCV(Mat src)
{
    //enter the kernel size
    int kerSize;
    cout << "Please enter the kernel size: \n";
    cin >> kerSize;

    Mat dst;
    //smooth the image in the "src" and save it to "dst"
    blur( src, dst, Size( kerSize, kerSize ) );

    return dst;
}

/*
 * Apply average filter in Plain domain
 */
Mat
AverageFilter(Mat src, vector < vector<mpz_class> > &filter)
{
    Mat dst = src.clone(); //Copy of src as initial value
    Scalar s;
    //Filter size
    int ro = filter.size();
    int co = filter.data()->size();

    //get summation of the filter elements to divide by
    mpz_class nFilter = 0;
    for (int i=0; i < ro; i++){
	     for (int j=0; j < co; j++){
            nFilter += filter[i][j];
       }
    }

    for (int i=1; i < src.rows - 2; i++){
	     for (int j=1; j < src.cols - 2; j++){
            mpz_class sum=0;
            for(int m=-1; m<=1; m++){
		            for(int n=-1; n<=1; n++){
                    s = src.at<uchar>(i+m,j+n);
                    sum += s.val[0] * filter[m+1][n+1];
                }
            }
            s.val[0] = sum.get_d()/nFilter.get_d();
            dst.at<uchar>(i,j) = s.val[0];
        }
    }
    return dst;
}

/*
 * Apply Sobel operator in Plain domain using open cv ready made function (Sobel)
 */
Mat
SobelFilterOpenCV(Mat src)
{
    int scale = 1;
    int delta = 0;
    int ddepth = CV_16S;

    Mat dst;

    // Generate grad_x and grad_y
    Mat grad_x, grad_y;
    Mat abs_grad_x, abs_grad_y;

    // Gradient X
    Sobel( src, grad_x, ddepth, 1, 0, 3, scale, delta, BORDER_DEFAULT );
    convertScaleAbs( grad_x, abs_grad_x );

    // Gradient Y
    Sobel( src, grad_y, ddepth, 0, 1, 3, scale, delta, BORDER_DEFAULT );
    convertScaleAbs( grad_y, abs_grad_y );

    // Total Gradient (approximate)
    addWeighted( abs_grad_x, 0.5, abs_grad_y, 0.5, 0, dst );

    return dst;
}

/*
 * Apply Sobel operator in Plain domain using indirect access
 */
Mat
SobelFilter(Mat src)
{
    int dx[3][3] = {{1,0,-1},{2,0,-2},{1,0,-1}};
    int dy[3][3] = {{1,2,1},{0,0,0},{-1,-2,-1}};

    Mat dst = src.clone(); //Copy of src as initial value
    //Vec3b s; // for colored images
    Scalar s;

    for (int i=1; i < src.rows - 2; i++){
	     for (int j=1; j < src.cols - 2; j++){
            // apply kernel in X and Y directions
            int sum_x=0;
            int sum_y=0;
            for(int m=-1; m<=1; m++){
		            for(int n=-1; n<=1; n++){
                    //s=cvGet2D(img,i+m,j+n); // get the (i,j) pixel value
                    s = src.at<uchar>(i+m,j+n);
                    sum_x += s.val[0] * dx[m+1][n+1];
                    sum_y += s.val[0] * dy[m+1][n+1];
                }
            }

            int sum=abs(sum_x)+abs(sum_y);
            s.val[0] = (sum>255)? 255:sum;
            dst.at<uchar>(i,j) = s.val[0];
            //dst.at<Vec2d>(i,j)[0] = s;
            //cvSet2D(dst,i,j,s); // set the (i,j) pixel value
        }
    }
    return dst;
}

/*
 * Apply convolution in Plain domain
 */
Mat
convolution(Mat src, vector < vector<mpz_class> > &filter)
{
    Mat dst = src.clone(); //Copy of src as initial value
    Scalar s;

    for (int i=1; i < src.rows - 1; i++){
	     for (int j=1; j < src.cols - 1; j++){
            mpz_class sum=0;
            for(int m=-1; m<=1; m++){
		            for(int n=-1; n<=1; n++){
                    if (filter[m+1][n+1] == 0) { //Not encrypted value
                        continue;
                    }
                    s = src.at<uchar>(i+m,j+n);
                    sum += s.val[0] * filter[m+1][n+1];
                }
            }
            s.val[0] = abs(sum.get_d());
            dst.at<uchar>(i,j) = s.val[0];
        }
    }
    return dst;
}

/*
 * used to Add/Sub constant value to the image brightness in plain domain
 */
Mat
AdjustBrightness(Mat src, int value)
{
    Mat dst = src.clone(); //Copy of src as initial value
    Scalar s;
    int newValue;
    for (int i=0; i < src.rows; i++){
	     for (int j=0; j < src.cols; j++){
            s = src.at<uchar>(i,j);
            newValue = s.val[0] + value;
            newValue = (newValue>255)? 255:
                       (newValue<0)?     0:newValue;
            dst.at<uchar>(i,j) = newValue;
        }
    }
    return dst;
}

/*
 * used to fix the image colors range after return from the homo domain
 */
vector < vector<mpz_class> >
FixRange(vector < vector<mpz_class> > &in)
{

    int rows = in.size();
    int cols = in.data()->size();
    int newValue;

    vector < vector<mpz_class> > out(rows, vector<mpz_class>(cols));
    for (int i=0; i < rows; i++){
	     for (int j=0; j < cols; j++){
            newValue = (int)in[i][j].get_d();
            newValue = (newValue>255)? 255:
                       (newValue<0)?     0:newValue;
            out[i][j] = (int)newValue;
        }
    }
    return out;
}

/*
 * used to get the negative image in plain domain
 */
Mat
NegativeImage(Mat src)
{
    Mat dst = src.clone(); //Copy of src as initial value
    Scalar s;
    int newValue;
    for (int i=0; i < src.rows; i++){
	     for (int j=0; j < src.cols; j++){
            s = src.at<uchar>(i,j);
            newValue = 255 - s.val[0];
            dst.at<uchar>(i,j) = newValue;
        }
    }
    return dst;
}

/*
 * used to add salt and pepper noise to the image
 */
Mat
AddSaltPepperNoise(Mat src)
{
    Mat saltpepper_noise = Mat::zeros(src.rows, src.cols,CV_8U);
    randu(saltpepper_noise,0,255);

    Mat black = saltpepper_noise < 10;
    Mat white = saltpepper_noise > 245;

    Mat saltpepper_img = src.clone();
    saltpepper_img.setTo(255,white);
    saltpepper_img.setTo(0,black);

    return saltpepper_img;
}

/*
 * convert an image from Mat type to mpz_class type used in encryption
 */
vector < vector<mpz_class> >
Mat2mpz(Mat src)
{
    vector < vector<mpz_class> > A(src.rows, vector<mpz_class>(src.cols));
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
 * convert an image from mpz_class type used in encryption to Mat type (opencv)
 */
Mat
mpz2Mat(vector < vector<mpz_class> > &A)
{
    int rows = A.size();
    int cols = A.data()->size();
    Mat dst(rows,cols,CV_8UC1);
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++){
            dst.at<uchar>(i,j) = (int)A[i][j].get_d(); //double to int
        }
    }
    return dst;
}

/*
 * used to divide by n in the average filter in time domain
 * would be removed later after implementing paillier encoding
 */
vector < vector<mpz_class> >
divideP (const vector < vector<mpz_class> > &src, const vector < vector<mpz_class> > &filter)
{
    //Src image size
    int rows = src.size();
    int cols = src.data()->size();
    //Filter size
    int ro = filter.size();
    int co = filter.data()->size();

    vector < vector<mpz_class> > dst (rows, vector<mpz_class>(cols));

    //get summation of the filter elements to divide by
    mpz_class sum = 0;
    for (int i=0; i < ro; i++){
	     for (int j=0; j < co; j++){
            sum += filter[i][j];
        }
    }

    for (int i=0; i < rows; i++){
	     for (int j=0; j < cols; j++){
            dst[i][j] = src[i][j]/sum;
        }
    }
    return dst;
}

static void
test_image()
{
    //create 2 empty windows
    //namedWindow( "Original Image" , CV_WINDOW_AUTOSIZE );
    //namedWindow( "Smoothed Image" , CV_WINDOW_AUTOSIZE );

    // Load an image from file
    Mat src = imread( "cameraman.JPG", 1 ); //0 --> gray, 1 --> RGB image

    //show the loaded image
    imshow( "Original Image", src );

    Mat dst;

    // Test Average filter by Opencv
    dst = AverageFilterOpenCV(src);
    imshow( "Opencv Smoothed Image", dst );

    //Add Salt and pepper noise:
    Mat saltpepper_img = AddSaltPepperNoise(src);
    imshow( "Noisy Image", saltpepper_img );

    dst = AverageFilterOpenCV(saltpepper_img);
    imshow( "Opencv Smoothed Image", dst );

    // Test Sobel filter by Me
    Mat src_gray;
    cvtColor(src, src_gray, CV_BGR2GRAY); //change the color image to grayscale image
    imshow( "Gray Image", src_gray );

    dst = SobelFilter(src_gray);
    imshow( "My Sobel Image", dst );

    // Test Sobel filter by Opencv
    dst = SobelFilterOpenCV(src_gray);
    imshow( "Opencv Sobel Image", dst );
}

static void
test_paillier()
{
    cout << "Test Paillier ...\n" << flush;

    gmp_randstate_t randstate;
    gmp_randinit_default(randstate);
    gmp_randseed_ui(randstate,time(NULL));

    auto sk = Paillier_priv::keygen(randstate,16,2); //600 ,256
    Paillier_priv pp(sk,randstate);

    auto pk = pp.pubkey();
    mpz_class n = pk[0];
    Paillier p(pk,randstate);

    //mpz_class pt0, pt1,m;
    //mpz_urandomm(pt0.get_mpz_t(),randstate,n.get_mpz_t());
    //mpz_urandomm(pt1.get_mpz_t(),randstate,n.get_mpz_t());
    //mpz_urandomm(m.get_mpz_t(),randstate,n.get_mpz_t());
    mpz_class pt0 = 2, pt1 = 3, m = 5; //instead of the random values
    cout << "pt0: "<< pt0 << endl;
    cout << "pt1: "<< pt1 << endl;

    mpz_class ct0 = p.encrypt(pt0);         cout << "ct0: "<< ct0 << endl;
    mpz_class ct1 = p.encrypt(pt1);         cout << "ct1: "<< ct1 << endl;
    mpz_class sum = p.add(ct0, ct1);        cout << "sum_enc: "<< sum << endl;
    mpz_class prod = p.constMult(m,ct0);    cout << "prod_enc: "<< prod << endl;
    mpz_class diff = p.sub(ct0, ct1);       cout << "diff_enc: "<< diff << endl;

    mpz_class ct0_dec = pp.decrypt(ct0);
    cout << "ct0_dec "<< ct0_dec.get_d() << endl;
    cout << "n " <<  n.get_d() <<endl;
    assert(pp.decrypt(ct0) == pt0);
    assert(pp.decrypt(ct1) == pt1);
    assert(pp.decrypt(sum) == (pt0+pt1)%n);
    mpz_class d = pt0 - pt1;
    if (d < 0) {
        d += n;
    }
    assert(pp.decrypt(diff) == d);
    assert(pp.decrypt(prod) == (m*pt0)%n);

    cout << "Test Paillier passed" << endl;

    cout << "Test Matrix Paillier ...\n" << flush;
    vector < vector<mpz_class> > A = {{1,2,3},{1,1,1},{4,5,6}};
    int Arows = A.size();
    int Acols = A.data()->size();
    vector < vector<mpz_class> > A_enc(Arows, vector<mpz_class>(Acols));
    A_enc = p.encryptMatrix(A);
    vector < vector<mpz_class> > A_dec(Arows, vector<mpz_class>(Acols));
    A_dec = pp.decryptMatrix(A_enc);
    cout << "Test Matrix Paillier passed" << endl;
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

void
testImagePaillierTime()
{
    struct timespec t0,t1;
    int n_iteration = 10;
    uint64_t t;

    //Generate Paillier key
    gmp_randstate_t randstate;
    gmp_randinit_default(randstate);
    gmp_randseed_ui(randstate,time(NULL));

    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t0);
    for (size_t i = 0; i < n_iteration; i++) {
        auto sk = Paillier_priv::keygen(randstate,32,2); //600 ,256
    }
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t1);
    t = (((uint64_t)t1.tv_sec) - ((uint64_t)t0.tv_sec) )* 1000000000 + (t1.tv_nsec - t0.tv_nsec);
    cerr << "Key Generation: "<<  ((double)t/1000000)/n_iteration <<"ms" << endl;
    auto sk = Paillier_priv::keygen(randstate,32,2); //600 ,256
    Paillier_priv pp(sk,randstate);

    auto pk = pp.pubkey();
    mpz_class n = pk[0];
    Paillier p(pk,randstate);

    //read image
    Mat src = imread( "cameraman.JPG", 0);
    imshow( "Original Image", src );

    int rows = src.rows;
    int cols = src.cols;
    vector < vector<mpz_class> > A(rows, vector<mpz_class>(cols));
    A = Mat2mpz(src);


    //Encrypt image
    vector < vector<mpz_class> > A_enc = p.encryptMatrix(A);
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t0);
    for (size_t i = 0; i < n_iteration; i++) {
        A_enc = p.encryptMatrix(A);
    }
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t1);
    //imshow( "Encrypted Image", mpz2Mat(A_enc) );
    t = (((uint64_t)t1.tv_sec) - ((uint64_t)t0.tv_sec) )* 1000000000 + (t1.tv_nsec - t0.tv_nsec);
    cerr << "public encryption: "<<  ((double)t/1000000)/n_iteration <<"ms" << endl;


    //Decrypt image
    vector < vector<mpz_class> > A_dec = pp.decryptMatrix(A_enc);
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t0);
    for (size_t i = 0; i < n_iteration; i++) {
        A_dec = pp.decryptMatrix(A_enc);
    }
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t1);
    t = (((uint64_t)t1.tv_sec) - ((uint64_t)t0.tv_sec) )* 1000000000 + (t1.tv_nsec - t0.tv_nsec);
    cerr << "private decryption: "<<  ((double)t/1000000)/n_iteration <<"ms" << endl;
    Mat dst;
    dst = mpz2Mat(A_dec);
    //imshow( "Decrypted Image", dst );


    //Negate image in plain domain
    dst = NegativeImage(src);
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t0);
    for (size_t i = 0; i < n_iteration; i++) {
        dst = NegativeImage(src);
    }
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t1);
    t = (((uint64_t)t1.tv_sec) - ((uint64_t)t0.tv_sec) )* 1000000000 + (t1.tv_nsec - t0.tv_nsec);
    cerr << "Negative Image PD: "<<  ((double)t/1000000)/n_iteration <<"ms" << endl;
    //imshow( "Negative Image", dst );

    //Negate image in encrypted domain
    vector < vector<mpz_class> > Neg_enc = p.NegativeImageH(A_enc);
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t0);
    for (size_t i = 0; i < n_iteration; i++) {
        vector < vector<mpz_class> > Neg_enc = p.NegativeImageH(A_enc);
    }
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t1);
    t = (((uint64_t)t1.tv_sec) - ((uint64_t)t0.tv_sec) )* 1000000000 + (t1.tv_nsec - t0.tv_nsec);
    cerr << "Negative Image ED: "<<  ((double)t/1000000)/n_iteration <<"ms" << endl;
    //imshow( "Encrypted Negative Image", mpz2Mat(Neg_enc) );

    vector < vector<mpz_class> > Neg_dec = pp.decryptMatrix(Neg_enc);
    //imshow( "Decrypted Negative Image", mpz2Mat(Neg_dec) );

    //Test Brightness
    int value = 50;
    //Increase brightness in plain domain
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t0);
    for (size_t i = 0; i < n_iteration; i++) {
        dst = AdjustBrightness(src, value);
    }
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t1);
    t = (((uint64_t)t1.tv_sec) - ((uint64_t)t0.tv_sec) )* 1000000000 + (t1.tv_nsec - t0.tv_nsec);
    cerr << "Adjust Brightness PD: "<<  ((double)t/1000000)/n_iteration <<"ms" << endl;
    dst = AdjustBrightness(src, value);
    //imshow( "Brightness Image", dst );

    //Increase brightness in encrypted domain
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t0);
    for (size_t i = 0; i < n_iteration; i++) {
        mpz_class value_enc = p.encrypt(value);
        vector < vector<mpz_class> > Bright_enc = p.AdjustBrightnessH(A_enc, value_enc);
    }
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t1);
    t = (((uint64_t)t1.tv_sec) - ((uint64_t)t0.tv_sec) )* 1000000000 + (t1.tv_nsec - t0.tv_nsec);
    cerr << "Adjust Brightness ED: "<<  ((double)t/1000000)/n_iteration <<"ms" << endl;
    mpz_class value_enc = p.encrypt(value);
    vector < vector<mpz_class> > Bright_enc = p.AdjustBrightnessH(A_enc, value_enc);
    //imshow( "Encrypted Brightness Image", mpz2Mat(Bright_enc) );

    vector < vector<mpz_class> > Bright_dec = pp.decryptMatrix(Bright_enc);
    //imshow( "Decrypted Brightness Image", mpz2Mat(Bright_dec) );


    //Test convolution
    vector < vector<mpz_class> > filter;
    //filter = {{1,0,-1},{2,0,-2},{1,0,-1}};//Vertical edges
    filter = {{1,2,1},{0,0,0},{-1,-2,-1}};  //Horizontal edges
    vector < vector<mpz_class> > conv_enc = p.convolutionH(A_enc, filter);
    //imshow( "Encrypted convoluted Image", mpz2Mat(conv_enc) );
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t0);
    for (size_t i = 0; i < n_iteration; i++) {
        vector < vector<mpz_class> > conv_enc = p.convolutionH(A_enc, filter);
    }
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t1);
    t = (((uint64_t)t1.tv_sec) - ((uint64_t)t0.tv_sec) )* 1000000000 + (t1.tv_nsec - t0.tv_nsec);
    cerr << "Convolution (3x3) Sobel ED: "<<  ((double)t/1000000)/n_iteration <<"ms" << endl;

    vector < vector<mpz_class> > conv_dec = pp.decryptMatrix(conv_enc);
    //imshow( "Decrypted convoluted Image", mpz2Mat(conv_dec) );

    //Test convolution in Plain domain
    dst = convolution(src, filter);
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t0);
    for (size_t i = 0; i < n_iteration; i++) {
        dst = convolution(src, filter);
    }
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t1);
    t = (((uint64_t)t1.tv_sec) - ((uint64_t)t0.tv_sec) )* 1000000000 + (t1.tv_nsec - t0.tv_nsec);
    cerr << "Convolution (3x3) Sobel PD: "<<  ((double)t/1000000)/n_iteration <<"ms" << endl;
    //imshow( "convoluted Image", dst );


    //Test Average filter
    filter = {{1,1,1},{1,1,1},{1,1,1}};
    vector < vector<mpz_class> > Avg_enc = p.convolutionH(A_enc, filter);
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t0);
    for (size_t i = 0; i < n_iteration; i++) {
        Avg_enc = p.convolutionH(A_enc, filter);
    }
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t1);
    t = (((uint64_t)t1.tv_sec) - ((uint64_t)t0.tv_sec) )* 1000000000 + (t1.tv_nsec - t0.tv_nsec);
    cerr << "Average (3x3) ED: "<<  ((double)t/1000000)/n_iteration <<"ms" << endl;
    //imshow( "Encrypted averaged Image", mpz2Mat(Avg_enc) );
    vector < vector<mpz_class> > Avg_dec = pp.decryptMatrix(Avg_enc);
    //Divide by filter summation in plain domain
    vector < vector<mpz_class> > Avg_decP = divideP(Avg_dec, filter);
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t0);
    for (size_t i = 0; i < n_iteration; i++) {
        Avg_decP = divideP(Avg_dec, filter);
    }
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t1);
    t = (((uint64_t)t1.tv_sec) - ((uint64_t)t0.tv_sec) )* 1000000000 + (t1.tv_nsec - t0.tv_nsec);
    cerr << "Average (3x3) ED Post Processing: "<<  ((double)t/1000000)/n_iteration <<"ms" << endl;
    //imshow( "Decrypted averaged Image", mpz2Mat(Avg_decP) );


    //Test average filter in plain domain
    dst = AverageFilter(src,filter);
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t0);
    for (size_t i = 0; i < n_iteration; i++) {
        dst = AverageFilter(src,filter);
    }
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t1);
    t = (((uint64_t)t1.tv_sec) - ((uint64_t)t0.tv_sec) )* 1000000000 + (t1.tv_nsec - t0.tv_nsec);
    cerr << "Average (3x3) PD: "<<  ((double)t/1000000)/n_iteration <<"ms" << endl;
    //imshow( "Averaged Image", dst );

}

void
testImagePaillier()
{
    //Generate Paillier key
    gmp_randstate_t randstate;
    gmp_randinit_default(randstate);
    gmp_randseed_ui(randstate,time(NULL));

    auto sk = Paillier_priv::keygen(randstate,32,2); //600 ,256
    Paillier_priv pp(sk,randstate);

    auto pk = pp.pubkey();
    mpz_class n = pk[0];
    Paillier p(pk,randstate);
    cout << "key  ready ...\n" << endl;

    //read image
    Mat src = imread( "cameraman.JPG", 0);
    imshow( "Original Image", src );

    int rows = src.rows;
    int cols = src.cols;
    vector < vector<mpz_class> > A(rows, vector<mpz_class>(cols));
    A = Mat2mpz(src);

    Img img1(A);
    img1.setPbKey(p);

    //Encrypt image
    vector < vector<mpz_class> > A_enc = img1.encryptMatrix(A);
    imshow( "Encrypted Image", mpz2Mat(A_enc) );

    //The client only can decrypt a matrix
    vector < vector<mpz_class> > A_dec = pp.decryptMatrix(A_enc);

    //Decrypt image
    Mat dst;
    dst = mpz2Mat(A_dec);
    imshow( "Decrypted Image", dst );
    //cout << "Root Mean Square Error in decryption: "<<getRMSE(src,dst)<<endl;

    //Negate image in plain domain
    dst = NegativeImage(src);
    imshow( "Negative Image", dst );

    //Negate image in encrypted domain
    vector < vector<mpz_class> > Neg_enc = img1.NegativeImageH(A_enc);
    imshow( "Encrypted Negative Image", mpz2Mat(Neg_enc) );

    vector < vector<mpz_class> > Neg_dec = pp.decryptMatrix(Neg_enc);
    imshow( "Decrypted Negative Image", mpz2Mat(Neg_dec) );
    //cout << "Root Mean Square Error in Negative image: "<<getRMSE(dst,mpz2Mat(Neg_dec))<<endl;

    //Test Brightness
    int value = 50;
    //Increase brightness in plain domain
    dst = AdjustBrightness(src, value);
    imshow( "Brightness Image", dst );

    //Increase brightness in encrypted domain
    mpz_class value_enc = p.encrypt(value);
    vector < vector<mpz_class> > Bright_enc = img1.AdjustBrightnessH(A_enc, value_enc);
    imshow( "Encrypted Brightness Image", mpz2Mat(Bright_enc) );

    vector < vector<mpz_class> > Bright_dec = pp.decryptMatrix(Bright_enc);
    imshow( "Decrypted Brightness Image", mpz2Mat(Bright_dec) );
    //cout << "Root Mean Square Error in Brightness image: "<<getRMSE(dst,mpz2Mat(Bright_dec))<<endl;
    //post-processing to fix the image
    //Bright_dec = FixRange(Bright_dec);
    //imshow( "Decrypted Brightness Image", mpz2Mat(Bright_dec) );

    //Test convolution
    vector < vector<mpz_class> > filter;
    //filter = {{1,0,-1},{2,0,-2},{1,0,-1}};//Vertical edges
    filter = {{1,2,1},{0,0,0},{-1,-2,-1}};  //Horizontal edges

    vector < vector<mpz_class> > conv_enc = img1.convolutionH(A_enc, filter);
    imshow( "Encrypted convoluted Image", mpz2Mat(conv_enc) );
    vector < vector<mpz_class> > conv_dec = pp.decryptMatrix(conv_enc);
    imshow( "Decrypted convoluted Image", mpz2Mat(conv_dec) );

    //Test convolution in Plain domain
    dst = convolution(src, filter);
    imshow( "convoluted Image", dst );
    //cout << "Root Mean Square Error in convoluted image: "<<getRMSE(dst,mpz2Mat(conv_dec))<<endl;

    //Test Average filter
    filter = {{1,1,1},{1,1,1},{1,1,1}};
    vector < vector<mpz_class> > Avg_enc = img1.convolutionH(A_enc, filter);
    imshow( "Encrypted averaged Image", mpz2Mat(Avg_enc) );
    vector < vector<mpz_class> > Avg_dec = pp.decryptMatrix(Avg_enc);
    //Divide by filter summation in plain domain
    vector < vector<mpz_class> > Avg_decP = divideP(Avg_dec, filter);
    imshow( "Decrypted averaged Image", mpz2Mat(Avg_decP) );

    //Test average filter in plain domain
    dst = AverageFilter(src,filter);
    imshow( "Averaged Image", dst );
    //cout << "Root Mean Square Error in Averaged image: "<<getRMSE(dst,mpz2Mat(Avg_decP))<<endl;

    //wait for a key press infinitely
    waitKey(0);
}

int
main ( int argc, char **argv )
{
    //SetSeed(to_ZZ(time(NULL)));
    //test_paillier();
    //test_image();

    testImagePaillier(); /*For generating images*/
    //testImagePaillierTime(); /*For calculating time*/

    return 0;
}
