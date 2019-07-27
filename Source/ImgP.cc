#include <ImgP.hh>
using namespace std;


ImgP::ImgP(){

}


/*
 * Apply average filter in Plain domain using open cv ready made function (blur)
 */
Mat
ImgP::AverageFilterOpenCV(Mat src)
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
ImgP::AverageFilter(Mat src, vector < vector<double> > &filter)
{
    Mat dst = src.clone(); //Copy of src as initial value
    Scalar s;
    //Filter size
    int ro = filter.size();
    int co = filter.data()->size();

    //get summation of the filter elements to divide by
    double nFilter = 0;
    for (int i=0; i < ro; i++){
	     for (int j=0; j < co; j++){
            nFilter += filter[i][j];
       }
    }

    for (int i=1; i < src.rows - 2; i++){
	     for (int j=1; j < src.cols - 2; j++){
            double sum=0;
            for(int m=-1; m<=1; m++){
		           for(int n=-1; n<=1; n++){
                    s = src.at<uchar>(i+m,j+n);
                    sum += s.val[0] * filter[m+1][n+1];
               }
            }
            s.val[0] = sum/nFilter;
            dst.at<uchar>(i,j) = s.val[0];
        }
    }
    return dst;
}

/*
 * Apply Sobel operator in Plain domain using open cv ready made function (Sobel)
 */
Mat
ImgP::SobelFilterOpenCV(Mat src)
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
ImgP::SobelFilter(Mat src)
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
ImgP::convolution(Mat src, vector < vector<double> > &filter)
{
    //Assuming all sizes are odd
    int rows = src.rows;
    int cols = src.cols;
    int frows = filter.size();
    int fcols = filter.data()->size();

    Mat dst = src.clone(); //Copy of src as initial value
    Scalar s;

    for (int i=frows/2; i < rows - frows/2; i++){
	     for (int j=fcols/2; j < cols - fcols/2; j++){
            double sum=0;
            for(int m=-frows/2; m<=frows/2; m++){
		            for(int n=-fcols/2; n<=fcols/2; n++){
                    if (filter[m+frows/2][n+fcols/2] == 0) { //Not encrypted value
                        continue;
                    }
                    s = src.at<uchar>(i+m,j+n);
                    sum += s.val[0] * filter[m+frows/2][n+fcols/2];
                }
            }
            s.val[0] = abs(sum);
            dst.at<uchar>(i,j) = s.val[0];
        }
    }
    return dst;
}

/*
 * used to Add/Sub constant value to the image brightness in plain domain
 */
Mat
ImgP::AdjustBrightness(Mat src, int value)
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
 * used to get the negative image in plain domain
 */
Mat
ImgP::NegativeImage(Mat src)
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
ImgP::AddSaltPepperNoise(Mat src)
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
 * used to perform Histogram equalization in plain domain
 * TODO: Need to be optimized
 */
Mat
ImgP::equalHist(Mat src)
{
    int rows = src.rows;
    int cols = src.cols;
    Scalar s;
    vector <int> Frequencies (256);
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++){
            s = src.at<uchar>(i,j);
            Frequencies[(int)s[0]]++;
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

    vector <int> cumulativeHist (G);
    cumulativeHist[0] = ColorFreq[0][1];
    for (int i = 1; i < G; i++) {
        cumulativeHist[i] = cumulativeHist[i-1] + ColorFreq[i][1];
    }

    double alpha = (1.0*G-1)/(rows*cols);
    vector < vector<double> > newFreqHist(G, vector<double> (2));
    for (int i = 0; i < G; i++) {
        newFreqHist[i][1] = alpha * cumulativeHist[i];
        newFreqHist[i][0] = ColorFreq[i][0]; //The enc intensity level
    }

    //updateHist
    vector <int> newColors (256);
    for (int i = 0; i<256; i++)
        newColors[i] = i;

    for (int i = 0; i<G; i++){
        newColors[newFreqHist[i][0]] = round(newFreqHist[i][1]);
    }

    Mat dst(rows,cols,CV_8UC1);
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++){
            s = src.at<uchar>(i,j);
            dst.at<uchar>(i,j) = newColors[(int)s[0]];
        }
    }
    return dst;
}

/*
 * Apply the high boost filter
 * output = A*Img - AvgImg
 */
Mat
ImgP::highBoostFilter(Mat src, double &A)
{
    int rows = src.rows;
    int cols = src.cols;
    Scalar s, s_avg;

    Mat dst(rows,cols,CV_8UC1);
    if (A < 0){
        cout << "A should not be a negative value" << endl;
        return src;
    }

    //Create box filter
    vector < vector <double> > filter = {{1.0/9,1.0/9,1.0/9},{1.0/9,1.0/9,1.0/9},{1.0/9,1.0/9,1.0/9}};
    //Mat AvgImg = AverageFilter(src,filter);
    //Manual Average to conserve the FP part of the number
    vector < vector <double> > AvgImg_d(rows, vector <double>(cols));
    for (int i=0; i < rows ; i++){
	     for (int j=0; j < cols ; j++){
            if( (i < 1) || (i >= rows - 1) || (j < 1) || (j >= cols - 1) ) {
                //For the borders --> Keep them as in the src
                s = src.at<uchar>(i,j);
                AvgImg_d[i][j] = (int)s[0];
            }
            else{
                double sum=0;
                for(int m=-1; m<=1; m++){
                    for(int n=-1; n<=1; n++){
                        s = src.at<uchar>(i+m,j+n);
                        sum += s.val[0] * filter[m+1][n+1];
                    }
                }
                AvgImg_d[i][j] = sum/9.0;
            }
        }
    }

    double temp = 0;
    //Multiply source by A and subtract the average
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++){
            s = src.at<uchar>(i,j);
            //s_avg = AvgImg.at<uchar>(i,j);
            //temp = A*(int)s[0] - (int)s_avg[0];
            temp = A*(int)s[0] - AvgImg_d[i][j];
            dst.at<uchar>(i,j) = (int)abs(temp);
        }
    }
    return dst;
}

/*
 * Apply Edge Detection Filter in Plain domain using one of the following methods
 * 'P' Prewitt, 'S' Sobel, 'R' Robinson, 'K' Kirsch
 * TODO: 'R' Roberts 2x2 filter
 */
Mat
ImgP::EdgeDetectionFilter(Mat src, char &type, char &direction)
{
    vector < vector <double> >  dx(3, vector <double> (3));
    vector < vector <double> >  dy(3, vector <double> (3));

    Mat dst = src.clone();

    switch (type){
        case 'P':
            dx = {{-1,0,1},{-1,0,1},{-1,0,1}};
            dy = {{1,1,1},{0,0,0},{-1,-1,-1}};
            break;
        case 'S':
            dx = {{-1,0,1},{-2,0,2},{-1,0,1}};
            dy = {{1,2,1},{0,0,0},{-1,-2,-1}};
            break;
        case 'R':
            dx = {{-1,1,1},{-1,-2,1},{-1,1,1}};
            dy = {{1,1,1},{1,-2,1},{-1,-1,-1}};
            break;
        case 'K':
            dx = {{-5,3,3},{-5,0,3},{-5,3,3}};
            dy = {{3,3,3},{3,0,3},{-5,-5,-5}};
            break;
        default:
            cout<<"Please enter a valid type" <<endl;
            return dst;
    }

    switch (direction){
        case 'H':
            dst = convolution(src, dx);
            break;
        case 'V':
            dst = convolution(src, dy);
            break;
        default:
            cout<<"Please enter a valid direction" <<endl;
            return dst;
    }
    return dst;
}
