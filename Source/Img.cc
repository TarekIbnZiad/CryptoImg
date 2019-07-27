/*
 * File:   Img.cc
 * Author: Mohamed TarekIbnZiad
 *
 * Created on September 16, 2015, 4:12 AM
 */

#include <cstdlib>
#include <Img.hh>
using namespace std;


Img::Img(vector < vector<encnum> > &_A){
    A = _A;
}


void
Img::setPbKey( Paillier &_p)
{
    p = &_p;
}

//more efficient implementation
vector < vector<encnum> >
Img::NegativeImageH(const vector < vector<encnum> > &src_enc)
{
    int rows = src_enc.size();
    int cols = src_enc.data()->size();

    double white = 255;
    encnum white_enc = p->encrypt_f(white); //encrypt it one time only instead of a matrix encryption
    vector < vector<encnum> > dst_enc (rows, vector<encnum>(cols));
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++){
            dst_enc[i][j] = p->sub_f(white_enc, src_enc[i][j]);
        }
    }
    return dst_enc;
}

vector < vector<encnum> >
Img::AdjustBrightnessH(const vector < vector<encnum> > &src_enc, const encnum &value_enc)
{
    int rows = src_enc.size();
    int cols = src_enc.data()->size();

    vector < vector<encnum> > dst_enc (rows, vector<encnum>(cols));
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++){
            dst_enc[i][j] = p->add_f(value_enc, src_enc[i][j]);
        }
    }
    return dst_enc;
}

vector < vector<encnum> >
Img::convolutionH(const vector < vector<encnum> > &src_enc, const vector < vector<double> > &filter)
{
    //Src image size
    int rows = src_enc.size();
    int cols = src_enc.data()->size();
    int frows = filter.size();
    int fcols = filter.data()->size();

    vector < vector<encnum> > dst_enc (rows, vector<encnum>(cols));

    encnum zero_enc = p->encrypt_f(0);

    for (int i=0; i < rows; i++){
	for (int j=0; j < cols; j++){
            if( (i < frows/2) || (i >= rows - frows/2) || (j < fcols/2) || (j >= cols - fcols/2) ) {
                //For the borders --> Keep them as in the src
                dst_enc[i][j] = src_enc[i][j];
            }
            else{
                //For the rest of image
                encnum sum = zero_enc;
                for(int m=-frows/2; m<=frows/2; m++){
                    for(int n=-fcols/2; n<=fcols/2; n++){
                        if (filter[m+frows/2][n+fcols/2] == 0) { //Not encrypted value
                            continue;
                        }
                        sum = p->add_f(sum, p->constMult_f(filter[m+frows/2][n+fcols/2],src_enc[i+m][j+n]));
                    }
                }
                dst_enc[i][j] = sum;
            }
        }
    }
    return dst_enc;
}

/*
 * Adds the values under the specified element mask
 */
vector < vector<encnum> >
Img::morphH(const vector < vector<encnum> > &src_enc, const vector < vector<double> > &element)
{
    //Src image size
    int rows = src_enc.size();
    int cols = src_enc.data()->size();
    int erows = element.size();
    int ecols = element.data()->size();

    vector < vector<encnum> > dst_enc (rows, vector<encnum>(cols));

    encnum zero_enc = p->encrypt_f(0);

    for (int i=erows/2; i < rows - erows/2; i++){
	for (int j=ecols/2; j < cols - ecols/2; j++){
            encnum sum = zero_enc;
            for(int m=-erows/2; m<=erows/2; m++){
		for(int n=-ecols/2; n<=ecols/2; n++){
                    sum = p->add_f(sum, src_enc[i+m][j+n]);
                    //sum = p->add_f(sum, p->constMult_f(1.0/255.0,src_enc[i+m][j+n]));
                }
            }
            dst_enc[i][j] = sum;
        }
    }
    return dst_enc;
}

vector < vector<encnum> >
Img::equalHistH(const vector < vector<encnum> > &src_enc, int &rows, int &cols)
{
    int G = src_enc.size(); //Intensity levels

    vector <encnum> cumulativeHist (G);
    cumulativeHist[0] = src_enc[0][1];
    for (int i = 1; i < G; i++) {
        cumulativeHist[i] = p->add_f(cumulativeHist[i-1], src_enc[i][1]);
    }

    double alpha = (1.0*G-1)/(rows*cols);
    //encrypted as 0.001 so we had to increase the precision to 0.000001 and key to 128 !
    vector < vector<encnum> > newColorFreq_enc(G, vector<encnum> (2));
    for (int i = 0; i < G; i++) {
        newColorFreq_enc[i][1] = p->constMult_f(alpha, cumulativeHist[i]);
        newColorFreq_enc[i][0] = src_enc[i][0]; //The enc intensity level
    }
    return newColorFreq_enc;
}

/*
 * Apply the high boost filter
 * output = A*Img - AvgImg
 */
vector < vector<encnum> >
Img::highBoostFilterH(const vector < vector<encnum> > &src_enc, double &A)
{
    int rows = src_enc.size();
    int cols = src_enc.data()->size();

    vector < vector<encnum> > dst_enc (rows, vector<encnum>(cols));

    if (A < 0){
        cout << "A should not be a negative value" << endl;
        return src_enc;
    }

    //Create box filter
    vector < vector <double> > filter = {{1.0/9,1.0/9,1.0/9},{1.0/9,1.0/9,1.0/9},{1.0/9,1.0/9,1.0/9}};

    vector < vector <encnum> > AvgImg = convolutionH(src_enc, filter);

    //Multiply source by A and subtract the average
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++){
            dst_enc[i][j] = p->sub_f(p->constMult_f(A,src_enc[i][j]), AvgImg[i][j]);
        }
    }
    return dst_enc;
}

vector < vector<encnum> >
Img::highBoostFilterH2(const vector < vector<encnum> > &src_enc, const vector < vector<encnum> > &Avg_enc, double &A)
{
    int rows = src_enc.size();
    int cols = src_enc.data()->size();

    vector < vector<encnum> > dst_enc (rows, vector<encnum>(cols));

    if (A < 0){
        cout << "A should not be a negative value" << endl;
        return src_enc;
    }

    //Multiply source by A and subtract the average
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++){
            dst_enc[i][j] = p->sub_f(p->constMult_f(A,src_enc[i][j]), Avg_enc[i][j]);
        }
    }
    return dst_enc;
}

/*
 * Apply Edge Detection Filter in encrypted domain using one of the following methods
 * 'P' Prewitt, 'S' Sobel, 'R' Robinson, 'K' Kirsch
 */
vector < vector<encnum> >
Img::EdgeDetectionFilterH(const vector < vector<encnum> > &src_enc, char &type, char &direction)
{
    int rows = src_enc.size();
    int cols = src_enc.data()->size();

    vector < vector <double> >  dx(3, vector <double> (3));
    vector < vector <double> >  dy(3, vector <double> (3));

    vector < vector<encnum> > dst_enc (rows, vector<encnum>(cols));

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
            return dst_enc;
    }

    switch (direction){
        case 'H':
            dst_enc = convolutionH(src_enc, dx);
            break;
        case 'V':
            dst_enc = convolutionH(src_enc, dy);
            break;
        default:
            cout<<"Please enter a valid direction" <<endl;
            return dst_enc;
    }
    return dst_enc;
}
