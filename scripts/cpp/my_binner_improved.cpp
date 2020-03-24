#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <cctype>
#include <locale>
#include <sstream>

using namespace std;

#define NUMBINS 5

int main(int argc, char** argv)
{

    // Check that a filename was given
    if ((argv[1] == NULL) || (argv[2] == NULL) || (argv[3] == NULL)){
        cout << "Usage: " << argv[0] << " <x_coord file> <y_coord file> <qa file>" << endl;
        exit(1);
    }

    // First calculate min and max of data
    int max = INT_MIN;
    int min = INT_MAX; 

    ifstream is(argv[2]);
    string str;
    while(getline(is, str))
    {
        int i_str = std::stoi(str);
        if (i_str > max){
            max = i_str;
        }
        if (i_str < min){
            min = i_str;
        }
    }
    cout << "Binning will be done based on " << argv[2] << endl;
    cout << "Max of " << argv[2] << ": " << max << endl;
    cout << "Min of " << argv[2] << ": " << min << endl;

    ifstream is1(argv[1]);
    ifstream is2(argv[2]);
    ifstream is3(argv[3]);
    string x_coord;
    string y_coord;
    string qa;


    ofstream outFileX[NUMBINS];
    ofstream outFileY[NUMBINS];
    ofstream outFileQA[NUMBINS];

    stringstream sstm;

    for (int i=0; i<NUMBINS; i++){
        sstm.str("");
        sstm << "bin" << i+1 << "_x_coord.txt";
        outFileX[i].open(sstm.str());

        sstm.str("");
        sstm << "bin" << i+1 << "_y_coord.txt";
        outFileY[i].open(sstm.str());

        sstm.str("");
        sstm << "bin" << i+1 << "_qa.txt";
        outFileQA[i].open(sstm.str());

    }

    // Bins are [2072, 41846, 81619, 121392, 161165, 200941]
    // Create bins

    int bin_size = ((max-min)/NUMBINS) + NUMBINS; // Not sure if this is right,
                                                  // Adding NUMBINS to correct for int round down

    int bins [NUMBINS+1];
    sstm.str("Bins are: ");
    for (int i = 0; i < NUMBINS+1; i++){
        bins[i] = min + i*bin_size;
        if (i != 0){
            sstm << "[ " << bins[i-1] << ", " << bins[i] << "], ";
        }
    } 
    cout << sstm.str() << endl;

    //int bins [NUMBINS+1] = {2072, 41846, 81619, 121392, 161165, 200941};


    while((getline(is1, x_coord) && (getline(is2, y_coord)) && (getline(is3, qa))))
    {
        //cout << x_coord << "  " << y_coord << endl;

        int i_x_coord = std::stoi(x_coord);
        int i_y_coord = std::stoi(y_coord); // Probably bad form to declare in a loop
        int i_qa      = std::stoi(qa);

        // Bins are [2072, 41846, 81619, 121392, 161165, 200941]


        for (int i = 0; i < NUMBINS; i++){    
                if ((i_y_coord >= bins[i]) && (i_y_coord < bins[i+1])){
                    outFileX[i] << i_x_coord << endl;
                    outFileY[i] << i_y_coord << endl;
                    outFileQA[i] << i_qa << endl;
                    break;
                }
        }

        /*
        if ((i_y_coord > 2072) && (i_y_coord <= 41846)){
            bin1x << i_x_coord << endl;
            bin1y << i_y_coord << endl;
            bin1qa << i_qa << endl;
        }
        if ((i_y_coord > 41846) && (i_y_coord <= 81619)){
            bin2x << i_x_coord << endl;
            bin2y << i_y_coord << endl;
            bin2qa << i_qa << endl;
        }
        if ((i_y_coord > 81619) && (i_y_coord <= 121392)){
            bin3x << i_x_coord << endl;
            bin3y << i_y_coord << endl;
            bin3qa << i_qa << endl;
        }
        if ((i_y_coord > 121392) && (i_y_coord <= 161165)){
            bin4x << i_x_coord << endl;
            bin4y << i_y_coord << endl;
            bin4qa << i_qa << endl;
        }
        if ((i_y_coord > 161165) && (i_y_coord <= 200941)){
            bin5x << i_x_coord << endl;
            bin5y << i_y_coord << endl;
            bin5qa << i_qa << endl;
        }*/

    }

    for (int i=0; i<NUMBINS; i++){
        outFileX[i].close();
        outFileY[i].close();
        outFileQA[i].close();
    }

    /*
    bin1x.close();
    bin1y.close();
    bin1qa.close();
    bin2x.close();
    bin2y.close();
    bin2qa.close();
    bin3x.close();
    bin3y.close();
    bin3qa.close();
    bin4x.close();
    bin4y.close();
    bin4qa.close();
    bin5x.close();
    bin5y.close();
    bin5qa.close();*/

}













