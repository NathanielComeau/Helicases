#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <cctype>
#include <locale>

using namespace std;


// Get frequencies of data. Assumes range of data is fairly small; this will work poorly for
// data with max-min > than roughly a few billion.

int main(int argc, char** argv)
{

    // Check that a filename was given
    if (argv[1] == NULL){
        cout << "Usage: " << argv[0] << " <input file>" << endl;
        exit(1);
    }



    // First calculate min and max of data
    int max = INT_MIN;
    int min = INT_MAX; 

    ifstream is(argv[1]);
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
    cout << "Final Max: " << max << endl;
    cout << "Final Min: " << min << endl;

    // Make an array to hold each possible datapoint
    int size = max-min+1;
    int *frequencies;
    frequencies = (int*) malloc((max-min+1)*sizeof(int));
    for (int i=0; i<(max-min+1); i++){
        frequencies[i] = 0;
    }

    ifstream is2(argv[1]);
    while(getline(is2, str))
    {
        int i_str = std::stoi(str);
        frequencies[i_str-min] += 1;
    }

    cout << "i  freq" << endl;
    for(int i = 0; i< (max-min+1); i++){
        if (frequencies[i] != 0) {
            cout << (i+min) << ": " << frequencies[i] << endl;
        }
    }

    free(frequencies);

}




