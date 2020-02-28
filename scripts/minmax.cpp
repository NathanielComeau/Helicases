#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <cctype>
#include <locale>

using namespace std;

int main()
{
    int max = INT_MIN;
    int min = INT_MAX; 

    ifstream is("tile_numbers.txt");
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
}




