#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <cctype>
#include <locale>

using namespace std;

int main()
{
    ifstream is("read_numbers.txt");
    int x = 0;
    string str;
    while(getline(is, str))
    {
    	x++;
    	int i_str = std::stoi(str);
    	if (x != i_str){
    		cout << str << " not equal " << x << endl;
    	}
    }
}




