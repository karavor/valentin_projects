#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>

using namespace std;

const char* AAGS_output = "./output/AAGS_result";
const char* AAGS_caller_stat = "./output/AAGS_caller_stat";

const int max_loop_size = 50;

#define AAGS_CALLER "\
    #/bin/bash \n\
    ./antiaff_antihairpin_gen_smart \n\
"

int int_pow (int n, int power)
{
    int result = 1;

    for (int i = 0; i < power; i++)
        result *= n;

    return result;
}

int int_numb_reader (string* s)
{
    int result = 0;
    int ten_pow = 0;
    int minus_flag = 1;

    if ( (*s)[0] == '-' )
    {
        minus_flag = -1;
        (*s).erase( (*s).begin(), (*s).begin() + 1 );
    }

    for (int i = (*s).length() - 1; i >= 0 ; i--)
    {
        result += ( (int)(*s)[i] - 48 ) * int_pow (10, ten_pow);
        ten_pow++;
    }
    return result * minus_flag;
}

float decimal_numb_reader( string* s )
{
    int minus_flag = 1;
    if ( (*s)[0] == '-' )
    {
        minus_flag = -1;
        (*s).erase( (*s).begin(), (*s).begin() + 1 );
    }
    size_t dot_position = (*s).find('.');
    if ( dot_position == (*s).npos )
    {
        return int_numb_reader(s) * minus_flag;
    }
    string int_part;

    for (int i = 0; i < (int)dot_position ; i++)
        int_part.push_back( (*s)[i] );

    (*s).erase( (*s).begin(), (*s).begin() + dot_position + 1 ); //in s now is real part

    return (float) minus_flag * ( (float)int_numb_reader(&int_part) + (float)int_numb_reader(s) / (float)int_pow (10, (*s).length() ) );
}


main()
{
    ofstream write;
    write.open(AAGS_caller_stat);

    ifstream read;
    string read_str;

    float dG, max_dG = -10000;
    string best_variant;
    string variant;

    for (int i = 0; i < max_loop_size; i++)
    {
        system(AAGS_CALLER);

        read.open("./error_flag");
        getline(read, read_str);
        read.close();

        if ( read_str[0] == '1' )
            continue;

        read.open(AAGS_output);

        getline(read, read_str);
        write << read_str << endl;

        dG = decimal_numb_reader(&read_str);

        getline(read, read_str);
        write << read_str << endl;

        variant = read_str;

        do
        {
            getline(read, read_str);
            write << read_str << endl;
        }
        while (read_str[0] != '!');

        if ( dG > max_dG && read_str[1] == '+' )
        {
            max_dG = dG;
            best_variant = variant;
        }

        read.close();
    }
    write << endl << endl << "best_variant: " << best_variant;
    write.close();
    return 0;
}
