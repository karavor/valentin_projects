#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <cmath>
#include <vector>

using namespace std;

const char* input_file    = "./output/OS_result";
const char* output_file    = "./output/best_olig";

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
        result += ( (int)(*s)[i] - (int)'0' ) * int_pow (10, ten_pow);
        ten_pow++;
    }
    return result * minus_flag;
}

float decimal_numb_reader (string* s)
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

float exp_numb_reader (string* s)
{
    size_t exp_position = (*s).find('e');

    if ( exp_position == (*s).npos )
    {
        exp_position = (*s).find('E');
        if ( exp_position == (*s).npos )
            return decimal_numb_reader(s);
    }
    int minus_exp_flag = 1;

    if ( (*s)[exp_position + 1]  == '-' )
        minus_exp_flag = -1;

    string mantissa;

    for (int i = 0; i < (int)exp_position ; i++)
        mantissa.push_back( (*s)[i] );

    (*s).erase( (*s).begin(), (*s).begin() + exp_position + 1 );

    if ( (*s)[0] == '+' || (*s)[0] == '-')
        (*s).erase( (*s).begin(), (*s).begin() + 1 );

    int exp_part = int_numb_reader(s);

    return decimal_numb_reader(&mantissa) * pow(10.0, exp_part * minus_exp_flag);
}

bool input_data_reader(float*  crit_dG,
                       string* left_olig_part,
                       string* right_olig_part,
                       vector<string>* antiaff_oligs)
{
    ifstream read;
    read.open(input_file);

    getline(read, (*left_olig_part));
    (*left_olig_part).erase( (*left_olig_part).begin(),
                             (*left_olig_part).begin() +
                             (*left_olig_part).find(' ') + 1);

    (*crit_dG) = decimal_numb_reader(left_olig_part);

    getline(read, (*left_olig_part));
    (*left_olig_part).erase( (*left_olig_part).begin(),
                             (*left_olig_part).begin() +
                             (*left_olig_part).find(' ') + 1);

    getline(read, (*right_olig_part));
    (*right_olig_part).erase( (*right_olig_part).begin(),
                              (*right_olig_part).begin() +
                              (*right_olig_part).find(' ') + 1);

    int i = 0;
    string read_str;

    while(!read.eof())
    {
        i++;
        getline(read, read_str);
        (*antiaff_oligs).push_back(read_str);

        ( (*antiaff_oligs)[i] ).erase( ( (*antiaff_oligs)[i] ).begin(),
                                       ( (*antiaff_oligs)[i] ).begin() +
                                       ( (*antiaff_oligs)[i] ).find(' ') + 1);
    }
    read.close();
}

void string_with_prmtr_reader(ifstream* read, float* prmtr)
{
    string read_str;
    getline( (*read), read_str);
    read_str.erase( read_str.begin(), read_str.begin() + read_str.find( ' ' ) + 1 );
    *prmtr = decimal_numb_reader(&read_str);
}

float arr_sum (float* arr, int arr_length)
{
    float result = 0;

    for (int i = 0; i < arr_length; i++)
        result += arr[i];

    return result;
}

int main()
{
    ifstream read;
    read.open(input_file);
    float prmtrs[9];
    float best_prmtrs_sum[4] = {1000, 1000, 1000, 1000},
          prmtrs_sum[4],
          defined_dG[4] = {-0.2, -0.3, -0.4, -0.5},
          dG;

    string read_str;
    string olig,
           best_olig[4];

    while(!read.eof())
    {
        getline( read, read_str);
        getline( read, read_str);

        olig = read_str;
        olig.erase( olig.begin() + olig.find( ' ' ), olig.begin() + olig.length() );

        read_str.erase( read_str.begin(), read_str.begin() + read_str.find( ' ' ) + 1 );
        dG = decimal_numb_reader(&read_str);

        getline( read, read_str);

        for (int i = 0; i < 5; i++)
            string_with_prmtr_reader( &read, &(prmtrs[i]) );

        getline( read, read_str);

        string_with_prmtr_reader( &read, &(prmtrs[5]) );
        string_with_prmtr_reader( &read, &(prmtrs[6]) );

        getline( read, read_str);

        string_with_prmtr_reader( &read, &(prmtrs[7]) );
        string_with_prmtr_reader( &read, &(prmtrs[8]) );

        for (int i = 0; i < 4; i++)
        {
            if (dG == defined_dG[i])
            {
                prmtrs_sum[i] = arr_sum(prmtrs, 9);
                if (prmtrs_sum[i] < best_prmtrs_sum[i])
                {
                    best_prmtrs_sum[i] = prmtrs_sum[i];
                    best_olig[i] = olig;
                }
            }
        }
    }

    read.close();

    ofstream write;
    write.open(output_file);
    for (int i = 0; i < 4; i++)
        write << best_olig[i] << ' ' << -0.2 - 0.1*i << endl;
    write.close();
    return 0;
}
