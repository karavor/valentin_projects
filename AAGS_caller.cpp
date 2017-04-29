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

void string_with_prmtr_reader(ifstream* read, float* prmtr)
{
    string read_str;
    getline( (*read), read_str);
    read_str.erase( read_str.begin(), read_str.begin() + read_str.find( ' ' ) + 1 );
    *prmtr = decimal_numb_reader(&read_str);
}

main()
{
    ofstream write;
    write.open(AAGS_caller_stat);

    float crit_dG,
          crit_c_dG,
          crit_cmplx_dG,
          crit_cross_aff,
          crit_c_cross_aff,
          reading_prmtr;

    string read_str;

    ifstream read;
    read.open("./critical_prmtrs");

    string_with_prmtr_reader(&read, &crit_dG);
    string_with_prmtr_reader(&read, &crit_c_dG);
    string_with_prmtr_reader(&read, &crit_cmplx_dG);
    string_with_prmtr_reader(&read, &crit_cross_aff);
    string_with_prmtr_reader(&read, &crit_c_cross_aff);

    read.close();

    int continue_flag;
    int pause;

    for (int i = 0; i < max_loop_size; i++)
    {
        continue_flag = 0;

        system(AAGS_CALLER);

        read.open(AAGS_output);

        string_with_prmtr_reader(&read, &reading_prmtr);
        if( reading_prmtr < crit_dG )
        {
            read.close();
            continue;
        }

        string_with_prmtr_reader(&read, &reading_prmtr);
        if( reading_prmtr < crit_c_dG )
        {
            read.close();
            continue;
        }

        string_with_prmtr_reader(&read, &reading_prmtr);
        if( reading_prmtr < crit_c_dG )
        {
            read.close();
            continue;
        }

        string_with_prmtr_reader(&read, &reading_prmtr);
        if( reading_prmtr > crit_cmplx_dG )
        {
            read.close();
            continue;
        }

        string_with_prmtr_reader(&read, &reading_prmtr);
        if( reading_prmtr > crit_cmplx_dG )
        {
            read.close();
            continue;
        }

        string_with_prmtr_reader(&read, &reading_prmtr);
        if( reading_prmtr > crit_c_cross_aff )
        {
            read.close();
            continue;
        }

        string_with_prmtr_reader(&read, &reading_prmtr);
        if( reading_prmtr > crit_c_cross_aff )
        {
            read.close();
            continue;
        }

        while (1)
        {
            getline(read, read_str);
            if (read_str[0] == '$')
                break;

            read_str.erase( read_str.begin(), read_str.begin() + read_str.find( ' ' ) + 1 );
            reading_prmtr = decimal_numb_reader(&read_str);

            if( reading_prmtr > crit_cross_aff )
            {
                continue_flag = 1;
                break;
            }
        }

        if(continue_flag)
        {
            read.close();
            continue;
        }

        read_str.erase(read_str.begin());
        write << read_str << ' ';

        getline(read, read_str);
        read_str.erase( read_str.begin(), read_str.begin() + read_str.find( ' ' ) + 1 );

        write << read_str << endl;
        read.close();
    }
    write.close();
    return 0;
}
