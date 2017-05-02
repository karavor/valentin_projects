#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <cmath>
#include <vector>

using namespace std;

const char* input_file    = "./input/OS_input";
const char* output_file    = "./output/OS_result";
const char* conc_input    = "./input/complex.in"; // file complex.con should exist
const char* conc_output   = "./input/complex.eq";
const char* input_mfe    = "./input/input_seq.in";
const char* output_mfe    = "./input/input_seq.mfe";

#define MFE_FUNC "\
    #/bin/bash \n\
    export NUPACKHOME=/home/valentin/work_backup/nupack \n\
    $NUPACKHOME/bin/mfe $PWD/input/input_seq \n\
"

#define CONCENTRATIONS_FUNC "\
    #/bin/bash \n\
    export NUPACKHOME=/home/valentin/work_backup/nupack \n\
    $NUPACKHOME/bin/complexes -T 23 -material rna -quiet $PWD/input/complex \n\
    $NUPACKHOME/bin/concentrations -quiet $PWD/input/complex \n\
"

#define CONCENTRATIONS_FUNC_DNA "\
    #/bin/bash \n\
    export NUPACKHOME=/home/valentin/work_backup/nupack \n\
    $NUPACKHOME/bin/complexes -T 23 -material dna $PWD/input/complex \n\
    $NUPACKHOME/bin/concentrations $PWD/input/complex \n\
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

float aff_reader (string* seq_1, string* seq_2, int dna_flag = 0)
{
    ofstream w;
    w.open(conc_input);
    w << 2        << endl
      << (*seq_1) << endl
      << (*seq_2) << endl
      << 2;
    w.close();
    if (dna_flag)
        system(CONCENTRATIONS_FUNC_DNA);
    else system(CONCENTRATIONS_FUNC);

    string temp_str;
    ifstream r;
    r.open(conc_output);

    do getline(r, temp_str);
    while (temp_str[0] == '%');

    while (temp_str[2] != temp_str[4])
        getline(r, temp_str);

    r.close();
    string aff_str;

    for(int i = 20; temp_str[i] != '\t'; i++)
        aff_str.push_back(temp_str[i]);

    return 1.0e8 * (float)exp_numb_reader(&aff_str);
}

float dG_and_struct_reader (string* seq, string* seq_struct)
{
    ofstream w;
    w.open(input_mfe);
    w << *seq;
    w.close();

    system(MFE_FUNC);

    string temporal_str;

    ifstream r;
    r.open(output_mfe);

    do    {getline(r, temporal_str);}
    while (temporal_str[0] != '-' &&
           temporal_str[0] != '0'    );

    getline( r, (*seq_struct) );

    r.close();

    return decimal_numb_reader( &temporal_str );
}

int left_decade_pop (int* left_decade_nmb,
                     int* right_decade_nmb,
                     int  max_nmb)
{
    if ( (*right_decade_nmb) > max_nmb )
    {
        (*right_decade_nmb) = 0;
        (*left_decade_nmb)++;
        return 0;
    }
    return 1;
}

void olig_from_int_arr_maker( int* int_arr,
                              int int_arr_length,
                              string* olig,
                              bool dna_flag      )
{
    (*olig).clear();

    for (int i = 0; i < int_arr_length; i++)
    {
        switch( int_arr[i] )
        {
            case 0:
                (*olig).push_back('A');
                break;

            case 1:
                if (dna_flag)
                     (*olig).push_back('T');
                else (*olig).push_back('U');
                break;

            case 2:
                (*olig).push_back('G');
                break;

            case 3:
                (*olig).push_back('C');
                break;
        }
    }
}

bool work_with_current_olig (float* best_dG,
                             string* best_olig,
                             bool dna_flag,
                             float   crit_dG,
                             ofstream* write,
                             string* current_olig,
                             string* left_olig_part,
                             string* right_olig_part,
                             vector<string>* antiaff_oligs)
{
    string bracket = (*left_olig_part) +
                     (*current_olig) +
                     (*right_olig_part);

    string bracket_struct;

    float br_dG = dG_and_struct_reader(&bracket, &bracket_struct);

    if ( br_dG > (*best_dG) )
    {
        (*best_dG) = br_dG;
        (*best_olig) = bracket;
    }

    if (br_dG > crit_dG)
    {
        (*write) << endl << endl << bracket << ' ' << br_dG << endl
                 << bracket_struct << endl;

        (*write) << "0) " << aff_reader(&bracket, &bracket, dna_flag) << endl
                 << "1) " << aff_reader(&bracket, &( (*antiaff_oligs)[1] ), dna_flag) << endl
                 << "4) " << aff_reader(&bracket, &( (*antiaff_oligs)[4] ), dna_flag) << endl
                 << "8) " << aff_reader(&bracket, &( (*antiaff_oligs)[8] ), dna_flag) << endl
                 << "9) " << aff_reader(&bracket, &( (*antiaff_oligs)[9] ), dna_flag) << endl;

        bracket = (*left_olig_part) + (*current_olig);
        (*write) << "   aff_of_B2O2:" << endl
                 << "2) " << aff_reader(&bracket, &( (*antiaff_oligs)[2] ), dna_flag) << endl
                 << "7) " << aff_reader(&bracket, &( (*antiaff_oligs)[7] ), dna_flag) << endl;

        bracket = (*current_olig) + (*right_olig_part);
        (*write) << "   aff_of_O2A2:" << endl
                 << "3) " << aff_reader(&bracket, &( (*antiaff_oligs)[3] ), dna_flag) << endl
                 << "6) " << aff_reader(&bracket, &( (*antiaff_oligs)[6] ), dna_flag);
    }
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

int main()
{
    int n = 4,
        k = 12;
    bool dna_flag = 1;
    int comb_nmb = (int)pow(n, k);
    int n_pow_k_arr [k];
    string current_olig;

    for (int i = 0; i < k; i++)
        n_pow_k_arr[i] = 0;

    string left_olig_part,
           right_olig_part;

    vector<string>antiaff_oligs(1);
    antiaff_oligs[0] = "A";

    float crit_dG,
          best_dG = -1000;
    string best_olig;

    ofstream write;
    write.open(output_file);

    input_data_reader(&crit_dG,
                      &left_olig_part,
                      &right_olig_part,
                      &antiaff_oligs   );

    for (int i = 0; i < comb_nmb; i++)
    {
        olig_from_int_arr_maker(n_pow_k_arr,
                                k,
                                &current_olig,
                                dna_flag      );

        work_with_current_olig( &best_dG,
                                &best_olig,
                                dna_flag,
                                crit_dG,
                                &write,
                                &current_olig,
                                &left_olig_part,
                                &right_olig_part,
                                &antiaff_oligs   );

        n_pow_k_arr[0]++;

        for (int l = 0; l < (k - 1); l++)
        {
            if ( left_decade_pop( &(n_pow_k_arr[l + 1]),
                                  &(n_pow_k_arr[l]),
                                  (n - 1)               ) )
                break;
        }
        /*if (i == 100)
            break;*/
    }
    write << endl << endl << "best_olig: " << best_olig << ' ' << best_dG;
    write.close();
    return 0;
}
