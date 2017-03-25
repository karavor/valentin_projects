#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <cmath>

using namespace std;
const char* input_file    = "./input/CAR_input";
const char* output_file   = "./output/CAR_result";
const char* conc_input    = "./input/complex.in"; // file complex.con should exist
const char* conc_output   = "./input/complex.eq";

#define CONCENTRATIONS_FUNC "\
    #/bin/bash \n\
    export NUPACKHOME=/home/valentin/work_backup/nupack \n\
    $NUPACKHOME/bin/complexes -T 23 -material rna -quiet $PWD/input/complex \n\
    $NUPACKHOME/bin/concentrations -quiet $PWD/input/complex \n\
"

#define CONCENTRATIONS_FUNC_DNA "\
    #/bin/bash \n\
    export NUPACKHOME=/home/valentin/work_backup/nupack \n\
    $NUPACKHOME/bin/complexes -T 23 -material dna -quiet $PWD/input/complex \n\
    $NUPACKHOME/bin/concentrations -quiet $PWD/input/complex \n\
"

int decimal_numb_reader( string* s, int str_position, float* numb, int minus_flag = 0 )
{
    (*numb) = 0;
    int i = str_position;

    if ( (*s)[i] == '-' )
        return decimal_numb_reader ( s, str_position + 1, numb, 1 );

    int s_length = (*s).length();
    int main_numb_part[10];

    while ( (*s)[i] != '\t' &&
            (*s)[i] != '.'  &&
            i < s_length       )
    {
        main_numb_part[i - str_position] = (*s)[i] - '0';
        i++;
    }
    i--;
    int j = 0;
    while ( i - str_position + 1 )
    {
        (*numb) += main_numb_part[j] * pow( 10, i - str_position );
        j++;
        i--;
    }
    i = str_position + j + 1;
    if ( (*s)[i - 1] == '.' )
    {
        j = -1;
        while ( (*s)[i] != '\t' &&
                (*s)[i] != '\n' &&
                i < s_length       )
        {
            (*numb) += int( (*s)[i] - '0' ) * pow( 10, j );
            j--;
            i++;
        }
        i++;
    }
    if ( minus_flag )
        (*numb) *= (-1);

    return i;
}

float exp_numb_reader(string* s)
{
    int i;
    string first_part_of_numb_str,
           second_part_of_numb_str;
    float first_part_of_numb,
          second_part_of_numb;

    for (i = 0; (*s)[i] != 'E' && (*s)[i] != 'e'; i++)
        first_part_of_numb_str.push_back((*s)[i]);

    i++;
    int minus_flag;
    if ((*s)[i] == '+') minus_flag = 1;
    else minus_flag = -1;

    for (i = i + 1; i < (*s).length(); i++)
        second_part_of_numb_str.push_back((*s)[i]);

    decimal_numb_reader(&first_part_of_numb_str, 0, &first_part_of_numb);
    decimal_numb_reader(&second_part_of_numb_str, 0, &second_part_of_numb);

    return first_part_of_numb * pow(10, (int)(minus_flag * second_part_of_numb) );
}


float aff_reader (string* seq_1, string* seq_2, int dna_flag)
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

    for(int i = 20; i < temp_str.length(); i++)
        aff_str.push_back(temp_str[i]);

    float aff = exp_numb_reader(&aff_str);
    return (aff/1e-06)*100;
}

int main ()
{
    char sep = '\t';

    int oligs_nmb;
    int dna_flag;

    float read_fl;

    string read_str;

    ifstream inp_r;
    inp_r.open(input_file);

    getline(inp_r, read_str);
    decimal_numb_reader(&read_str, (int)read_str.find(sep) + 1, &read_fl);
    oligs_nmb = (int)read_fl;

    getline(inp_r, read_str);
    decimal_numb_reader(&read_str, (int)read_str.find(sep) + 1, &read_fl);
    dna_flag = (int)read_fl;

    string oligs [oligs_nmb];

    for (int i = 0; i < oligs_nmb; i++)
    {
        getline(inp_r, read_str);
        read_str.erase(read_str.begin(), read_str.begin() + read_str.find(sep) + 1);
        oligs[i] = read_str;
    }
    inp_r.close();

    float cross_aff [oligs_nmb][oligs_nmb];

    ofstream out_w;
    out_w.open(output_file);

    for (int i = 1; i <= oligs_nmb; i++)

        out_w << sep << i;

    for (int i = 0; i < oligs_nmb; i++)
    {
        out_w << endl;
        out_w << (i + 1);
        for (int j = i; j < oligs_nmb; j++)
        {
            cross_aff [i][j] = aff_reader( &(oligs[i]), &(oligs[j]), dna_flag);
            cross_aff [j][i] = cross_aff [i][j];
        }
        for (int j = 0; j < oligs_nmb; j++)
        {
            out_w << sep << cross_aff [i][j];
        }
    }
    out_w.close();
    return 0;
}
