#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <cmath>
#include <vector>
#include <iomanip>

using namespace std;

const char* conc_input    = "./input/complex.in"; // file complex.con should exist
const char* conc_output   = "./input/complex.eq";
const char* input_mfe    = "./input/input_seq.in";
const char* output_mfe    = "./input/input_seq.mfe";
const char* test_input    = "./test_input";
const char* test_output    = "./test_output";

#define MFE_FUNC_RNA "\
    #/bin/bash \n\
    export NUPACKHOME=/home/valentin/work_backup/nupack \n\
    export CODDING_HOME=/home/valentin/work_backup/codding \n\
    $NUPACKHOME/bin/mfe -material rna $CODDING_HOME/input/input_seq \n\
"

#define MFE_FUNC_DNA "\
    #/bin/bash \n\
    export NUPACKHOME=/home/valentin/work_backup/nupack \n\
    export CODDING_HOME=/home/valentin/work_backup/codding \n\
    $NUPACKHOME/bin/mfe -material dna $CODDING_HOME/input/input_seq \n\
"

#define DESIGN_FUNC "\
    #/bin/bash \n\
    export NUPACKHOME=/home/valentin/work_backup/nupack \n\
    export CODDING_HOME=/home/valentin/work_backup/codding \n\
    $NUPACKHOME/bin/design -loadseed $CODDING_HOME/input/input_seq \n\
"

#define CONCENTRATIONS_FUNC_RNA "\
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

struct A_n_k_elmnt_with_nmb
{
    unsigned short int nmb;
    vector<unsigned short int> elmnt;
};

bool A_n_k_elmnt_test ( A_n_k_elmnt_with_nmb* A_n_k_elmnt,
                        vector< vector<unsigned short int> >* forbidden_A_n_k_elmnts_nmbs_stat)
{
    if ( (*A_n_k_elmnt).elmnt[ (*A_n_k_elmnt).elmnt.size() ] == 2) //for example
    {
        (*forbidden_A_n_k_elmnts_nmbs_stat)
        [(*forbidden_A_n_k_elmnts_nmbs_stat).size() - 1].
        push_back( (*A_n_k_elmnt).nmb );

        return 0;
    }

    return 1;
}

bool forbidden_A_n_k_elmnt_nmb_test( unsigned short int n,
                                     unsigned short int current_A_n_k_elmnt_nmb,
                                     vector< vector<unsigned short int> >* forbidden_A_n_k_elmnts_nmbs_stat )
{
    for (int i = 0; i < (*forbidden_A_n_k_elmnts_nmbs_stat).size(); i++)
    {
        for (int j = 0; j < (*forbidden_A_n_k_elmnts_nmbs_stat)[i].size(); j++)
        {
            if ( ( current_A_n_k_elmnt_nmb -
                   (*forbidden_A_n_k_elmnts_nmbs_stat)[i][j] ) % int_pow(n, i + 1) == 0 )
            {
                return 0;
            }
        }
    }
    return 1;
}

bool A_n_k_maker ( unsigned short int n,
                   unsigned short int k,
                   vector< A_n_k_elmnt_with_nmb >* A_n_k )
{
    if ( n < k )
        return 1;

    A_n_k_elmnt_with_nmb current_A_n_k_elmnt;
    vector< vector<unsigned short int> > forbidden_A_n_k_elmnts_nmbs_stat;
    vector<unsigned short int> init_vect;

    forbidden_A_n_k_elmnts_nmbs_stat.push_back(init_vect);

    for (unsigned short int i = 0; i < n; i++) // arr initialisation
    {
        current_A_n_k_elmnt.elmnt.push_back(i);
        current_A_n_k_elmnt.nmb = i;

        if ( A_n_k_elmnt_test( &current_A_n_k_elmnt,
                               &forbidden_A_n_k_elmnts_nmbs_stat ) )
        {
            (*A_n_k).push_back(current_A_n_k_elmnt);
        }

        current_A_n_k_elmnt.elmnt.clear();
    }

    int fix_arr_size;
    int insert_nmb;
    unsigned short int last_erased_elmnt_nmb;

    for (int itertn_nmb = 1; itertn_nmb < k; itertn_nmb++)
    {
        insert_nmb = 0;
        current_A_n_k_elmnt.nmb = 0;
        last_erased_elmnt_nmb = 0;

        fix_arr_size = (*A_n_k).size();

        forbidden_A_n_k_elmnts_nmbs_stat.push_back(init_vect);

        for (int i = 0; i < fix_arr_size; i++)
        {
            current_A_n_k_elmnt.elmnt = ((*A_n_k)[insert_nmb]).elmnt;

            current_A_n_k_elmnt.nmb = n * ((*A_n_k)[insert_nmb]).nmb;

            if ( ((*A_n_k)[insert_nmb]).nmb - 1 > last_erased_elmnt_nmb )
            {
                for (int j = (last_erased_elmnt_nmb + 1) * n;
                         j < ((*A_n_k)[insert_nmb]).nmb  * n; j++)
                {
                    forbidden_A_n_k_elmnts_nmbs_stat
                    [forbidden_A_n_k_elmnts_nmbs_stat.size() - 1].
                    push_back( j );
                }
            }

            (*A_n_k).erase( (*A_n_k).begin() + insert_nmb );

            last_erased_elmnt_nmb = ((*A_n_k)[insert_nmb]).nmb;

            for (int j = 0; j < n; j++)
            {
                if ( forbidden_A_n_k_elmnt_nmb_test( n,
                                                     current_A_n_k_elmnt.nmb,
                                                     &forbidden_A_n_k_elmnts_nmbs_stat ) )
                {
                    if ( A_n_k_elmnt_test( &current_A_n_k_elmnt,
                                           &forbidden_A_n_k_elmnts_nmbs_stat ) )
                    {
                        current_A_n_k_elmnt.elmnt.push_back(j);

                        (*A_n_k).insert( (*A_n_k).begin() + insert_nmb,
                                  current_A_n_k_elmnt);
                        insert_nmb ++;

                        current_A_n_k_elmnt.elmnt.
                        erase( current_A_n_k_elmnt.elmnt.end() - 1 );
                    }
                    else forbidden_A_n_k_elmnts_nmbs_stat
                         [forbidden_A_n_k_elmnts_nmbs_stat.size() - 1].
                         push_back( current_A_n_k_elmnt.nmb );
                }

                current_A_n_k_elmnt.nmb ++;
            }
        }
    }
    return 0;
}

int main()
{   // DESIGN_FUNC
    //MFE_FUNC_RNA MFE_FUNC_DNA
    //CONCENTRATIONS_FUNC_RNA CONCENTRATIONS_FUNC_DNA
    //system (MFE_FUNC_DNA);

    return 0;
}
