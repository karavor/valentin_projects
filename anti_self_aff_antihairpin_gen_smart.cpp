#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <cmath>
#include <vector>
#include <iomanip>

using namespace std;
const char* AAGS_input_file    = "./input/AAGS_input"; // text file input should exist
const char* AAGS_result_file   = "./output/AAGS_result";
const char* input_mfe     = "./input/input_seq.in";
const char* output_mfe    = "./input/input_seq.mfe";

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

float exp_numb_reader(string* s)
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

#define MFE_FUNC_DNA "\
    #/bin/bash \n\
    export NUPACKHOME=/home/valentin/work_backup/nupack \n\
    $NUPACKHOME/bin/mfe -material dna $PWD/input/input_seq \n\
"

#define MFE_FUNC_RNA "\
    #/bin/bash \n\
    export NUPACKHOME=/home/valentin/work_backup/nupack \n\
    $NUPACKHOME/bin/mfe -material rna $PWD/input/input_seq \n\
"

#define CLEAN_FUNC "\
    #/bin/bash \n\
    rm $PWD/input/input_seq.in \n\
    rm $PWD/input/input_seq.mfe \n\
"

float dG_reader (string* seq, bool dna_flag)
{
    ofstream w;
    w.open(input_mfe);
    w << *seq;
    w.close();

    if (dna_flag)
         system(MFE_FUNC_DNA);
    else system(MFE_FUNC_RNA);

    string temporal_str;

    ifstream r;
    r.open(output_mfe);

    for (int i = 0; i < 15; i++)
    {
        getline(r, temporal_str);
    }

    if (temporal_str.length() == 0)
        return 0;

    return decimal_numb_reader( &temporal_str);
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

char ncltd_generator (int ncltd_nmb, bool dna_flag) //ncltd_nmb can be 0, 1, 2 or 3
{
    char ncltds [4] = {'A', 'U', 'G', 'C'};

    if (dna_flag)
        ncltds[1] = 'T';

    return ncltds[ncltd_nmb];
}

struct pttrn_with_self_dG_score
{
    string pttrn;
    float self_dG_score;
};

float double_pttrn_score_calc (string* pttrn_1,
                               string* pttrn_2,
                               bool dna_flag,
                               float max_pttrn_cmplx_dG)
{
    string pttrn_self_cmplx = (*pttrn_1) + "+" + (*pttrn_2);
    return abs( dG_reader(&pttrn_self_cmplx, dna_flag)/max_pttrn_cmplx_dG );
}

void pttrn_base_generator (vector<pttrn_with_self_dG_score>* pttrn_base,
                           int pttrn_length,
                           bool dna_flag,
                           float max_pttrn_cmplx_dG,
                           float crit_score)
{
    pttrn_with_self_dG_score current_pttrn;
    int n = 4,
        k = pttrn_length;

    int comb_nmb = (int)pow(n, k);
    int n_pow_k_arr [k];

    for (int i = 0; i < k; i++)
        n_pow_k_arr[i] = 0;

    for (int i = 0; i < comb_nmb; i++)
    {
        for (int j = 0; j < k; j++)
        {
            (current_pttrn.pttrn).push_back( ncltd_generator (n_pow_k_arr[j], dna_flag) );
        }
        current_pttrn.self_dG_score = double_pttrn_score_calc (&(current_pttrn.pttrn),
                                                               &(current_pttrn.pttrn),
                                                               dna_flag,
                                                               max_pttrn_cmplx_dG);

        if ( current_pttrn.self_dG_score < crit_score )
        {
            (*pttrn_base).push_back(current_pttrn);
        }
        (current_pttrn.pttrn).clear();
        n_pow_k_arr[0]++;

        for (int l = 0; l < (k - 1); l++)
        {
            if ( left_decade_pop( &(n_pow_k_arr[l + 1]),
                                  &(n_pow_k_arr[l]),
                                  (n - 1)               ) )
                break;
        }
    }
}

float string_with_prmtr_reader(ifstream* read, char sep)
{
    string read_str;
    getline( (*read), read_str);
    read_str.erase( read_str.begin(), read_str.begin() + read_str.find( sep ) + 1 );
    return decimal_numb_reader(&read_str);
}

float max_pttrn_cmplx_dG_calc (int pttrn_length, bool dna_flag)
{
    string worst_pttrn_cmplx;

    for (int i = 0; i < pttrn_length; i++)
    {
        worst_pttrn_cmplx.push_back('C');
    }
    worst_pttrn_cmplx.push_back('+');

    for (int i = 0; i < pttrn_length; i++)
    {
        worst_pttrn_cmplx.push_back('G');
    }
    return dG_reader(&worst_pttrn_cmplx, dna_flag);
}

struct two_pttrn_gluing
{
    float cmplx_dG_score;
    vector<int> appeared_pttrns_nmbs; //in order from 5' to 3'
    float gluing_dG_score;
};

int  string_in_base_finder (string* pttrn,
                            vector<pttrn_with_self_dG_score>* pttrn_base)
{
    for (int i = 0; i < (int)(*pttrn_base).size(); i++)
    {
        if ( (*pttrn) == ((*pttrn_base)[i]).pttrn )
            return i;
    }
    //cout << endl << "There's no such pattern in base" << endl;
    return -1;
}

bool appeared_pttrns_of_gluing_nmbs (vector<pttrn_with_self_dG_score>* pttrn_base,
                                     two_pttrn_gluing*  gluing_matrix_element,
                                     int i,
                                     int j)
{
    ( (*gluing_matrix_element).appeared_pttrns_nmbs ).push_back( i );
    ( (*gluing_matrix_element).appeared_pttrns_nmbs ).push_back( j );

    string appeared_pttrn;
    string gluing = ( ((*pttrn_base)[i]).pttrn ) + ( ((*pttrn_base)[j]).pttrn );
    gluing.erase( gluing.begin() );
    gluing.erase( gluing.begin() + gluing.length() - 1 );

    int pttrn_length = (int)( ((*pttrn_base)[0]).pttrn ).length();
    int appeared_pttrn_nmb;

    for (int k = 0; k < pttrn_length - 1; k++)
    {
        for (int l = k; l < k + pttrn_length; l++)
        {
            appeared_pttrn.push_back( gluing[l] );
        }
        appeared_pttrn_nmb = string_in_base_finder (&appeared_pttrn, pttrn_base);

        if (appeared_pttrn_nmb < 0)
        {
            (*gluing_matrix_element).gluing_dG_score = 1;//this score we shouldn't change lower
            //return 0; // it will be more rationally, but now it's not convenient for analyzing
            ( (*gluing_matrix_element).appeared_pttrns_nmbs ).push_back( -1 );
            appeared_pttrn.clear();
            continue;
        }
        ( (*gluing_matrix_element).appeared_pttrns_nmbs ).push_back( appeared_pttrn_nmb );
        appeared_pttrn.clear();
    }
    return 0;
}

void gluing_matrix_maker (vector<pttrn_with_self_dG_score>* pttrn_base,
                          vector< vector<two_pttrn_gluing> >* gluing_matrix,
                          bool dna_flag,
                          float max_pttrn_cmplx_dG)
{
    two_pttrn_gluing gluing_matrix_element;
    vector<two_pttrn_gluing> gluing_matrix_string;

    for (int i = 0; i < (*pttrn_base).size(); i++) // string in gluing_matrix
    {
        for (int j = 0; j < (*pttrn_base).size(); j++) //column in gluing_matrix
        {
            gluing_matrix_element.gluing_dG_score = 0;

            if (i == j)
                gluing_matrix_element.cmplx_dG_score = ((*pttrn_base)[i]).self_dG_score;
            else
                gluing_matrix_element.cmplx_dG_score
                =
                double_pttrn_score_calc ( &( ((*pttrn_base)[i]).pttrn ),
                                          &( ((*pttrn_base)[j]).pttrn ),
                                          dna_flag,
                                          max_pttrn_cmplx_dG);

            appeared_pttrns_of_gluing_nmbs (pttrn_base,
                                            &gluing_matrix_element,
                                            i,
                                            j);

            gluing_matrix_string.push_back(gluing_matrix_element);
        }
        (*gluing_matrix).push_back(gluing_matrix_string);
        gluing_matrix_string.clear();
    }
}

template<class _T_>

bool C_n_k_generator (vector<_T_>* elements_arr,
                      vector< vector<_T_> >* result_arr,
                      int k)
{
    int n = (*elements_arr).size();

    if (n < k)
        return 1;

    vector<_T_> current_res;

    for (int j = 0; j <= n - k; j++)
    {
        for (int l = j; l < j + k; l++)
            current_res.push_back( (*elements_arr)[l] );

        (*result_arr).push_back(current_res);

        for (int i = j + k; i < n; i++)
        {
            current_res.erase(current_res.begin() + k - 1);
            current_res.push_back( (*elements_arr)[i] );
            (*result_arr).push_back(current_res);
        }
        current_res.clear();
    }
    return 0;
}

bool gluing_dG_score_calc (vector< vector<two_pttrn_gluing> >* gluing_matrix)
{
    vector< vector<int> > C_n_k_arr;
    float max_cmplx_dG_score;
    float current_cmplx_dG_score;

    for (int i = 0; i < (*gluing_matrix).size(); i++) // string in gluing_matrix
    {
        for (int j = 0; j < (*gluing_matrix).size(); j++) //column in gluing_matrix
        {
            if ( ((*gluing_matrix)[i][j]).gluing_dG_score == 0 )
            {
                if ( C_n_k_generator<int> (&( ((*gluing_matrix)[i][j]).appeared_pttrns_nmbs ),
                                           &C_n_k_arr,
                                           2                                                  ) )
                {
                    cout << endl << "C_n_k_generator fail: n < k" << endl;
                    return 1;
                }
                max_cmplx_dG_score = 0;

                for (int l = 0; l < C_n_k_arr.size(); l++)
                {
                    current_cmplx_dG_score = ((*gluing_matrix)[ C_n_k_arr[l][0] ]
                                                              [ C_n_k_arr[l][1] ]).cmplx_dG_score;

                    if ( current_cmplx_dG_score > max_cmplx_dG_score )
                    {
                        max_cmplx_dG_score = current_cmplx_dG_score;
                    }
                }
                ((*gluing_matrix)[i][j]).gluing_dG_score = max_cmplx_dG_score;
            }
        }
    }
    return 0;
}

void cmplx_dG_score_string_writer (char sep,
                                   ofstream* write,
                                   vector<two_pttrn_gluing>* gluing_matrix_str)
{
    (*write) << sep;

    for (int i = 0; i < (*gluing_matrix_str).size(); i++)
    {
        (*write) << sep << ( (*gluing_matrix_str)[i] ).cmplx_dG_score;
    }
    (*write) << endl;
}

void appeared_pttrns_string_writer (char sep,
                                    string* str_pttrn,
                                    ofstream* write,
                                    vector<two_pttrn_gluing>* gluing_matrix_str,
                                    vector<pttrn_with_self_dG_score>* pttrn_base)
{
    int current_pttrn_nmb;
    int pttrn_length = ( ( (*pttrn_base)[0] ).pttrn ).length();

    for (int k = 0; k < pttrn_length - 1; k++)
    {
        if (k == (int)( (pttrn_length - 1) / 2))
            (*write) << (*str_pttrn);
        else
            (*write) << sep;

        for (int i = 0; i < (*gluing_matrix_str).size(); i++)
        {
            (*write) << sep;
            current_pttrn_nmb = ( ((*gluing_matrix_str)[i]).appeared_pttrns_nmbs )[k + 2];

            if ( current_pttrn_nmb >= 0 )
                (*write) << ( (*pttrn_base)[ current_pttrn_nmb ] ).pttrn;
            else
                (*write) << '-' << sep;
        }
        (*write) << endl;
    }
}

void gluing_dG_score_string_writer (char sep,
                                    ofstream* write,
                                    vector<two_pttrn_gluing>* gluing_matrix_str)
{
    (*write) << sep;

    for (int i = 0; i < (*gluing_matrix_str).size(); i++)
    {
        (*write) << sep << ( (*gluing_matrix_str)[i] ).gluing_dG_score;
    }
    (*write) << endl << endl;
}

void gluing_matrix_writer( ofstream* write,
                           vector< vector<two_pttrn_gluing> >* gluing_matrix,
                           vector<pttrn_with_self_dG_score>* pttrn_base      )
{
    int n = (*gluing_matrix).size();
    char sep = '\t';

    (*write) << fixed << setprecision(2);
    (*write) << sep << (int)2;

    for (int i = 0; i < n; i++)
    {
        (*write) << sep << ( (*pttrn_base)[i] ).pttrn;
    }
    (*write) << endl << endl << (int)1;

    for (int i = 0; i < n; i++) // string in gluing_matrix
    {
        cmplx_dG_score_string_writer (sep,
                                      write,
                                      &((*gluing_matrix)[i]));

        appeared_pttrns_string_writer (sep,
                                       &( ((*pttrn_base)[i]).pttrn ),
                                       write,
                                       &((*gluing_matrix)[i]),
                                       pttrn_base                   );

        gluing_dG_score_string_writer (sep,
                                       write,
                                       &((*gluing_matrix)[i]));
    }
}

int main()
{
    ifstream read;
    read.open(AAGS_input_file);
    char sep = ' ';

    bool  dna_flag      = string_with_prmtr_reader(&read, sep);
    int   pttrn_length  = string_with_prmtr_reader(&read, sep),
          min_pttrn_nmb = string_with_prmtr_reader(&read, sep),
          max_pttrn_nmb = string_with_prmtr_reader(&read, sep);
    float crit_score    = string_with_prmtr_reader(&read, sep);

    read.close();

    float max_pttrn_cmplx_dG = max_pttrn_cmplx_dG_calc (pttrn_length, dna_flag);

    vector<pttrn_with_self_dG_score> pttrn_base;
    pttrn_base_generator (&pttrn_base,
                          pttrn_length,
                          dna_flag,
                          max_pttrn_cmplx_dG,
                          crit_score);

    vector< vector<two_pttrn_gluing> > gluing_matrix; //like column of the strings

    gluing_matrix_maker (&pttrn_base,
                         &gluing_matrix,
                         dna_flag,
                         max_pttrn_cmplx_dG);

    if ( gluing_dG_score_calc (&gluing_matrix) )
    {
        cout << endl << "Fail in function gluing_dG_score_calc" << endl;
        return 1;
    }

    ofstream write;

    write.open(AAGS_result_file);

    gluing_matrix_writer( &write,
                          &gluing_matrix,
                          &pttrn_base     );
    write.close();
/*
    ofstream write;

    write.open(AAGS_result_file);

    for (int i = 0; i < pttrn_base.size(); i++)
    {
        write << (pttrn_base[i]).pttrn << '\t' << (pttrn_base[i]).self_dG_score << endl;
    }

    write.close();

    //system(CLEAN_FUNC);
*/
    return 0;
}
