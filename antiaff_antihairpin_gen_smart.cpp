#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <cmath>
#include <vector>

using namespace std;
const char* input_file    = "./input/AAGS_input"; // text file input should exist
const char* output_file   = "./output/AAGS_result";
const char* statistic     = "./output/statistic";
const char* input_mfe     = "./input/input_seq.in";
const char* output_mfe    = "./input/input_seq.mfe";
const char* conc_input    = "./input/complex.in"; // file complex.con should exist
const char* conc_output   = "./input/complex.eq";
const char* rndm_seed_file   = "./input/rndm_seed";

int pause;
int rndm_seed = -1;
float max_aff_level = 1.0;

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

#define CLEAN_FUNC "\
    #/bin/bash \n\
    rm $PWD/input/input_seq.in \n\
    rm $PWD/input/input_seq.mfe \n\
    rm $PWD/input/complex.in \n\
    rm $PWD/input/complex.eq \n\
    rm $PWD/input/complex.cx \n\
"

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

    for(int i = 20; i < temp_str.length() - 1; i++)
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

char complement_ncltd (char ncltd, int dna_flag)
{
    switch(ncltd)
    {
        case 'C':
            return 'G';
        case 'G':
            return 'C';
        case 'A':
            if (dna_flag)
                return 'T';
            return 'U';
        case 'U':
            return 'A';
        case 'T':
            return 'A';
    }
    return 'E';
}

void string_inverser (string* s)
{
    string s_copy = (*s);
    int l = (*s).length();
    for (int i = 0; i < l; i++ )
    {
        (*s)[i] = s_copy[l-1-i];
    }
}

int c_seq_maker (string* seq, string* c_seq, int dna_flag)
{
    int l = (*seq).length();
    for (int i = 0; i < l; i++ )
    {
        switch((*seq)[i])
        {
            case 'C':
                (*c_seq).push_back('G');
                break;
            case 'G':
                (*c_seq).push_back('C');
                break;
            case 'A':
                if (dna_flag)
                     (*c_seq).push_back('T');
                else (*c_seq).push_back('U');
                break;
            case 'U':
                (*c_seq).push_back('A');
                break;
            case 'T':
                (*c_seq).push_back('A');
                break;
        }
    }
    string_inverser(c_seq);
    return 0;
}

char random_ncltd (int dna_flag, int rnd_seed)
{
    srand (rnd_seed);
    int random_int = rand() % 6;

    if (random_int <= 1)
        return 'A';

    if (random_int <= 3)
    {
        if (dna_flag)
            return 'T';
        return 'U';
    }
    if (random_int == 4)
        return 'C';

    return 'G';
}

void olig_generator (int olig_length,
                     string* olig_seq,
                     int dna_flag,
                     int rndm_modifier = 0)
{
    for (int i = 0; i < olig_length; i++)

        (*olig_seq).push_back( random_ncltd(dna_flag, time(NULL) + i + rndm_modifier) );
}

void input_data_reader (int* dna_flag,
                        int* oligs_nmb,
                        int* olig_length,
                        string* initial_olig,
                        int* pattern_length,
                        vector<string>* oligs_group)
{
    string temp_str;

    ifstream read;
    read.open(input_file);

    getline(read, temp_str);
    temp_str.erase(temp_str.begin(),
                   temp_str.begin() + temp_str.find(' ') + 1);
    *dna_flag = int_numb_reader(&temp_str);

    getline(read, temp_str);
    if ( temp_str[temp_str.length() - 1] != '-')
    {
        temp_str.erase(temp_str.begin(),
                       temp_str.begin() + temp_str.find(' ') + 1);
        rndm_seed = int_numb_reader(&temp_str);
    }

    getline(read, temp_str);
    temp_str.erase(temp_str.begin(),
                   temp_str.begin() + temp_str.find(' ') + 1);
    *olig_length = int_numb_reader(&temp_str);

    getline(read, temp_str);
    if ( temp_str[temp_str.length() - 1] != '-')
    {
        temp_str.erase(temp_str.begin(),
                       temp_str.begin() + temp_str.find(' ') + 1);

        (*initial_olig) = temp_str;
    }

    getline(read, temp_str);
    temp_str.erase(temp_str.begin(),
                   temp_str.begin() + temp_str.find(' ') + 1);
    *pattern_length = int_numb_reader(&temp_str);

    getline(read, temp_str);
    temp_str.erase(temp_str.begin(),
                   temp_str.begin() + temp_str.find(' ') + 1);
    *oligs_nmb = int_numb_reader(&temp_str);

    (*oligs_group).resize(*oligs_nmb);
    for (int i = 0; i < (*oligs_nmb); i++)
    {
        getline( read, (*oligs_group)[i] );
        (*oligs_group)[i].erase((*oligs_group)[i].begin(),
                                (*oligs_group)[i].begin() + (*oligs_group)[i].find(' ') + 1);
    }
    read.close();
}

bool string_in_base_finder (string*         str,
                            vector<string>* string_base)
{
    for (int i = 0; i < (int)(*string_base).size(); i++)
    {
        if ( (*str) == (*string_base)[i] )
            return 1;
    }
    return 0;
}

void prevent_pattern_base_maker (vector<string>* oligs_group,
                                 int     oligs_nmb,
                                 int     pattern_length,
                                 vector<string>* prevent_pattern_base)
{
    int pattern_oligs_nmb;
    string new_pattern;

    for (int i = 0; i < oligs_nmb; i++)
    {
        pattern_oligs_nmb = ((*oligs_group)[i]).length() - pattern_length + 1;

        for (int j = 0; j < pattern_oligs_nmb; j++)
        {
            for (int k = 0; k < pattern_length; k++)
            {
                new_pattern.push_back((*oligs_group)[i][j+k]);
            }
            if ( string_in_base_finder(&new_pattern, prevent_pattern_base) == 0)
            {
                (*prevent_pattern_base).push_back(new_pattern);
            }
            new_pattern.clear();
        }
    }
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

char ncltd_generator (int ncltd_nmb, int dna_flag) //except_ncltd_nmb can be 0, 1, 2 or 3
{
    char ncltds [4] = {'A', 'U', 'G', 'C'};

    if (dna_flag)
        ncltds[1] = 'T';

    return ncltds[ncltd_nmb];
}

template <class data_type>

void vector_element_deleter (data_type* element,
                             vector<data_type>* vect)
{
    for (int i = 0; i < (int)(*vect).size(); i++)
    {
        if ( (*element) == (*vect)[i] )
        {
            (*vect).erase( (*vect).begin() + i );
            break;
        }
    }
}

void complement_pair_deleter (int dna_flag,
                              string* pattern,
                              vector<string>* allowed_pattern_base,
                              vector<string>* prevent_pattern_base )
{
    string c_pattern;
    c_seq_maker(pattern, &c_pattern, dna_flag);

    vector_element_deleter<string> (&c_pattern, allowed_pattern_base);
    (*prevent_pattern_base).push_back( c_pattern );
}

bool self_aff_pattern_test (string* pattern,
                            int dna_flag    )
{
    string c_pattern;
    c_seq_maker(pattern, &c_pattern, dna_flag);

    if ( (*pattern) == c_pattern )
    {
        return 1;
    }
    return 0;
}

void allowed_pattern_base_maker (int pattern_length,
                                 int dna_flag,
                                 vector<string>* prevent_pattern_base,
                                 vector<string>* allowed_pattern_base)
{
    string current_pattern;
    int n = 4,
        k = pattern_length;

    int comb_nmb = (int)pow(n, k);
    int n_pow_k_arr [k];

    for (int i = 0; i < k; i++)
        n_pow_k_arr[i] = 0;

    for (int i = 0; i < comb_nmb; i++)
    {
        for (int j = 0; j < k; j++)
        {
            current_pattern.push_back( ncltd_generator (n_pow_k_arr[j], dna_flag) );
        }
        if ( string_in_base_finder (&current_pattern, prevent_pattern_base) == 0 )
        {
            if ( self_aff_pattern_test (&current_pattern, dna_flag) )

                (*prevent_pattern_base).push_back(current_pattern);
            else
                (*allowed_pattern_base).push_back(current_pattern);
        }
        current_pattern.clear();
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

bool rndm_int_without_repeat_generator (vector<int>* int_vect,
                                        int int_nmb,
                                        int min_int,
                                        int max_int)
{
    int elements_nmb = max_int - min_int + 1;

    if (int_nmb > elements_nmb)
    {
        cout << endl << "Error in function rndm_int_without_repeat_generator";
        return 1;
    }

    vector<int> new_int_vect(elements_nmb);

    for (int i = 0; i < elements_nmb; i++)
        new_int_vect[i] = i;

    int rndm_int;
    string read_str;
    float read_fl;
    ifstream read;

    if (rndm_seed == -1)
    {
        rndm_seed = (int)time(NULL);
    }

    for (int i = 0; i < int_nmb; i++)
    {
        srand ( rndm_seed + i*17 );

        rndm_int = rand() % (elements_nmb - i);
        (*int_vect).push_back( min_int + new_int_vect[rndm_int] );

        new_int_vect[rndm_int] =
        new_int_vect[ (int)new_int_vect.size() - 1 - i ];
    }
    return 0;
}

bool initial_olig_maker (int dna_flag,
                         int olig_length,
                         int pattern_length,
                         string* initial_olig,
                         vector<string>* allowed_pattern_base,
                         vector<string>* prevent_pattern_base )
{
    int patterns_nmb = 1 + (int)(olig_length/pattern_length);
    vector<int> rndm_int_vect (0);
    string building_olig;

    if (rndm_int_without_repeat_generator ( &rndm_int_vect,
                                            patterns_nmb,
                                            0,
                                            -1 - 2*patterns_nmb + (int)(*allowed_pattern_base).size() ) )
    {
        cout << endl << "There are too few patterns in allowed_pattern_base" << endl;
        return 1;
    }

    for (int i = 0; i < patterns_nmb; i++)
    {
        building_olig = (*allowed_pattern_base)[ rndm_int_vect[i] ];
        (*initial_olig) += building_olig;

        (*allowed_pattern_base).erase( (*allowed_pattern_base).begin() + rndm_int_vect[i] );

        complement_pair_deleter (dna_flag,
                                 &building_olig,
                                 allowed_pattern_base,
                                 prevent_pattern_base );

        //(*prevent_pattern_base).push_back( building_olig ); //This may lead to definite areas become more specific
    }
    if ((*initial_olig).length() > olig_length)
        (*initial_olig).erase( (*initial_olig).begin() + olig_length, (*initial_olig).end() );

    return 0;
}

void current_pattern_base_maker (int olig_length,
                                 int pattern_length,
                                 string* initial_olig,
                                 vector<string>* current_pattern_base)
{
    string current_pattern;

    for (int i = 0; i < (olig_length - pattern_length + 1); i++)
    {
        for (int j = 0; j < pattern_length; j++)
        {
            current_pattern.push_back( (*initial_olig)[i+j] );
        }
        (*current_pattern_base).push_back(current_pattern);
        current_pattern.clear();
    }
}

void current_pattern_base_corrector (int dna_flag,
                                     vector<string>* prevent_pattern_base,
                                     vector<string>* allowed_pattern_base,
                                     vector<string>* current_pattern_base )
{
    int pattern_length = ( (*prevent_pattern_base)[0] ).length();
    int prevent_neighbour_patterns_nmb;
    int min_prevent_neighbour_patterns_nmb;
    int new_patterns_group_size = 2*( pattern_length - 1 ) + 1;
    string best_new_patterns [new_patterns_group_size];
    bool break_flag;

    for (int i = 0; i < (int)(*current_pattern_base).size(); i++)
    {
        if ( string_in_base_finder(&((*current_pattern_base)[i]),
                                   prevent_pattern_base      ) )
        {
            break_flag = 0;
            min_prevent_neighbour_patterns_nmb = 2 * (int)(*current_pattern_base).size();

            for (int j = 0; j < (int)(*allowed_pattern_base).size(); j++)
            {
                (*current_pattern_base)[i] = (*allowed_pattern_base)[j];
                prevent_neighbour_patterns_nmb = 0;

                for (int k = pattern_length - 1; k > 0; k--) //in this cycle we change neighbour dependent patterns
                {
                    if (i - k >= 0)
                    {
                        for (int l = 0; l < pattern_length - k; l++)
                        {
                            (*current_pattern_base)[i - k][k + l] = (*allowed_pattern_base)[j][l];
                        }
                        if ( string_in_base_finder(&((*current_pattern_base)[i - k]),
                                                   prevent_pattern_base)             )
                            prevent_neighbour_patterns_nmb++;
                    }
                    if (i + k < (int)(*current_pattern_base).size())
                    {
                        for (int l = 0; l < pattern_length - k; l++)
                        {
                            (*current_pattern_base)[i + k][pattern_length - 1 - k - l] =
                            (*allowed_pattern_base)[j][pattern_length - 1 - l];
                        }
                        if ( string_in_base_finder(&((*current_pattern_base)[i + k]),
                                                   prevent_pattern_base)             )
                            prevent_neighbour_patterns_nmb++;
                    }
                }
                if (prevent_neighbour_patterns_nmb < min_prevent_neighbour_patterns_nmb)
                {
                    min_prevent_neighbour_patterns_nmb = prevent_neighbour_patterns_nmb;

                    if (min_prevent_neighbour_patterns_nmb == 0)
                    {
                        break_flag = 1;
                        break;
                    }

                    for (int k = 0; k < new_patterns_group_size; k++)
                    {
                        if ( i - (pattern_length - 1) + k < 0 ||
                             i - (pattern_length - 1) + k >= (int)(*current_pattern_base).size() )
                            continue;
                        best_new_patterns[k] = (*current_pattern_base)[i - (pattern_length - 1) + k];
                    }
                }
            }
            if (break_flag == 0)
            {
                for (int k = 0; k < new_patterns_group_size; k++)
                {
                    if ( i - (pattern_length - 1) + k < 0 ||
                         i - (pattern_length - 1) + k >= (int)(*current_pattern_base).size() )
                        continue;
                    (*current_pattern_base)[i - (pattern_length - 1) + k] = best_new_patterns[k];
                }
            }
            complement_pair_deleter ( dna_flag, &( (*current_pattern_base)[i] ), allowed_pattern_base, prevent_pattern_base );
            vector_element_deleter<string> (    &( (*current_pattern_base)[i] ), allowed_pattern_base );
            (*prevent_pattern_base).push_back( (*current_pattern_base)[i] ); //This may lead to definite areas become more specific
        }
    }
}

void current_patterns_integrator (string* result,
                                  vector<string>* current_pattern_base)
{
    int i;
    for (i = 0; i < (int)(*current_pattern_base).size() - 1; i++)
    {
        (*result).push_back( (*current_pattern_base)[i][0] );
    }
    (*result) += (*current_pattern_base)[i];
}

int main()
{
    int dna_flag,
        oligs_nmb,
        olig_length,
        pattern_length;

    vector<string> oligs_group(0);

    string initial_olig;

    input_data_reader (&dna_flag,
                       &oligs_nmb,
                       &olig_length,
                       &initial_olig,
                       &pattern_length,
                       &oligs_group);

    vector<string> prevent_pattern_base (0);

    prevent_pattern_base_maker (&oligs_group,
                                oligs_nmb,
                                pattern_length,
                                &prevent_pattern_base);

    vector<string> allowed_pattern_base (0);

    allowed_pattern_base_maker (pattern_length,
                                dna_flag,
                                &prevent_pattern_base,
                                &allowed_pattern_base);

    ofstream error_flag;

    if (initial_olig.length() == 0)
    {
        if ( initial_olig_maker (dna_flag,
                                 olig_length,
                                 pattern_length,
                                 &initial_olig,
                                 &allowed_pattern_base,
                                 &prevent_pattern_base ) )
        {

            error_flag.open("./error_flag");
            error_flag << 1;
            error_flag.close();
            return 0;
        }
    }

    vector<string> current_pattern_base (0);
    current_pattern_base_maker (olig_length,
                                pattern_length,
                                &initial_olig,
                                &current_pattern_base);

    current_pattern_base_corrector (dna_flag,
                                    &prevent_pattern_base,
                                    &allowed_pattern_base,
                                    &current_pattern_base );

    ofstream write;
    write.open(statistic);

    write << "prevent_pattern_base: " << endl;
    for (int i = 0; i < (int)prevent_pattern_base.size(); i++)
        write << prevent_pattern_base[i] << endl;

    write << endl << "current_pattern_base: " << endl;
    write << initial_olig << endl;
    for (int i = 0; i < (int)current_pattern_base.size(); i++)
        write << endl << current_pattern_base[i];

    write << endl << endl << "allowed_pattern_base: ";
    for (int i = 0; i < (int)allowed_pattern_base.size(); i++)
        write << endl << allowed_pattern_base[i];

    write.close();

    string result;
    current_patterns_integrator(&result, &current_pattern_base);

    string c_result;
    c_seq_maker(&result, &c_result, dna_flag);

    string result_struct;
    float dG = dG_and_struct_reader(&result, &result_struct);

    string c_result_struct;
    float c_dG = dG_and_struct_reader(&c_result, &c_result_struct);

    write.open(output_file);
    write << c_dG + dG << endl;
    write << result << endl
          << result_struct << endl
          << "dG = " << dG << " kkal/mol" << endl << endl;

    write << "Aff to oligs, %: " << endl;

    vector<float> aff_vect (0);

    aff_vect.push_back( aff_reader(&result, &result, dna_flag) );
    write << "0) " << aff_vect[0] << endl;

    for (int i = 1; i <= oligs_nmb; i++)
    {
        aff_vect.push_back( aff_reader(&result, &(oligs_group[i - 1]), dna_flag) );
        write << i << ") " << aff_vect [i] << endl;
    }

    write << endl << "for complement olig: " << endl
                  << c_result << endl
                  << c_result_struct << endl
                  << "dG = " << c_dG << " kkal/mol" << endl << endl;

    write << "Aff to oligs, %: " << endl;

    aff_vect.push_back( aff_reader(&c_result, &c_result, dna_flag) );
    write << "0) " << aff_vect[oligs_nmb + 1] << endl;

    for (int i = 1; i <= oligs_nmb; i++)
    {
        aff_vect.push_back( aff_reader(&c_result, &(oligs_group[i - 1]), dna_flag) );
        write << i << ") " << aff_vect[oligs_nmb + 1 + i] << endl;
    }

    write << endl << "rndm_seed: " << rndm_seed << endl;

    int break_flag = 0;

    for (int i = 0; i < (int)aff_vect.size(); i++)
    {
        if ( aff_vect[i] > max_aff_level )
        {
            write << "!-";
            break_flag = 1;
            break;
        }
    }
    if (break_flag == 0)
        write << "!+";

    write.close();
    system(CLEAN_FUNC);

    error_flag.open("./error_flag");
    error_flag << 0;
    error_flag.close();

    return 0;
}
