#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <cmath>

using namespace std;
const char* input_file    = "./input/input"; // ./input directory and text file input should exist
const char* input_mfe     = "./input/input_seq.in";
const char* output_mfe    = "./input/input_seq.mfe";
const char* input_design  = "./input/input_seq.fold";
const char* output_design = "./input/input_seq.summary";
const char* output_file   = "./output/result";
const char* statistic     = "./output/statistic";
const char* conc_input    = "./input/complex.in"; // file complex.con should exist
const char* conc_output   = "./input/complex.eq";

const int GEN_NEXT  = 0; //ok, print and continue
const int GEN_TERM  = 1 ;//ok, terminate
const int GEN_EMPTY = 2; //ok, print EMPTY SET and continue
const int GEN_ERROR = 3;

const int max_loop_repetition = 10;

#define MFE_FUNC "\
    #/bin/bash \n\
    export NUPACKHOME=/home/valentin/work_backup/nupack \n\
    $NUPACKHOME/bin/mfe $PWD/input/input_seq \n\
"
#define DESIGN_FUNC "\
    #/bin/bash \n\
    export NUPACKHOME=/home/valentin/work_backup/nupack \n\
    $NUPACKHOME/bin/design -loadseed $PWD/input/input_seq \n\
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
    rm $PWD/input/input_seq.fold \n\
    rm $PWD/input/input_seq.summary \n\
    rm $PWD/input/input_seq.mfe \n\
    rm $PWD/input/complex.in \n\
    rm $PWD/input/complex.eq \n\
    rm $PWD/input/complex.cx \n\
"
struct anti_aff_olig
{
    string seq;
    string cmplx_struct_part;
    float aff;
};

void string_inverser (string* s)
{
    string s_copy = (*s);
    int l = (*s).length();
    for (int i = 0; i < l; i++ )
    {
        (*s)[i] = s_copy[l-1-i];
    }
}

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

void design_maker (string* target_struct, //read
                   string* undefined_rna, //read
                   string* defined_rna,   //write
                   int     dna_flag     )
{
    string undef_dna = (*undefined_rna);
    int l = undef_dna.length();

    if (dna_flag == 1)
    {
        for (int i = 0; i < l; i++)
        {
            if (undef_dna[i] == 'T')
                undef_dna[i]  = 'U';
        }
    }

    ofstream input_w;
    input_w.open(input_design);
    input_w << *target_struct << endl;
    input_w << undef_dna;
    input_w.close();

    system (DESIGN_FUNC);

    string defined_dna;

    ifstream output_r;
    output_r.open(output_design);
    do    {getline(output_r, defined_dna);}
    while ( defined_dna[0] == '%' );
    output_r.close();

    if (dna_flag == 1)
    {
        for (int i = 0; i < l; i++)
        {
            if (defined_dna[i] == 'U')
                defined_dna[i]  = 'T';
        }
    }
    (*defined_rna) = defined_dna;
}

float dG_reader (string* seq)
{
    ofstream w;
    w.open(input_mfe);
    w << *seq;
    w.close();

    system(MFE_FUNC);

    string temporal_str;
    float free_energy;

    ifstream r;
    r.open(output_mfe);

    do    {getline(r, temporal_str);}
    while (temporal_str[0] != '-' &&
           temporal_str[0] != '0'    );

    decimal_numb_reader( &temporal_str, 0, &free_energy );
    return free_energy;
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

    for(int i = 20; i < temp_str.length(); i++)
        aff_str.push_back(temp_str[i]);

    float aff = exp_numb_reader(&aff_str);
    return (aff/1e-06)*100;
}

double fctrl (int n)
{
    if (n == 1 || n == 0) return 1;
    return n*fctrl(n-1);
}

int C_n_k  (int n, int k)
{
    return ( (int)(fctrl(n)/fctrl(n-k)/fctrl(k)) );
}

int gen_comb_norep_lex_init(int *arr, const int n, const int k)
{
int j; //index

//test for special cases
if(k > n)
 return(GEN_ERROR);

if(k == 0)
 return(GEN_EMPTY);

//initialize: arr[0, ..., k - 1] are 0, ..., k - 1
for(j = 0; j < k; j++)
 arr[j] = j;

return(GEN_NEXT);
}

int gen_comb_norep_lex_next(int *arr, const int n, const int k)
{
int j; //index

//easy case, increase rightmost element
if(arr[k - 1] < n - 1)
 {
 arr[k - 1]++;
 return(GEN_NEXT);
 }

//find rightmost element to increase
for(j = k - 2; j >= 0; j--)
 if(arr[j] < n - k + j)
  break;

//terminate if arr[0] == n - k
if(j < 0)
 return(GEN_TERM);

//increase
arr[j]++;

//set right-hand elements
while(j < k - 1)
 {
 arr[j + 1] = arr[j] + 1;
 j++;
 }

return(GEN_NEXT);
}

int comb_maker(int** arr_, int n, int k, int nmb)
{
    int arr[k];
    int           gen_result;         //return value of generation functions
    int           x;                  //iterator
    //initialize
    gen_result = gen_comb_norep_lex_init(arr, n, k);
    for (int i = 0; i < k; i++)
        arr_[0][i] = arr[i];

    if(gen_result == GEN_ERROR)
        return(EXIT_FAILURE);

    int n_numb = 1;
    while(gen_result == GEN_NEXT)
     {
         //for(x = 0; x < k; x++)
            //cout << arr[x] << ' ';

         //cout << endl;

         gen_result = gen_comb_norep_lex_next(arr, n, k); if(n_numb == nmb) return 0;
         for (int i = 0; i < k; i++)
            arr_[n_numb][i] = arr[i];
         n_numb++;
     }

    return(EXIT_SUCCESS);
}

char except_ncltd_symb (char ncltd)
{
    switch(ncltd)
    {
        case 'C':
            return 'D';
        case 'G':
            return 'H';
        case 'A':
            return 'B';
        case 'U':
            return 'V';
        case 'T':
            return 'V';
    }
    return 'E';
}

void out_of_complex_right_brackets_deleter (string* complex_struct)
{
    int l = (*complex_struct).length();
    int brakets_counter;
    for (int i = 0; i < l; i++)
    {
        if ((*complex_struct)[i] == '(')
        {
            brakets_counter = 0;
            for (int j = i; j < l; j++)
            {
                switch((*complex_struct)[j])
                {
                    case '(':
                        brakets_counter++;
                        break;
                    case ')':
                        brakets_counter--;
                }
                if (brakets_counter == 0)
                {
                    (*complex_struct)[j] = '.';
                    break;
                }
            }
        }
    }
}

void target_olig_structure_part_in_complex_reader (string* new_defined_rna,
                                                   string* input_rna,
                                                   string* complex_struct)
{
    ifstream inp_r;
    ofstream inp_w;
    inp_w.open(input_mfe);
    inp_w << (*input_rna) << "+" << (*new_defined_rna);
    inp_w.close();

    system(MFE_FUNC);

    inp_r.open(output_mfe);
    do getline(inp_r, *complex_struct);
    while ((*complex_struct)[0] != '.' &&
           (*complex_struct)[0] != '('   );
    inp_r.close();

    (*complex_struct).erase((*complex_struct).begin(),
                            (*complex_struct).begin() +
                            (*complex_struct).find('+') + 1);
    out_of_complex_right_brackets_deleter (complex_struct);
}

int max_overlap_pstns_aff_to_gr_rdcr (string* best_new_def_rna,
                                      int     oligs_nmb,
                                      float   aff_min_level,
                                      //float   min_aff_to_target_olig,
                                      //string* target_rna,
                                      string* input_rna,
                                      int     rna_length,
                                      anti_aff_olig* AAO_table,
                                      int*    max_overlap_positions,
                                      int     max_overlap_positions_nmb,
                                      string* target_struct,
                                      ofstream* stat_w,
                                      int     dna_flag             )
{
    //string cmplx = (*target_rna) + "+" + (*input_rna);
    float aff;
    float aff_to_target_olig;
    //float initial_cmplx_dG = dG_reader(&cmplx);
    //float cmplx_dG;// = initial_cmplx_dG;
    float antiaff_index;
    float old_best_antiaff_index = 100*oligs_nmb;
    float best_antiaff_index;
    int best_antiaff_index_update_flag;
    //float lost_aff_index;// = (initial_cmplx_dG - cmplx_dG)/initial_cmplx_dG; //lower is better
    //float evaluate_index;// = abs (lost_aff_index * antiaff_index);
    //float max_evaluate_index;
    int mut_variants_nmb;
    int temp_var;
    int*** arr = new int**[1];

    string undefined_rna = (*input_rna);
    //int big_aff_oligs_nmb;
    int break_flag = 0;
    string new_def_rna;
    string old_best_def_rna = (*input_rna);
    for(int mutant_ncltd_numb = 1; mutant_ncltd_numb < max_overlap_positions_nmb; mutant_ncltd_numb++)
    {
        //max_evaluate_index = 0;
        best_antiaff_index = 100*oligs_nmb;
        best_antiaff_index_update_flag = 0;
        mut_variants_nmb = C_n_k(max_overlap_positions_nmb, mutant_ncltd_numb);
        arr[0] = new int* [mut_variants_nmb];
        for (int i = 0; i < mut_variants_nmb; i++)
            arr[0][i] = new int [mutant_ncltd_numb];

        comb_maker (arr[0], max_overlap_positions_nmb, mutant_ncltd_numb, mut_variants_nmb);

        for (int i = 0; i < mut_variants_nmb; i++)
        {
            for (int j = 0; j < mutant_ncltd_numb; j++)
            {
                temp_var = max_overlap_positions[ arr[0][i][j] ];
                undefined_rna[temp_var] = except_ncltd_symb (undefined_rna[temp_var]);// if big cycle CHANGES NEED
            }
            design_maker(target_struct, &undefined_rna, &new_def_rna, dna_flag);
            //aff_to_target_olig = aff_reader(target_rna, &new_def_rna, dna_flag);
            /*
            if (aff_to_target_olig < min_aff_to_target_olig)
            {
                for (int j = 0; j < mutant_ncltd_numb; j++)
                {
                    temp_var = max_overlap_positions[ arr[0][i][j] ];
                    undefined_rna[temp_var] = (*input_rna)[temp_var];
                }
                continue;
            }
            */
            antiaff_index = 0;
            //big_aff_oligs_nmb = 0;
            for (int j = 0; j < oligs_nmb; j++)
            {
                aff = aff_reader(&new_def_rna, &(AAO_table[j].seq), dna_flag);
                if (aff > aff_min_level)
                {
                    antiaff_index += aff; //lower is better
                    //big_aff_oligs_nmb++;
                }
            }
            //if (big_aff_oligs_nmb == 0)
                //big_aff_oligs_nmb = 1;

            if (antiaff_index == 0)
            {
                break_flag = 1;
                (*best_new_def_rna) = new_def_rna;
                break; //we consider this variant as good
            }

            //antiaff_index = (100 - antiaff_index/big_aff_oligs_nmb)/100; //higher is better
            //cmplx = (*target_rna) + "+" + new_def_rna;
            //cmplx_dG = dG_reader(&cmplx);
            //lost_aff_index = (initial_cmplx_dG - cmplx_dG)/initial_cmplx_dG;// it's considered that cmplx_dG > initial_cmplx_dG, lower is better
            //evaluate_index = (1 - lost_aff_index) / antiaff_index; //higher is better

            if (antiaff_index < best_antiaff_index)
            {
                best_antiaff_index = antiaff_index;
                best_antiaff_index_update_flag = 1;
                (*best_new_def_rna) = new_def_rna;
            }
            //if (mutant_ncltd_numb < 3)
                //(*stat_w) << undefined_rna << ' ' << lost_aff_index << ' ' << antiaff_index << ' ' << evaluate_index << endl << new_def_rna << endl;

            for (int j = 0; j < mutant_ncltd_numb; j++)
            {
                temp_var = max_overlap_positions[ arr[0][i][j] ];
                undefined_rna[temp_var] = (*input_rna)[temp_var]; //if big cycle best undefined_rna should be saved
            }
        }
        //(*stat_w) << endl << best_antiaff_index << endl << max_evaluate_index << endl << endl;
        if (break_flag)
        {
            for (int i = 0; i < mut_variants_nmb; i++)
                delete [] arr[0][i];
            break;
        }
        if ( (old_best_antiaff_index < best_antiaff_index) ||
             (best_antiaff_index_update_flag == 0)           )
        {
            for (int i = 0; i < mut_variants_nmb; i++)
                delete [] arr[0][i];
            (*best_new_def_rna) = old_best_def_rna;
            break_flag = 1;
            break;
        }
        old_best_def_rna = (*best_new_def_rna);
        old_best_antiaff_index = best_antiaff_index;

        for (int i = 0; i < mut_variants_nmb; i++)
            delete [] arr[0][i];
    }
    delete [] arr[0];
    delete [] arr;
    return break_flag;
}

void data_to_gr_reader_and_max_overlap_counter (int     oligs_nmb,
                                                float   aff_min_level,
                                                string* input_rna,
                                                int     rna_length,
                                                anti_aff_olig* AAO_table,
                                                int*    max_overlap_positions,
                                                int*    max_overlap_positions_nmb,
                                                ofstream* stat_w,
                                                int     dna_flag)
{
    int overlap_arr [rna_length];
    int max_overlap_counter = 0;
    int overlap_counter = 0;

    for (int i = 0; i < oligs_nmb; i++)
    {
        AAO_table[i].aff = aff_reader(input_rna, &(AAO_table[i].seq), dna_flag);
        target_olig_structure_part_in_complex_reader( input_rna,
                                                      &(AAO_table[i].seq),
                                                      &(AAO_table[i].cmplx_struct_part) );//here will be information about what ncltds in input_rna are connected with ncltds of AAO_table[i].seq
        //(*stat_w) << AAO_table[i].cmplx_struct_part << ' ' << AAO_table[i].aff << endl;
    }
    //(*stat_w) << endl;
    for (int j = 0; j < rna_length; j++)
    {
        overlap_counter = 0;
        for (int i = 0; i < oligs_nmb; i++)
        {
            if (AAO_table[i].aff > aff_min_level        &&
                AAO_table[i].cmplx_struct_part[j] == ')'  )
            {
                overlap_counter++;
            }
        }
        overlap_arr[j] = overlap_counter;
        //(*stat_w) << overlap_arr[j] << ' ';

        if (overlap_counter > max_overlap_counter)
        {
            max_overlap_counter = overlap_counter;
        }
    }
    //(*stat_w) << endl;
    (*max_overlap_positions_nmb) = 0;
    for (int j = 0; j < rna_length; j++)
    {
        if (overlap_arr[j] == max_overlap_counter)
        {
            max_overlap_positions [(*max_overlap_positions_nmb)] = j;
            //(*stat_w) << j << ' ';
            (*max_overlap_positions_nmb)++;
        }
    }
    //(*stat_w) << endl << (*max_overlap_positions_nmb) << endl;
}

int main()
{
    string temp_str;
    float temp_fl;

    ifstream inp_r;
    inp_r.open(input_file);

    getline(inp_r, temp_str);
    decimal_numb_reader(&temp_str, (int)temp_str.find(' ') + 1, &temp_fl);
    int dna_flag = (int)temp_fl;
/*
    string target_rna;
    getline(inp_r, target_rna);
    target_rna.erase(target_rna.begin(),
                     target_rna.begin() +
                     target_rna.find(' ') + 1);
*/
    string input_rna;
    getline(inp_r, input_rna);
    input_rna.erase(input_rna.begin(),
                    input_rna.begin() +
                    input_rna.find(' ') + 1);//cout << target_rna << endl;
/*
    float min_aff_to_target_olig;
    getline(inp_r, temp_str);
    decimal_numb_reader(&temp_str, (int)temp_str.find(' ') + 1, &min_aff_to_target_olig);
*/
    getline(inp_r, temp_str);
    decimal_numb_reader(&temp_str, (int)temp_str.find(' ') + 1, &temp_fl);
    int oligs_nmb = (int)temp_fl;

    float aff_min_level;
    getline(inp_r, temp_str);
    decimal_numb_reader(&temp_str, (int)temp_str.find(' ') + 1, &aff_min_level);

    int rna_length = input_rna.length();//
    anti_aff_olig AAO_table [oligs_nmb];//

    for (int i = 0; i < oligs_nmb; i++)
        getline(inp_r, AAO_table[i].seq);

    inp_r.close();

    int max_overlap_positions [rna_length];//+
    int max_overlap_positions_nmb = 0;//+

    string target_struct;//+
    for (int i = 0; i < rna_length; i++)
        target_struct.push_back('.');

    string best_new_def_rna = input_rna;
    int loop_counter = 0;
    int temp_int;

    ofstream stat_w;
    //stat_w.open(statistic);

    do
    {
        /*
        if(loop_counter)
        {
            cout << "Enter any number to continue" << endl;
            cin >> temp_int;
        }
        */
        data_to_gr_reader_and_max_overlap_counter (oligs_nmb,
                                                   aff_min_level,
                                                   &best_new_def_rna,
                                                   rna_length,
                                                   AAO_table,
                                                   max_overlap_positions,
                                                   &max_overlap_positions_nmb,
                                                   &stat_w,
                                                   dna_flag);
        //cout << "Enter any number to continue" << endl;
        //cin >> temp_int;
        loop_counter++;
        if (loop_counter > max_loop_repetition)
            break;
    }
    while ( max_overlap_pstns_aff_to_gr_rdcr (&best_new_def_rna,
                                              oligs_nmb,
                                              aff_min_level,
                                              //min_aff_to_target_olig,
                                              //&target_rna,
                                              &input_rna,
                                              rna_length,
                                              AAO_table,
                                              max_overlap_positions,
                                              max_overlap_positions_nmb,
                                              &target_struct,
                                              &stat_w,
                                              dna_flag                   ) == 0 );
    stat_w.close();

    //float aff_to_target_olig = aff_reader(&best_new_def_rna, &target_rna);
    ofstream result;
    result.open(output_file);
    result << best_new_def_rna << endl
           << "dG (kkal/mol) = " << dG_reader(&best_new_def_rna) << endl;
           //<< "aff to target olig = " << aff_to_target_olig;
    //if (aff_to_target_olig > min_aff_to_target_olig)
         //result << " +" << endl;
    //else result << " -" << endl;
    result << endl << "aff to every olig from group:" << endl;
    float aff_to_olig_from_group;
    for (int i = 0; i < oligs_nmb; i++)
    {
        aff_to_olig_from_group = aff_reader(&best_new_def_rna, &(AAO_table[i].seq));
        result << (i + 1) << ") " << aff_to_olig_from_group;

        if (aff_to_olig_from_group < aff_min_level)
             result << " +" << endl;
        else result << " -" << endl;
    }
    result.close();
    system(CLEAN_FUNC);
    return 0;
}
