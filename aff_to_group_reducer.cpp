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

const int   absolute_aff_changing = 1;
#define MFE_FUNC "\
    #/bin/bash \n\
    export NUPACKHOME=/home/valentin/work_backup/nupack \n\
    export CODDING_HOME=/home/valentin/work_backup/codding \n\
    $NUPACKHOME/bin/mfe $CODDING_HOME/input/input_seq \n\
"
#define DESIGN_FUNC "\
    #/bin/bash \n\
    export NUPACKHOME=/home/valentin/work_backup/nupack \n\
    export CODDING_HOME=/home/valentin/work_backup/codding \n\
    $NUPACKHOME/bin/design -loadseed $CODDING_HOME/input/input_seq \n\
"

#define CONCENTRATIONS_FUNC "\
    #/bin/bash \n\
    export NUPACKHOME=/home/valentin/work_backup/nupack \n\
    export CODDING_HOME=/home/valentin/work_backup/codding \n\
    $NUPACKHOME/bin/complexes -T 23 -material rna -quiet $CODDING_HOME/input/complex \n\
    $NUPACKHOME/bin/concentrations -quiet $CODDING_HOME/input/complex \n\
"
struct anti_aff_olig
{
    string seq;
    string cmplx_struct_part;
    float aff;
};

struct hairpin_nucleotide_links
{
    float first_nucleotide_number_in_link;
    float second_nucleotide_number_in_link;
    float distance_between_nucleotides;
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
                   string* defined_rna   )//write
{
    ofstream input_w;
    input_w.open(input_design);
    input_w << *target_struct << endl;
    input_w << *undefined_rna;
    input_w.close();

    system (DESIGN_FUNC);

    ifstream output_r;
    output_r.open(output_design);
    do    {getline(output_r, *defined_rna);}
    while ( (*defined_rna)[0] == '%' );
    output_r.close();
}

void links_reader(hairpin_nucleotide_links* HNL,
                  string*                   link_string)
{
    decimal_numb_reader(link_string,
                        decimal_numb_reader(link_string, 0, &(*HNL).first_nucleotide_number_in_link),
                        &(*HNL).second_nucleotide_number_in_link);

    (*HNL).distance_between_nucleotides = (*HNL).second_nucleotide_number_in_link -
                                          (*HNL).first_nucleotide_number_in_link;
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


float dG_and_struct_reader (string* seq,
                            string* seq_struct,
                            hairpin_nucleotide_links*  HNL,
                            int*    number_of_links)
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
    int i = 0;
    getline(r, *seq_struct);
    do
    {
        getline(r, temporal_str);
        links_reader(&HNL[i],&temporal_str);
        i++;
    }
    while (temporal_str[0] != '%');
    *number_of_links = i - 1;

    r.close();
    return free_energy;
}


float aff_reader (string* seq_1, string* seq_2)
{
    ofstream w;
    w.open(conc_input);
    w << 2        << endl
      << (*seq_1) << endl
      << (*seq_2) << endl
      << 2;
    w.close();
    system(CONCENTRATIONS_FUNC);

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
    }
    return 'E';
}

int careful_aff_reducer( string* complex_struct,
                         string* undefined_rna,
                         string* new_defined_rna,
                         string* target_struct,
                         string* input_rna,
                         float   target_aff,
                         float   target_aff_accuracy,
                         float*  new_aff)
{
    int ncltds_in_cmplx_nmbs [(*complex_struct).length()];
    int cmplx_nmb_of_links = 0;
    for (int i = 0; i < (*complex_struct).length(); i++)
    {
        if ( (*complex_struct)[i] == ')' )
        {
            ncltds_in_cmplx_nmbs[cmplx_nmb_of_links] = i;
            cmplx_nmb_of_links++;
        }
    }
    int*** arr = new int**[1];
    int mut_variants_nmb;
    int temp_var;
    string new_def_rna;
    string best_def_rna;
    string old_best_def_rna;
    float d_aff;
    float min_d_aff = target_aff;
    float old_min_d_aff = target_aff;
    int break_flag = 0;

    for(int mutant_ncltd_numb = 1; ; mutant_ncltd_numb++)
    {
        mut_variants_nmb = C_n_k(cmplx_nmb_of_links, mutant_ncltd_numb);
        arr[0] = new int* [mut_variants_nmb];
        for (int i = 0; i < mut_variants_nmb; i++)
            arr[0][i] = new int [mutant_ncltd_numb];

        comb_maker (arr[0], cmplx_nmb_of_links, mutant_ncltd_numb, mut_variants_nmb);

        for (int i = 0; i < mut_variants_nmb; i++)
        {
            for (int j = 0; j < mutant_ncltd_numb; j++)
            {
                temp_var = ncltds_in_cmplx_nmbs[ arr[0][i][j] ];
                (*undefined_rna)[temp_var] = except_ncltd_symb ((*undefined_rna)[temp_var]);
            }
            design_maker(target_struct, undefined_rna, &new_def_rna);
            d_aff = aff_reader(input_rna, &new_def_rna)  - target_aff;

            if ( abs(d_aff) < abs(min_d_aff) )
            {
                min_d_aff = d_aff;
                best_def_rna = new_def_rna;
                if ( abs(min_d_aff) < target_aff_accuracy)
                {
                    break_flag = 1;
                    break;
                }
            }
            for (int j = 0; j < mutant_ncltd_numb; j++)
            {
                temp_var = ncltds_in_cmplx_nmbs[ arr[0][i][j] ];
                (*undefined_rna)[temp_var] = (*new_defined_rna)[temp_var];
            }
        }
        if (break_flag)
        {
            for (int i = 0; i < mut_variants_nmb; i++)
                delete [] arr[0][i];
            (*new_defined_rna) = best_def_rna;
            (*new_aff) = target_aff + min_d_aff;
            break;
        }
        if ( abs(min_d_aff) > abs(old_min_d_aff) )
        {
            for (int i = 0; i < mut_variants_nmb; i++)
                delete [] arr[0][i];
            (*new_defined_rna) = old_best_def_rna;
            (*new_aff) = target_aff + old_min_d_aff;
            break;
        }
        old_best_def_rna = best_def_rna;
        old_min_d_aff = min_d_aff;
        min_d_aff = target_aff;
        for (int i = 0; i < mut_variants_nmb; i++)
            delete [] arr[0][i];
    }
    delete [] arr[0];
    delete [] arr;
    return 1;
}

int aff_reducer_2(string* input_rna,
                  string* complex_struct,
                  string* new_defined_rna,
                  string* undefined_rna,
                  string* target_struct,
                  float*  new_aff,
                  float   aff)
{
    string last_variant = (*new_defined_rna);
    for (int i = 0; i < (*complex_struct).length(); i++)
    {
        if ( (*complex_struct)[i] == ')' )
        {
            (*undefined_rna)[i] = except_ncltd_symb ((*undefined_rna)[i]);
            design_maker(target_struct,
                         undefined_rna,
                         new_defined_rna);
            *new_aff = aff_reader (input_rna, new_defined_rna);
            if ( aff - *new_aff > absolute_aff_changing ) // aff started change significantly
            {
                (*new_defined_rna) = last_variant;
                return 1;
            }
            last_variant = (*new_defined_rna);
        }
    }
    return 0; // aff didn't start change significantly
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

int main()
{   // DESIGN_FUNC MFE_FUNC CONCENTRATIONS_FUNC
    //system(MFE_FUNC);
    int oligs_nmb = 5;
    float aff_min_level = 20;
    float min_aff_to_target_olig = 80;
    string target_rna = "CAAGCCACAGAUCCCCACGUUCAUUAAGCCGCC";
    string input_rna  = "GCCGGCUCAAUGAGCGUGGGGAUCCCUGGCGUG";
    anti_aff_olig AAO_table [oligs_nmb];
    ifstream inp_r;
    inp_r.open("./input/new_input")

    for (int i = 0; i < oligs_nmb; i++)
    {
        getline(inp_r, AAO_table[i].seq);
        AAO_table[i].aff = aff_reader(&input_rna, &(AAO_table[i].seq));
        target_olig_structure_part_in_complex_reader( &(AAO_table[i].seq),
                                                      &input_rna,
                                                      &(AAO_table[i].cmplx_struct_part) );
    }
    inp_r.close();

    int rna_length = input_rna.length();
    int overlap_arr [rna_length];
    int max_overlap_positions [rna_length];
    int max_overlap_positions_nmb = 0;
    int max_overlap_counter = 0;
    int overlap_counter = 0;
    string cmplx = target_rna + "+" + input_rna;
    float aff;
    float aff_to_target_olig;
    float initial_cmplx_dG = dG_reader(&cmplx);
    float cmplx_dG = initial_cmplx_dG;
    float antiaff_index = 0; //lower is better
    float lost_aff_index = (initial_cmplx_dG - cmplx_dG)/initial_cmplx_dG; //lower is better
    float evaluate_index = abs (lost_aff_index * antiaff_index);//lower is better
    int mut_variants_nmb;
    int*** arr = new int**[1];

    string undefined_rna = input_rna;
    string new_def_rna;
    string target_struct;
    for (int i = 0; i < rna_length; i++)
        target_struct.push_back('.');
/*
    for (int i = 0; i < oligs_nmb; i++)
    {
        aff = aff_reader(&new_def_rna, &(AAO_table[i].seq) );
        if (aff > aff_min_level)
        {
            antiaff_index += aff;
        }
    }
*/
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

        if (overlap_counter > max_overlap_counter)
        {
            max_overlap_counter = overlap_counter;
        }
    }
    for (int j = 0; j < rna_length; j++)
    {
        if (overlap_arr[j] == max_overlap_counter)
        {
            max_overlap_positions [max_overlap_positions_nmb] = j;
            max_overlap_positions_nmb++;
        }
    }
    for(int mutant_ncltd_numb = 1; mutant_ncltd_numb < max_overlap_positions_nmb; mutant_ncltd_numb++)
    {
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
                undefined_rna[temp_var] = except_ncltd_symb (undefined_rna[temp_var]);
            }
            design_maker(&target_struct, &undefined_rna, &new_def_rna);
            aff_to_target_olig = aff_reader(&target_rna, &new_def_rna);
            if (aff_to_target_olig < min_aff_to_target_olig)
                continue;
            antiaff_index = 0;
            for (int j = 0; j < oligs_nmb; j++)
            {
                aff = aff_reader(&new_def_rna, &(AAO_table[j].seq) );
                if (aff > aff_min_level)
                {
                    antiaff_index += aff;
                }
            }
            cmplx = target_rna + "+" + input_rna;
            cmplx_dG = dG_reader(&cmplx);
            lost_aff_index = (initial_cmplx_dG - cmplx_dG)/initial_cmplx_dG;
            evaluate_index = abs (lost_aff_index * antiaff_index);

            if ( abs(d_aff) < abs(min_d_aff) )
            {
                min_d_aff = d_aff;
                best_def_rna = new_def_rna;
                if ( abs(min_d_aff) < target_aff_accuracy)
                {
                    break_flag = 1;
                    break;
                }
            }
            for (int j = 0; j < mutant_ncltd_numb; j++)
            {
                temp_var = ncltds_in_cmplx_nmbs[ arr[0][i][j] ];
                (*undefined_rna)[temp_var] = (*new_defined_rna)[temp_var];
            }
        }
        if (break_flag)
        {
            for (int i = 0; i < mut_variants_nmb; i++)
                delete [] arr[0][i];
            (*new_defined_rna) = best_def_rna;
            (*new_aff) = target_aff + min_d_aff;
            break;
        }
        if ( abs(min_d_aff) > abs(old_min_d_aff) )
        {
            for (int i = 0; i < mut_variants_nmb; i++)
                delete [] arr[0][i];
            (*new_defined_rna) = old_best_def_rna;
            (*new_aff) = target_aff + old_min_d_aff;
            break;
        }
        old_best_def_rna = best_def_rna;
        old_min_d_aff = min_d_aff;
        min_d_aff = target_aff;
        for (int i = 0; i < mut_variants_nmb; i++)
            delete [] arr[0][i];
    }
    undefined_rna[max_link_overlap_ncltd_nmb] = except_ncltd_symb[]

    return 0;
}
