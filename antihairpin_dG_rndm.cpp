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

const int   max_loop_size = 10;
const int   absolute_aff_changing = 1;
int pause;
//const int   min_level_of_aff_relative_changing = 0.03; //should be < 1 and influence the way how diff variants are compared
//const float diff_in_energy_limit = 0.001;

int   rna_length = 0;

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
    return 0;
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

struct hairpin_nucleotide_links
{
    float first_nucleotide_number_in_link;
    float second_nucleotide_number_in_link;
    float distance_between_nucleotides;
};

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
    $NUPACKHOME/bin/complexes -T 23 -material rna $PWD/input/complex \n\
    $NUPACKHOME/bin/concentrations $PWD/input/complex \n\
"

#define CONCENTRATIONS_FUNC_DNA "\
    #/bin/bash \n\
    export NUPACKHOME=/home/valentin/work_backup/nupack \n\
    $NUPACKHOME/bin/complexes -T 23 -material dna $PWD/input/complex \n\
    $NUPACKHOME/bin/concentrations $PWD/input/complex \n\
"

#define CLEAN_FUNC "\
    #/bin/bash \n\
    rm $PWD/input/var* \n\
    rm $PWD/input/input_seq.in \n\
    rm $PWD/input/input_seq.fold \n\
    rm $PWD/input/input_seq.summary \n\
    rm $PWD/input/input_seq.mfe \n\
    rm $PWD/input/complex.in \n\
    rm $PWD/input/complex.eq \n\
    rm $PWD/input/complex.cx \n\
"

void links_reader(hairpin_nucleotide_links* HNL,
                  string*                   link_string)
{
    decimal_numb_reader(link_string,
                        decimal_numb_reader(link_string, 0, &(*HNL).first_nucleotide_number_in_link),
                        &(*HNL).second_nucleotide_number_in_link);

    (*HNL).distance_between_nucleotides = (*HNL).second_nucleotide_number_in_link -
                                          (*HNL).first_nucleotide_number_in_link;
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
    getline(r, *seq_struct);

    if (free_energy == 0)
    {
        *number_of_links = 0;
        return free_energy;
    }

    int i = 0;
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


int max_arr_elem_numb (float* arr, int arr_length)
{
    int max_elem_numb = 0;
    float max_elem = arr[0];
    for (int i = 1; i < arr_length; i++)
    {
        if (arr[i] > max_elem)
        {
            max_elem = arr[i];
            max_elem_numb = i;
        }
    }
    return max_elem_numb;
}

int link_with_max_distance (hairpin_nucleotide_links* HNL,
                            int*                      number_of_links)
{
    float max_link_distance  = 0.0;
    int   max_element_number = 0;

    for (int i = 0; i < *number_of_links; i++)
    {
        if ( max_link_distance < HNL[i].distance_between_nucleotides )
        {
             max_link_distance = HNL[i].distance_between_nucleotides;
             max_element_number = i;
        }
    }
    HNL[max_element_number].distance_between_nucleotides = 0.0;
    return max_element_number;
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

int cmplx_dG_reducer_1(string* input_rna,      //it's like simple total antihairpin (additional)
                  string* new_defined_rna,
                  string* undefined_rna,
                  string* target_struct,
                  int*    number_of_links,
                  float   target_cmplx_dG,
                  hairpin_nucleotide_links*  HNL,
                  int     dna_flag)
{
    string last_variant = (*new_defined_rna);
    string cmplx;
    float new_cmplx_dG;

    for (int i = 0; i < *number_of_links; i++)
    {
        if ( (*undefined_rna)[HNL[i].first_nucleotide_number_in_link - 1] == 'N' )
            (*undefined_rna)[HNL[i].second_nucleotide_number_in_link - 1] = 'N';
        else (*undefined_rna)[HNL[i].first_nucleotide_number_in_link - 1] = 'N';

        design_maker(target_struct,
                     undefined_rna,
                     new_defined_rna,
                     dna_flag);

        cmplx = (*input_rna) + "+" + (*new_defined_rna);
        new_cmplx_dG = dG_reader(&cmplx);
        if ( new_cmplx_dG > target_cmplx_dG )
        {
            (*new_defined_rna) = last_variant;
            return 1;
        }
        if ( dG_reader(new_defined_rna) == 0 )
            return 1;
        last_variant = (*new_defined_rna);
    }

    string temp_str;
    if (dG_and_struct_reader (new_defined_rna,
                              &temp_str,
                              HNL,
                              number_of_links)  < 0)
    {
        if (cmplx_dG_reducer_1(input_rna,
                          new_defined_rna,
                          undefined_rna,
                          target_struct,
                          number_of_links,
                          target_cmplx_dG,
                          HNL,
                          dna_flag)              )
        return 1;
    }
    return 0; // cmplx_dG didn't start change significantly
}

float cmplx_dG_reducer(float   target_cmplx_dG,
                       float   target_cmplx_dG_accuracy,
                       string* defined_rna,          // here is resulting rna
                       string* input_rna)
{
    string cmplx = (*input_rna) + "+" + (*defined_rna);
    string new_cmplx;
    string best_new_cmplx;

    float new_cmplx_dG;
    float best_cmplx_dG;
    float diff_cmplx_dG;
    float min_diff_cmplx_dG;

    int inp_rna_length = (*defined_rna).length();

    for (int cut_ncltd_nmb = 1; cut_ncltd_nmb <= inp_rna_length; cut_ncltd_nmb++)
    {
        min_diff_cmplx_dG = abs (100 * target_cmplx_dG) + 100;
        for (int i = 0; i <= cut_ncltd_nmb; i++)
        {
            new_cmplx = cmplx;
            new_cmplx.erase( new_cmplx.begin() + inp_rna_length + 1,
                             new_cmplx.begin() + inp_rna_length + 1 + i);

            new_cmplx.erase( new_cmplx.end() - cut_ncltd_nmb + i,
                             new_cmplx.end()                     );

            new_cmplx_dG  = dG_reader(&new_cmplx);
            diff_cmplx_dG = abs (new_cmplx_dG - target_cmplx_dG);

            if (diff_cmplx_dG < min_diff_cmplx_dG)
            {
                min_diff_cmplx_dG = diff_cmplx_dG;
                best_cmplx_dG = new_cmplx_dG;
                best_new_cmplx = new_cmplx;
            }
        }
        if (min_diff_cmplx_dG <= target_cmplx_dG_accuracy ||
            best_cmplx_dG > target_cmplx_dG)
        {
            break;
        }
    }
    best_new_cmplx.erase( best_new_cmplx.begin(), best_new_cmplx.begin() + inp_rna_length + 1);

    (*defined_rna) = best_new_cmplx;
    return best_cmplx_dG;
}

int careful_anti_hairpin(string* input_rna,
                         float   target_hairpin_dG,
                         float   target_dG_accuracy,
                         float   target_cmplx_dG,
                         float   target_cmplx_dG_accuracy,
                         string* old_defined_rna,
                         string* old_undefined_rna,
                         int     number_of_links,
                         hairpin_nucleotide_links* HNL,
                         int     dna_flag,
                         int     max_dist_last_HNL_element_number,
                         int     max_dist_HNL_element_number,
                         float   initial_dG)
{
    int ncltds_nmb_in_links_gr = 2 * (max_dist_last_HNL_element_number - max_dist_HNL_element_number + 1);
    int ncltds_nmbs_in_links_gr [ncltds_nmb_in_links_gr];
    int l = 0;
    for (int k = max_dist_last_HNL_element_number; k >= max_dist_HNL_element_number; k--)
    {
        ncltds_nmbs_in_links_gr [l]     = HNL[k].first_nucleotide_number_in_link  - 1;
        ncltds_nmbs_in_links_gr [l + 1] = HNL[k].second_nucleotide_number_in_link - 1;
        l += 2;
    }
    string target_struct;

    for (int i = 0; i < (*old_undefined_rna).length(); i++)

        target_struct.push_back('.');

    int*** arr = new int**[1];
    int mut_variants_nmb;
    int temp_var;
    int break_flag = 0;

    string new_undef_rna;
    string new_def_rna;
    string best_def_rna;
    string old_best_def_rna = (*old_defined_rna);
    string cmplx;
    float dG;
    float best_dG;
    float old_best_dG = initial_dG;
    float cmplx_dG;
    float best_cmplx_dG;

    //ofstream output_w;
    //output_w.open("./output/careful_antihairpin_stat");

    for(int mutant_ncltd_numb = 1; mutant_ncltd_numb <= ncltds_nmb_in_links_gr; mutant_ncltd_numb++)
    {
        best_cmplx_dG = 1;
        best_dG = initial_dG;
        mut_variants_nmb = C_n_k( ncltds_nmb_in_links_gr, mutant_ncltd_numb );
        arr[0] = new int* [mut_variants_nmb];
        for (int i = 0; i < mut_variants_nmb; i++)
            arr[0][i] = new int [mutant_ncltd_numb];

        comb_maker (arr[0], ncltds_nmb_in_links_gr, mutant_ncltd_numb, mut_variants_nmb);

        //output_w << endl << "mutant_ncltd_numb = " << mutant_ncltd_numb << endl << endl;
        for (int i = 0; i < mut_variants_nmb; i++)
        {
            new_undef_rna = (*old_undefined_rna);
            for (int j = 0; j < mutant_ncltd_numb; j++)
            {
                temp_var = ncltds_nmbs_in_links_gr[ arr[0][i][j] ];
                new_undef_rna[temp_var] = 'N';
            }
            design_maker(&target_struct, &new_undef_rna, &new_def_rna, dna_flag);

            cmplx = new_def_rna + "+" + (*input_rna);
            cmplx_dG = dG_reader(&cmplx);

            if ( cmplx_dG - target_cmplx_dG > target_cmplx_dG_accuracy )
                continue;

            dG = dG_reader(&new_def_rna);

            //output_w << '\t' << new_undef_rna << '\t' << cmplx_dG << '\t' << dG << endl;

            if (dG > best_dG)
            {
                //best_evaluate_index = evaluate_index;
                best_def_rna = new_def_rna;
                best_cmplx_dG = cmplx_dG;
                best_dG = dG;
                if ( abs (best_dG - target_hairpin_dG) < target_dG_accuracy )
                {
                    break_flag = 1;
                    break;
                }
            }
        }
        //output_w << endl << best_def_rna << '\t' << best_cmplx_dG << '\t' << best_dG << endl;

        if ( best_dG < old_best_dG || //previous set of mutations gave better result than current
             best_cmplx_dG == 1         || //there are no mutations in current set that save required cmplx_dG level
             break_flag == 1         )//success
        {
            for (int i = 0; i < mut_variants_nmb; i++)
                delete [] arr[0][i];

            if (break_flag == 1)
                old_best_def_rna = best_def_rna;

            break;
        }
        old_best_def_rna = best_def_rna;
        old_best_dG = best_dG;

        for (int i = 0; i < mut_variants_nmb; i++)
            delete [] arr[0][i];
    }
    (*old_defined_rna) = old_best_def_rna;

    //output_w.close();
    delete [] arr[0];
    delete [] arr;
    return 1;
}

void result_writer (ofstream* output_w,
                    string*   olig_seq,
                    float     cmplx_dG,
                    float     dG,
                    float     aff,
                    char      plus_or_minus)
{
    (*output_w) << (*olig_seq) << '\t'
                <<     cmplx_dG     << '\t'
                <<     dG      << '\t'
                <<     aff      << '\t'
                <<plus_or_minus<< endl;
}

int end_program_check (float     dG,
                       float     target_hairpin_dG,
                       float     target_dG_accuracy,
                       float     cmplx_dG,
                       float     target_cmplx_dG,
                       float     target_cmplx_dG_accuracy,
                       string*   defined_rna,
                       string*   old_defined_rna,
                       ofstream* output_w,
                       string*   input_rna,
                       int       new_number_of_links,
                       int       number_of_links,
                       hairpin_nucleotide_links*  new_HNL,
                       hairpin_nucleotide_links*  HNL,
                       int       dna_flag,
                       int       max_dist_last_HNL_element_number,
                       int       max_dist_HNL_element_number,
                       float     initial_dG)
{
    float new_cmplx_dG,
          new_dG;
    char  plus_or_minus;
    string cmplx,
           undefined_rna,
           target_struct;
    for (int i = 0; i < (*defined_rna).length(); i++)
        target_struct.push_back('.');

    int nmb_of_links = new_number_of_links;

    if ( cmplx_dG - target_cmplx_dG > target_cmplx_dG_accuracy )
    {
        /*
        careful_anti_hairpin(input_rna,
                             target_hairpin_dG,
                             target_dG_accuracy,
                             target_cmplx_dG,
                             target_cmplx_dG_accuracy,
                             old_defined_rna,
                             old_undefined_rna,
                             number_of_links,
                             HNL,
                             dna_flag,
                             max_dist_last_HNL_element_number,
                             max_dist_HNL_element_number,
                             initial_dG);*/

        cmplx = (*input_rna) + "+" + (*old_defined_rna);
        new_cmplx_dG = dG_reader(&cmplx);
        new_dG  = dG_reader(old_defined_rna);

        if ( abs(new_cmplx_dG - target_cmplx_dG) < target_cmplx_dG_accuracy &&
             ( abs(new_dG  - target_hairpin_dG) <= target_dG_accuracy ||
               new_dG > target_hairpin_dG                               )  )

             plus_or_minus = '+';
        else plus_or_minus = '-';

        result_writer (output_w,
                       old_defined_rna,
                       new_cmplx_dG,
                       new_dG,
                       aff_reader(old_defined_rna, input_rna, dna_flag),
                       plus_or_minus);
        return 1;
    }
    if ( abs(dG - target_hairpin_dG) <= target_dG_accuracy
                                                          ||
         dG > target_hairpin_dG                             )
    {
        if ( target_cmplx_dG - cmplx_dG > target_cmplx_dG_accuracy )
        {
/*            undefined_rna = (*defined_rna);
// THIS CODE
            cmplx_dG_reducer_1(input_rna,
                                     defined_rna,
                                     &undefined_rna,
                                     &target_struct,
                                     &nmb_of_links,
                                     target_cmplx_dG,
                                     new_HNL,
                                     dna_flag);

            cmplx = (*input_rna) + "+" + (*defined_rna);
            new_cmplx_dG = dG_reader(&cmplx);*/
// MAY CAUSE PROBLEMS
            if ( target_cmplx_dG - new_cmplx_dG > target_cmplx_dG_accuracy )
            {
                new_cmplx_dG = cmplx_dG_reducer(target_cmplx_dG,
                                                target_cmplx_dG_accuracy,
                                                defined_rna,
                                                input_rna);
            }
            new_dG = dG_reader(defined_rna);

            if ( abs(new_cmplx_dG - target_cmplx_dG) < target_cmplx_dG_accuracy )
                 plus_or_minus = '+';
            else plus_or_minus = '-';

            result_writer (output_w,
                           defined_rna,
                           new_cmplx_dG,
                           new_dG,
                           aff_reader(defined_rna, input_rna, dna_flag),
                           plus_or_minus);
        }
        else
        {
            if ( abs(cmplx_dG - target_cmplx_dG) < target_cmplx_dG_accuracy )
                 plus_or_minus = '+';
            else plus_or_minus = '-';

            result_writer (output_w,
                           defined_rna,
                           cmplx_dG,
                           dG,
                           aff_reader(defined_rna, input_rna, dna_flag),
                           plus_or_minus);
        }
        return 1;
    }
    return 0;
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

char except_defined_ncltd (char ncltd, int except_ncltd_nmb, int dna_flag) //except_ncltd_nmb can be 0, 1 or 2
{
    char ncltds [4] = {'A', 'U', 'G', 'C'};

    if (dna_flag)
        ncltds[1] = 'T';

    char except_ncltds [3];
    int j = 0;
    for (int i = 0; i < 4; i++)
    {
        if (ncltd != ncltds[i])
        {
            except_ncltds[j] = ncltds[i];
            j++;
        }
    }
    return except_ncltds[except_ncltd_nmb];
}

int anti_hairpin  (string* defined_rna,
                   string* input_rna,
                   float   initial_dG,
                   float   initial_cmplx_dG,
                   hairpin_nucleotide_links*  HNL,
                   int*    number_of_links,
                   float   target_hairpin_dG,
                   float   target_dG_accuracy,
                   float   target_cmplx_dG,
                   float   target_cmplx_dG_accuracy,
                   ofstream* output_w,
                   ofstream* stat_w,
                   float break_flag,
                   int     dna_flag)
                   //ofstream* bigstat_w,
                   //int loop_counter)
{
    int max_dist_HNL_element_number = 0;
    int lack_of_link_groups = *number_of_links; //*bigstat_w << "number_of_links" << '\t' << *number_of_links << endl;
    int i = 0;
    float cmplx_dG_, dG_;
    string complex_seq;

    while (i - max_dist_HNL_element_number == 0) //number of sequence links (one after another)
    {
        if (lack_of_link_groups <= 1)
        {
            complex_seq = (*input_rna) + "+" + (*defined_rna);
            cmplx_dG_ = dG_reader(&complex_seq);
            dG_ = dG_reader(defined_rna);

            (*output_w) << (*defined_rna) << '\t'
                        << cmplx_dG_ << '\t'
                        << dG_  << '\t'
                        << aff_reader(defined_rna, input_rna, dna_flag)  << '\t';
            if ( (abs( cmplx_dG_ - target_cmplx_dG ) < target_dG_accuracy     ||   // HERE CORRECTIONS NEED
                 cmplx_dG_ < target_cmplx_dG)                                 &&
                 (abs(dG_ - target_hairpin_dG) < target_dG_accuracy ||
                 dG_ > target_hairpin_dG)                              )
                 (*output_w) << '\t' << '+' << endl;
            else (*output_w) << '\t' << '-' << endl;

            return 1;
        }
        max_dist_HNL_element_number =
        link_with_max_distance(HNL, number_of_links);
        for (i = max_dist_HNL_element_number;

             ( (HNL[i].first_nucleotide_number_in_link  + 1) ==
               HNL[i+1].first_nucleotide_number_in_link        )
                                                                &&
             ( (HNL[i].second_nucleotide_number_in_link - 1) ==
               HNL[i+1].second_nucleotide_number_in_link       );

             i++) // i shows the last HNL elem numb in current group of max distant links
         {}
         lack_of_link_groups--;
    } //*bigstat_w << "i" << '\t' << i << endl;
    string new_defined_rna,
           best_defined_rna,
           new_rna_struct;
    int chain_flag;
    float evaluate_index,
          best_evaluate_index = 0,
          antihairpin_index,
          lost_affinity_index;

    hairpin_nucleotide_links* new_HNL;
    new_HNL = new hairpin_nucleotide_links [(int)(rna_length / 2) + 1];
    int new_number_of_links;

    ofstream input_w;
    ifstream output_r;
    //*bigstat_w << "variants" << '\t';
    int max_dist_last_HNL_element_number = i;
    int mut_pstns [(*defined_rna).length()];
    int mut_nmb;
    int comb_nmb;
    int* n_pow_k_arr;
    int mut_pstn;

    for (int j = 0; j < 4; j++)
    {
        chain_flag = 1;
        mut_nmb = 0;
        for (int k = i - j; k >= max_dist_HNL_element_number; k -= 2)
        {
            if (chain_flag == 1)
                mut_pstns[mut_nmb] = HNL[k].first_nucleotide_number_in_link  - 1;
            else
                mut_pstns[mut_nmb] = HNL[k].second_nucleotide_number_in_link  - 1;

            mut_nmb++;
            chain_flag *= -1;
        }
        if (j >= 2)
        {
            mut_pstns[mut_nmb] = HNL[i - (j%2)].second_nucleotide_number_in_link - 1;
            mut_nmb++;
        }
        comb_nmb  = (int)pow(3, mut_nmb);
        n_pow_k_arr = new int [mut_nmb];

        for (int p = 0; p < mut_nmb; p++)
            n_pow_k_arr[p] = 0;

        for (int p = 0; p < comb_nmb; p++)
        {
            new_defined_rna = (*defined_rna);

            for (int m = 0; m < mut_nmb; m++)
            {
                mut_pstn = mut_pstns[m];
                new_defined_rna[mut_pstn] = except_defined_ncltd (new_defined_rna[mut_pstn],
                                                                  n_pow_k_arr[m],
                                                                  dna_flag                  );
            }

            complex_seq = (*input_rna) + "+" + new_defined_rna;

            antihairpin_index   = (dG_reader(&new_defined_rna) - initial_dG) / abs(initial_dG);
            lost_affinity_index = (dG_reader(&complex_seq) - initial_cmplx_dG) / abs( initial_cmplx_dG );
            evaluate_index = abs (antihairpin_index/lost_affinity_index);

            if (evaluate_index > best_evaluate_index)
            {
                best_evaluate_index = evaluate_index;
                best_defined_rna = new_defined_rna;
            }

            n_pow_k_arr[0]++;

            for (int l = 0; l < (mut_nmb - 1); l++)
            {
                if ( left_decade_pop( &(n_pow_k_arr[l + 1]),
                                      &(n_pow_k_arr[l]),
                                      2                    ) )
                    break;
            }
        }
        delete []n_pow_k_arr;
    }
    //cout << "ENTER any numb to continue" << endl;
    //cin >> pause;
    //delete [] new_HNL;

    complex_seq = (*input_rna) + "+" + best_defined_rna;
    float new_cmplx_dG = dG_reader(&complex_seq),
          new_dG = dG_and_struct_reader(&best_defined_rna,
                                        &new_rna_struct,
                                        new_HNL,
                                        &new_number_of_links);

    *stat_w << endl << endl;
    *stat_w << best_defined_rna << '\t'
            << new_dG << '\t'
            << new_cmplx_dG << '\t'
            << best_evaluate_index << endl;
    *stat_w << new_rna_struct << endl;

    if (end_program_check(new_dG,
                          target_hairpin_dG,
                          target_dG_accuracy,
                          new_cmplx_dG,
                          target_cmplx_dG,
                          target_cmplx_dG_accuracy,
                          &best_defined_rna,
                          defined_rna,
                          output_w,
                          input_rna,
                          new_number_of_links,
                          *number_of_links,
                          new_HNL,
                          HNL,
                          dna_flag,
                          max_dist_last_HNL_element_number,
                          max_dist_HNL_element_number,
                          initial_dG)           )
        return 1;
    if (break_flag)
    {
        cout << "MAX LOOP number is reached" << endl;
        (*output_w) << best_defined_rna << '\t'
                    << new_cmplx_dG << '\t'
                    << new_dG  << '\t'
                    << aff_reader(&best_defined_rna, input_rna, dna_flag)  << '\t'
                    << '\t' << '-' << endl;
        return 1;
    }
    *number_of_links = new_number_of_links;

    for (int l = 0; l < *number_of_links; l++)
    {
        HNL[l] = new_HNL[l];
    }

    (*defined_rna) = best_defined_rna;
    delete [] new_HNL;
    return 0;
}

int main_work  (string*   input_rna,
                string*   initial_rna,
                float*    target_hairpin_dG,
                float*    target_dG_accuracy,
                float*    target_cmplx_dG,
                float*    target_cmplx_dG_accuracy,
                ofstream* output_w,
                int       dna_flag            )
{
    rna_length = (*initial_rna).length();

    string defined_rna = (*initial_rna);

    hairpin_nucleotide_links* HNL;
    HNL = new hairpin_nucleotide_links [(int)(rna_length / 2) + 1];
    int number_of_links;

    string defined_rna_struct;
    float initial_dG = dG_and_struct_reader(&defined_rna,
                                            &defined_rna_struct,
                                            HNL,
                                            &number_of_links);

    string cmplx = (*input_rna) + "+" + defined_rna;
    float initial_cmplx_dG = dG_reader(&cmplx);

    int end_prog_flag = 0;
    char plus_or_minus;
    float new_cmplx_dG;

    if ( (*target_cmplx_dG) - initial_cmplx_dG > (*target_cmplx_dG_accuracy) &&
         ( abs ( initial_dG - (*target_hairpin_dG) ) < (*target_dG_accuracy) ||
           initial_dG >= (*target_hairpin_dG)                                   ) )
    {
        new_cmplx_dG = cmplx_dG_reducer (*target_cmplx_dG,
                                         *target_cmplx_dG_accuracy,
                                         &defined_rna,
                                         input_rna);

        if ( abs(new_cmplx_dG - (*target_cmplx_dG)) < (*target_cmplx_dG_accuracy) )

             plus_or_minus = '+';
        else plus_or_minus = '-';

        result_writer(output_w,
                      &defined_rna,
                      new_cmplx_dG,
                      dG_reader(&defined_rna),
                      aff_reader(&defined_rna, input_rna, dna_flag),
                      plus_or_minus);
        delete [] HNL;
        return 0;
    }

    if ( initial_cmplx_dG - (*target_cmplx_dG) >  (*target_cmplx_dG_accuracy) )
    {
        end_prog_flag = 1;
        plus_or_minus = '-';
    }
    else
    {
        if ( abs(initial_dG  - (*target_hairpin_dG)) <= (*target_dG_accuracy ) ||
            initial_dG >= (*target_hairpin_dG)                                    )
        {
            end_prog_flag = 1;
            plus_or_minus = '+';
        }
    }
    if (end_prog_flag)
    {
        result_writer(output_w,
                      &defined_rna,
                      initial_cmplx_dG,
                      initial_dG,
                      aff_reader(&defined_rna, input_rna, dna_flag),
                      plus_or_minus);
        delete [] HNL;
        return 0;
    }

    ofstream stat_w;
    stat_w.open(statistic);
    stat_w << (*input_rna) << '\t'
           << (*target_hairpin_dG) << '\t'
           << (*target_dG_accuracy) << '\t'
           << (*target_cmplx_dG) << '\t'
           << (*target_cmplx_dG_accuracy) << endl;
    stat_w << "structure"      << '\t'
      << "hairpin_dG"          << '\t'
      << "cmplx_dG"            << '\t'
      << "evaluate_index" << endl;

    stat_w << defined_rna    << '\t'
      << initial_dG          << '\t'
      << initial_cmplx_dG         << '\t'
      << 0 << endl;
    stat_w << defined_rna_struct  << endl;

    int break_flag = 0;
    //ofstream bigstat_w;
    //bigstat_w.open("./output/bigstat");

    for ( int i = 0; i < max_loop_size; i++ )
    {
        if (anti_hairpin(&defined_rna,
                         input_rna,
                         initial_dG,
                         initial_cmplx_dG,
                         HNL,
                         &number_of_links,
                         *target_hairpin_dG,
                         *target_dG_accuracy,
                         *target_cmplx_dG,
                         *target_cmplx_dG_accuracy,
                         output_w,
                         &stat_w,
                         break_flag,
                         dna_flag)       )
                         //&bigstat_w,
                         //i)             )
             break;
        if (i == max_loop_size - 2)
            break_flag = 1;
    }
    stat_w.close();
    delete [] HNL;
    return 0;
}


int string_before_separator_reader (string* parent_s,
                                    string* s,
                                    char sep,
                                    int start_position)
{
    int i;
    for (i = start_position; (*parent_s)[i] != sep && i < (*parent_s).length(); i++)
        (*s).push_back( (*parent_s)[i] );
    return i + 1;
}

void input_rna_reader(string*   input_str,
                      string*   input_rna,
                      float*    target_hairpin_dG,
                      float*    target_dG_accuracy,
                      float*    target_cmplx_dG,
                      float*    target_cmplx_dG_accuracy,
                      int*      dna_flag)
{
    char sep = ' ';
    int i = string_before_separator_reader (input_str, input_rna, sep, 0);

    string temp_string;
    i = string_before_separator_reader (input_str, &temp_string, sep, i);
    decimal_numb_reader( &temp_string, 0, target_hairpin_dG);
    temp_string.clear();

    i = string_before_separator_reader (input_str, &temp_string, sep, i);
    decimal_numb_reader( &temp_string, 0, target_dG_accuracy);
    temp_string.clear();

    i = string_before_separator_reader (input_str, &temp_string, sep, i);
    decimal_numb_reader( &temp_string, 0, target_cmplx_dG);
    temp_string.clear();

    i = string_before_separator_reader (input_str, &temp_string, sep, i);
    decimal_numb_reader( &temp_string, 0, target_cmplx_dG_accuracy);
    temp_string.clear();

    float temp_var;
    string_before_separator_reader (input_str, &temp_string, sep, i);
    decimal_numb_reader( &temp_string, 0, &temp_var);
    (*dna_flag) = (int) temp_var;
}

void initial_rna_maker(string* input_rna, string* c_input_rna, int dna_flag)
{
    c_seq_maker( input_rna, c_input_rna, dna_flag );
    string_inverser( c_input_rna );
}


int main()
{
    float target_hairpin_dG,
          target_dG_accuracy,
          target_cmplx_dG,
          target_cmplx_dG_accuracy;

    string input_rna,
           initial_rna,
           input_str;

    ifstream r;
    r.open(input_file);

    ofstream w;
    w.open(output_file);
    w    << "final_olig"   << '\t'
         << "cmplx_dG, kcal/mol"  << '\t'
         << "hairpin_dG, kcal/mol" << '\t'
         << "affinity, %"  << endl;

    int pause;
    int dna_flag = 0;
    while ( getline(r, input_str) )
    {
        input_rna_reader(&input_str,
                         &input_rna,
                         &target_hairpin_dG,
                         &target_dG_accuracy,
                         &target_cmplx_dG,
                         &target_cmplx_dG_accuracy,
                         &dna_flag);
        initial_rna_maker(&input_rna, &initial_rna, dna_flag);
        //initial_rna = "TTTTTGATTTGTGGCCGTTGTAGG";
        main_work(&input_rna,
                  &initial_rna,
                  &target_hairpin_dG,
                  &target_dG_accuracy,
                  &target_cmplx_dG,
                  &target_cmplx_dG_accuracy,
                  &w,
                  dna_flag);
        input_rna.clear();
        initial_rna.clear();
        /*
        cout << endl << "Enter any number to continue";
        cin >> pause;
        */
    }
    r.close();
    w.close();
    //system(CLEAN_FUNC);
return 0;
}
