#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <cmath>

using namespace std;

const char* input_file    = "./input/input_b"; // ./input directory and text file input should exist
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

const char* mfe_variants[4] = {"./input/var_1",
                               "./input/var_2",
                               "./input/var_3",
                               "./input/var_4"};

const int   max_loop_size = 15;
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

void cmplxes_dG_reader (string* AOB,
                        string* seq_a,
                        string* seq_b,
                        int seq_A_length,
                        int seq_O_length,
                        float* Aa_dG,
                        float* Bb_dG)
{
    string complex_seq;

    string seq_A = (*AOB);
    seq_A.erase(seq_A.begin() + seq_A_length, seq_A.end());
    complex_seq = (*seq_a) + "+" + seq_A;
    (*Aa_dG) = dG_reader(&complex_seq);

    string seq_B = (*AOB);
    seq_B.erase(seq_B.begin(), seq_B.end() - seq_A_length - seq_O_length + 1);
    complex_seq = (*seq_b) + "+" + seq_B;
    (*Bb_dG) = dG_reader(&complex_seq);
}

int careful_anti_hairpin(int     seq_A_length,
                         int     seq_O_length,
                         string* seq_a,
                         string* seq_b,
                         float   target_hairpin_dG,
                         float   target_dG_accuracy,
                         float     min_Aa_dG,
                         float     min_Aa_dG_accuracy,
                         float     min_Bb_dG,
                         float     min_Bb_dG_accuracy,
                         string* old_defined_rna,
                         string* old_undefined_rna,
                         int     number_of_links,
                         hairpin_nucleotide_links* HNL,
                         int     dna_flag,
                         int     max_dist_last_HNL_element_number,
                         int     max_dist_HNL_element_number,
                         float   initial_AOB_dG)
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
    float dG;
    float best_dG;
    float old_best_dG = initial_AOB_dG;
    float Aa_dG, Bb_dG;
    int best_dG_update_flag;

    //ofstream output_w;
    //output_w.open("./output/careful_antihairpin_stat");

    for(int mutant_ncltd_numb = 1; mutant_ncltd_numb <= ncltds_nmb_in_links_gr; mutant_ncltd_numb++)
    {
        best_dG_update_flag = 0;
        best_dG = initial_AOB_dG;
        mut_variants_nmb = C_n_k( ncltds_nmb_in_links_gr, mutant_ncltd_numb );
        arr[0] = new int* [mut_variants_nmb];
        for (int i = 0; i < mut_variants_nmb; i++)
            arr[0][i] = new int [mutant_ncltd_numb];

//cout << "n = " << ncltds_nmb_in_links_gr << endl
     //<< "k = " << mutant_ncltd_numb << endl
     //<< mut_variants_nmb << endl;

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

            cmplxes_dG_reader(&new_def_rna,
                              seq_a,
                              seq_b,
                              seq_A_length,
                              seq_O_length,
                              &Aa_dG,
                              &Bb_dG);

            if ( Aa_dG - min_Aa_dG > min_Aa_dG_accuracy ||
                 Bb_dG - min_Bb_dG > min_Bb_dG_accuracy   )
                continue;

            dG = dG_reader(&new_def_rna);
            //output_w << '\t' << new_undef_rna << '\t' << cmplx_dG << '\t' << dG << endl;

            if (dG > best_dG)
            {
                //best_evaluate_index = evaluate_index;
                best_def_rna = new_def_rna;
                best_dG_update_flag = 1;
                best_dG = dG;
                if ( abs (best_dG - target_hairpin_dG) < target_dG_accuracy )
                {
                    break_flag = 1;
                    break;
                }
            }
        }
        //output_w << endl << best_def_rna << '\t' << best_cmplx_dG << '\t' << best_dG << endl;

        if ( best_dG < old_best_dG    || //previous set of mutations gave better result than current
             best_dG_update_flag == 0 || //there are no mutations in current set that save required cmplx_dG level
             break_flag == 1            )//success
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
                    float     AOB_dG,
                    float     Aa_dG,
                    float     Bb_dG,
                    char      plus_or_minus)
{
    (*output_w) << (*olig_seq) << '\t'
                <<     AOB_dG     << '\t'
                <<     Aa_dG      << '\t'
                 <<     Bb_dG      << '\t'
                <<plus_or_minus<< endl;
}

int end_program_check (int     seq_A_length,
                       int     seq_O_length,
                       string* seq_a,
                       string* seq_b,
                       float     dG,
                       float     target_hairpin_dG,
                       float     target_dG_accuracy,
                       float     Aa_dG,
                       float     min_Aa_dG,
                       float     min_Aa_dG_accuracy,
                       float     Bb_dG,
                       float     min_Bb_dG,
                       float     min_Bb_dG_accuracy,
                       string*   defined_rna,
                       string*   old_defined_rna,
                       ofstream* output_w,
                       int       number_of_links,
                       hairpin_nucleotide_links*  HNL,
                       string*   old_undefined_rna,
                       int       dna_flag,
                       int       max_dist_last_HNL_element_number,
                       int       max_dist_HNL_element_number,
                       float     initial_AOB_dG)
{
    float new_cmplx_dG,
          new_dG,
          new_Aa_dG,
          new_Bb_dG;
    char  plus_or_minus;

    if ( Aa_dG - min_Aa_dG > min_Aa_dG_accuracy ||
         Bb_dG - min_Bb_dG > min_Bb_dG_accuracy   )
    {
        careful_anti_hairpin(seq_A_length,
                             seq_O_length,
                             seq_a,
                             seq_b,
                             target_hairpin_dG,
                             target_dG_accuracy,
                             min_Aa_dG,
                             min_Aa_dG_accuracy,
                             min_Bb_dG,
                             min_Bb_dG_accuracy,
                             old_defined_rna,
                             old_undefined_rna,
                             number_of_links,
                             HNL,
                             dna_flag,
                             max_dist_last_HNL_element_number,
                             max_dist_HNL_element_number,
                             initial_AOB_dG);

        cmplxes_dG_reader(old_defined_rna,
                              seq_a,
                              seq_b,
                              seq_A_length,
                              seq_O_length,
                              &new_Aa_dG,
                              &new_Bb_dG);
        new_dG  = dG_reader(old_defined_rna);

        if ( abs(new_Aa_dG - min_Aa_dG) < min_Aa_dG_accuracy &&
             abs(new_Bb_dG - min_Bb_dG) < min_Bb_dG_accuracy &&
             ( abs(new_dG  - target_hairpin_dG) <= target_dG_accuracy ||
               new_dG > target_hairpin_dG                               )  )

             plus_or_minus = '+';
        else plus_or_minus = '-';

        result_writer (output_w,
                       old_defined_rna,
                       new_dG,
                       new_Aa_dG,
                       new_Bb_dG,
                       plus_or_minus);
        return 1;
    }
    if ( abs(dG - target_hairpin_dG) <= target_dG_accuracy
                                                          ||
         dG > target_hairpin_dG                             )
    {
        result_writer (output_w,
                       defined_rna,
                       new_dG,
                       new_Aa_dG,
                       new_Bb_dG,
                       '+');
        return 1;
    }
    return 0;
}

int anti_hairpin  (int     seq_A_length,
                   int     seq_O_length,
                   string* undefined_rna,
                   string* defined_rna,
                   string* seq_a,
                   string* seq_b,
                   float   initial_AOB_dG,
                   float   initial_Aa_dG,
                   float   initial_Bb_dG,
                   hairpin_nucleotide_links*  HNL,
                   int*    number_of_links,
                   string* target_struct,
                   float*  antihairpin_index_,
                   float*  lost_affinity_index_,
                   float   target_hairpin_dG,
                   float   target_dG_accuracy,
                   float   min_Aa_dG,
                   float   min_Aa_dG_acrcy,
                   float   min_Bb_dG,
                   float   min_Bb_dG_acrcy,
                   ofstream* output_w,
                   ofstream* stat_w,
                   float   break_flag,
                   int     dna_flag)
                   //ofstream* bigstat_w,
                   //int loop_counter)
{
    int max_dist_HNL_element_number = 0;
    int lack_of_link_groups = *number_of_links; //*bigstat_w << "number_of_links" << '\t' << *number_of_links << endl;
    int i = 0;
    float Aa_dG_, Bb_dG_, dG_;

    while (i - max_dist_HNL_element_number == 0) //number of sequence links (one after another)
    {
        if (lack_of_link_groups <= 1)
        {
            cmplxes_dG_reader(defined_rna,
                              seq_a,
                              seq_b,
                              seq_A_length,
                              seq_O_length,
                              &Aa_dG_,
                              &Bb_dG_);

            dG_ = dG_reader(defined_rna);

            (*output_w) << (*defined_rna) << '\t'
                        << dG_ << '\t'
                        << Aa_dG_  << '\t'
                        << Bb_dG_  << '\t';

            if ( (abs( Aa_dG_ - min_Aa_dG ) < min_Aa_dG_acrcy     ||   // HERE CORRECTIONS NEED
                  Aa_dG_ < min_Aa_dG)                                 &&
                 (abs( Bb_dG_ - min_Bb_dG ) < min_Bb_dG_acrcy     ||
                  Bb_dG_ < min_Bb_dG)                                 &&
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
    string new_undefined_rna,
           new_defined_rna,
           new_rna_struct;
    int chain_flag;
    float mutation_comparison_table[4],
          antihairpin_index[4],
          lost_affinity_index[4],
          new_dG,
          new_Aa_dG,
          new_Bb_dG;
    hairpin_nucleotide_links* new_HNL;
    new_HNL = new hairpin_nucleotide_links [(int)(rna_length / 2) + 1];
    int new_number_of_links;

    ofstream input_w;
    ifstream output_r;
    //*bigstat_w << "variants" << '\t';
    int max_dist_last_HNL_element_number = i;

    for (int j = 0; j < 4; j++)
    {
        new_undefined_rna = *undefined_rna;
        chain_flag = 1;
        for (int k = i - j; k >= max_dist_HNL_element_number; k -= 2)
        {
            if (chain_flag == 1)
            {
                if (new_undefined_rna[ HNL[k].first_nucleotide_number_in_link  - 1] != 'N')
                    new_undefined_rna[ HNL[k].first_nucleotide_number_in_link  - 1] =  'N';
                else
                    new_undefined_rna[ HNL[k].second_nucleotide_number_in_link - 1] =  'N';
            } else
            {
                if (new_undefined_rna[ HNL[k].second_nucleotide_number_in_link - 1] != 'N')
                    new_undefined_rna[ HNL[k].second_nucleotide_number_in_link - 1] =  'N';
                else
                    new_undefined_rna[ HNL[k].first_nucleotide_number_in_link  - 1] =  'N';
            }
            chain_flag *= -1;
        }
        if (j >= 2)
        new_undefined_rna[ HNL[i - (j%2)].second_nucleotide_number_in_link - 1] = 'N';

        design_maker (target_struct,
                      &new_undefined_rna,
                      &new_defined_rna,
                      dna_flag   );
        //*bigstat_w << new_undefined_rna << '\t';
        //*bigstat_w << new_defined_rna << '\t';

        new_dG = dG_and_struct_reader(&new_defined_rna,
                                      &new_rna_struct,
                                       new_HNL,
                                      &new_number_of_links);//*bigstat_w << new_dG << '\t';

        cmplxes_dG_reader(&new_defined_rna,
                          seq_a,
                          seq_b,
                          seq_A_length,
                          seq_O_length,
                          &new_Aa_dG,
                          &new_Bb_dG);

        antihairpin_index[j]   = (new_dG - initial_AOB_dG) / abs(initial_AOB_dG); //*bigstat_w << antihairpin_index[j] << '\t';// higher is better
        lost_affinity_index[j] = (new_Aa_dG - initial_Aa_dG) / abs( initial_Aa_dG ) +
                                 (new_Bb_dG - initial_Bb_dG) / abs( initial_Bb_dG );   //*bigstat_w << lost_affinity_index[j] << '\t';// lower is better
        //if (lost_affinity_index[j] < abs(antihairpin_index[j]/10))
            //mutation_comparison_table[j] = abs (antihairpin_index[j]);
        //else
            mutation_comparison_table[j] =                     // how many percents of hairpin destruction
            abs (antihairpin_index[j]/lost_affinity_index[j]); // we have on one percent of affinity loss
        input_w.open(mfe_variants[j]); //*bigstat_w << mutation_comparison_table[j] << endl << '\t';
        input_w << new_undefined_rna << endl <<
             new_defined_rna   << endl <<
             new_rna_struct << endl <<
             new_dG            << endl <<
             new_Aa_dG            << endl <<
             new_Bb_dG            << endl <<
             new_number_of_links;

        for (int l = 0; l < new_number_of_links; l++)
            input_w << endl << new_HNL[l].first_nucleotide_number_in_link << '\t' <<
                         new_HNL[l].second_nucleotide_number_in_link;
        input_w.close();
    }
    //cout << "ENTER any numb to continue" << endl;
    //cin >> pause;
    //delete [] new_HNL;
    i = max_arr_elem_numb( mutation_comparison_table, 4 ); //*bigstat_w << endl << "i" << '\t' <<  i << endl;

    *antihairpin_index_   = antihairpin_index[i];
    *lost_affinity_index_ = lost_affinity_index[i];

    output_r.open(mfe_variants[i]);
    getline(output_r, new_undefined_rna);
    getline(output_r, new_defined_rna);
    getline(output_r, new_rna_struct);

    string temp_str;
    getline(output_r, temp_str);
    decimal_numb_reader(&temp_str, 0, &new_dG);

    getline(output_r, temp_str);
    decimal_numb_reader(&temp_str, 0, &new_Aa_dG);

    getline(output_r, temp_str);
    decimal_numb_reader(&temp_str, 0, &new_Bb_dG);

    *stat_w << new_undefined_rna << endl << endl;
    *stat_w << new_defined_rna          << '\t'
          << new_dG                 << '\t'
          << new_Aa_dG         << '\t'
          << new_Bb_dG         << '\t'
          << *antihairpin_index_   << '\t'
          << *lost_affinity_index_ << endl;
    *stat_w << new_rna_struct << endl;

    float temp_f;
    getline(output_r, temp_str);
    decimal_numb_reader(&temp_str, 0, &temp_f);
    new_number_of_links = (int)temp_f;

    for (int l = 0; l < new_number_of_links; l++)
    {
        getline(output_r, temp_str);
        links_reader(&new_HNL[l], &temp_str);
    }

    int break_f = 1;
    for (int f = 0; f < new_undefined_rna.length(); f++)
    {
        if ( new_undefined_rna[f] != (*undefined_rna)[f] )
        {
            break_f = 0;
            break;
        }
    }
    if (break_f)
    {
        for (int k = max_dist_last_HNL_element_number; k >= max_dist_HNL_element_number; k--)
        {
            new_undefined_rna[ HNL[k].first_nucleotide_number_in_link  - 1] =  'N';
            new_undefined_rna[ HNL[k].second_nucleotide_number_in_link  - 1] =  'N';
        }
    }

    if (end_program_check(seq_A_length,
                          seq_O_length,
                          seq_a,
                          seq_b,
                          new_dG,
                          target_hairpin_dG,
                          target_dG_accuracy,
                          new_Aa_dG,
                          min_Aa_dG,
                          min_Aa_dG_acrcy,
                          new_Bb_dG,
                          min_Bb_dG,
                          min_Bb_dG_acrcy,
                          &new_defined_rna,
                          defined_rna,
                          output_w,
                          *number_of_links,
                          HNL,
                          undefined_rna,
                          dna_flag,
                          max_dist_last_HNL_element_number,
                          max_dist_HNL_element_number,
                          initial_AOB_dG)           )
        return 1;
    if (break_flag)
    {
        (*output_w) << new_defined_rna << '\t'
                    << new_dG  << '\t'
                    << new_Aa_dG << '\t'
                    << new_Bb_dG << '\t'
                    << '\t' << '-' << endl;
        return 1;
    }
    *number_of_links = new_number_of_links;

    for (int l = 0; l < *number_of_links; l++)
    {
        HNL[l] = new_HNL[l];
    }

    (*undefined_rna) = new_undefined_rna;
    (*defined_rna) = new_defined_rna;
    delete [] new_HNL;
    return 0;
}

int main_work  (  string* seq_a,
                  string* seq_b,
                  string* seq_A,
                  string* seq_O,
                  string* seq_B,
                  float  min_Aa_dG,
                  float  min_Aa_dG_acrcy,
                  float  min_Bb_dG,
                  float  min_Bb_dG_acrcy,
                  float  target_AOB_dG,
                  float  target_AOB_dG_acrcy,
                  int    dna_flag,
                  ofstream* output_w)
{
    string defined_rna = (*seq_A) + (*seq_O) + (*seq_B),
           undefined_rna = defined_rna;

    rna_length = defined_rna.length();

    string target_struct;
    for (int i = 0; i <  rna_length; i++)
        target_struct.push_back('.');

    hairpin_nucleotide_links* HNL;
    HNL = new hairpin_nucleotide_links [(int)(rna_length / 2) + 1];
    int number_of_links;
//defined_rna = "TCGTGTTATAAGTTAGTATTGGTTGGACGCTATTCAGGAAGATTGTTGGTTAAGGG";
    string defined_rna_struct;
    float initial_AOB_dG = dG_and_struct_reader(&defined_rna,
                                            &defined_rna_struct,
                                            HNL,
                                            &number_of_links);

    string cmplx = (*seq_A) + "+" + (*seq_a);
    float initial_Aa_dG = dG_reader(&cmplx);

    cmplx = (*seq_B) + "+" + (*seq_b);
    float initial_Bb_dG = dG_reader(&cmplx);

    float antihairpin_index = 0,
          lost_affinity_index = 0;

    int end_prog_flag = 0;
    char plus_or_minus;
    float new_cmplx_dG;

    if ( initial_Aa_dG - min_Aa_dG >  min_Aa_dG_acrcy ||
         initial_Bb_dG - min_Bb_dG >  min_Bb_dG_acrcy)
    {
        end_prog_flag = 1;
        plus_or_minus = '-';
    }
    else
    {
        if ( abs(initial_AOB_dG  - target_AOB_dG) <= target_AOB_dG_acrcy ||
            initial_AOB_dG >= target_AOB_dG                                    )
        {
            end_prog_flag = 1;
            plus_or_minus = '+';
        }
    }
    if (end_prog_flag)
    {
        result_writer(output_w,
                      &defined_rna,
                      initial_AOB_dG,
                      initial_Aa_dG,
                      initial_Bb_dG,
                      plus_or_minus);
        delete [] HNL;
        return 0;
    }
    ofstream stat_w;
    stat_w.open(statistic);
    stat_w << "structure"      << '\t'
      << "hairpin_dG"          << '\t'
      << "Aa_dG"            << '\t'
      << "Bb_dG"            << '\t'
      << "antihairpin_index"   << '\t'
      << "lost_affinity_index" << endl;

    stat_w << defined_rna    << '\t'
      << initial_AOB_dG          << '\t'
      << initial_Aa_dG         << '\t'
      << initial_Bb_dG         << '\t'
      << antihairpin_index   << '\t'
      << lost_affinity_index << endl;
    stat_w << defined_rna_struct  << endl;

    int break_flag = 0;
    //ofstream bigstat_w;
    //bigstat_w.open("./output/bigstat");

//undefined_rna = "TCGTGTTNTNNGNTANNATTGGTTGGNCGCTATTCAGGAAGATNGTTGGNNAANNG";
    for ( int i = 0; i < max_loop_size; i++ )
    {
        if (anti_hairpin((*seq_A).length(),
                         (*seq_O).length(),
                         &undefined_rna,
                         &defined_rna,
                         seq_a,
                         seq_b,
                         initial_AOB_dG,
                         initial_Aa_dG,
                         initial_Bb_dG,
                         HNL,
                         &number_of_links,
                         &target_struct,
                         &antihairpin_index,
                         &lost_affinity_index,
                         target_AOB_dG,
                         target_AOB_dG_acrcy,
                         min_Aa_dG,
                         min_Aa_dG_acrcy,
                         min_Bb_dG,
                         min_Bb_dG_acrcy,
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

void seq_after_sep_reader (ifstream* inp_r,
                           string* seq,
                           char sep)
{
    string temp_str;
    getline((*inp_r), temp_str);
    temp_str.erase(temp_str.begin(),
                   temp_str.begin() + temp_str.find(sep) + 1);
    (*seq) = temp_str;
}

void numb_after_sep_reader (ifstream* inp_r,
                            float* numb,
                            char sep)
{
    string temp_str;
    getline((*inp_r), temp_str);
    temp_str.erase(temp_str.begin(),
                   temp_str.begin() + temp_str.find(sep) + 1);
    decimal_numb_reader(&temp_str, 0, numb);
}

void help_reader (ifstream* inp_r,
                  string* seq,
                  float* min_dG,
                  float* min_dG_acc,
                  char sep)
{
    seq_after_sep_reader (inp_r, seq, sep);

    numb_after_sep_reader (inp_r, min_dG, sep);
    numb_after_sep_reader (inp_r, min_dG_acc, sep);
}


void input_reader(string* seq_a,
                  string* seq_b,
                  string* seq_A,
                  string* seq_O,
                  string* seq_B,
                  float*  min_Aa_dG,
                  float*  min_Aa_dG_acrcy,
                  float*  min_Bb_dG,
                  float*  min_Bb_dG_acrcy,
                  float*  target_AOB_dG,
                  float*  target_AOB_dG_acrcy,
                  int*    dna_flag)
{
    char sep = ' ';
    ifstream inp_r;
    inp_r.open(input_file);

    string temp_str;

    help_reader (&inp_r, seq_a, min_Aa_dG, min_Aa_dG_acrcy, sep);

    getline(inp_r, temp_str);

    help_reader (&inp_r, seq_b, min_Bb_dG, min_Bb_dG_acrcy, sep);

    getline(inp_r, temp_str);

    seq_after_sep_reader (&inp_r, seq_A, sep);
    seq_after_sep_reader (&inp_r, seq_O, sep);
    seq_after_sep_reader (&inp_r, seq_B, sep);


    getline(inp_r, temp_str);

    numb_after_sep_reader (&inp_r, target_AOB_dG, sep);
    numb_after_sep_reader (&inp_r, target_AOB_dG_acrcy, sep);

    getline(inp_r, temp_str);

    float dna_flag_f;
    numb_after_sep_reader (&inp_r, &dna_flag_f, sep);
    (*dna_flag) = (int)dna_flag_f;

    inp_r.close();
}

int main()
{
    int dna_flag;

    float min_Aa_dG,
          min_Aa_dG_acrcy,
          min_Bb_dG,
          min_Bb_dG_acrcy,
          target_AOB_dG,
          target_AOB_dG_acrcy;

    string seq_a,
           seq_b,
           seq_A,
           seq_O,
           seq_B;



    ofstream w;
    w.open(output_file);
    w    << "seq_AOB"   << '\t'
         << "AOB_dG, kcal/mol"  << '\t'
         << "Aa_dG, kcal/mol"  << '\t'
         << "Bb_dG, kcal/mol" << endl;

    input_reader( &seq_a,
                  &seq_b,
                  &seq_A,
                  &seq_O,
                  &seq_B,
                  &min_Aa_dG,
                  &min_Aa_dG_acrcy,
                  &min_Bb_dG,
                  &min_Bb_dG_acrcy,
                  &target_AOB_dG,
                  &target_AOB_dG_acrcy,
                  &dna_flag);

    main_work(    &seq_a,
                  &seq_b,
                  &seq_A,
                  &seq_O,
                  &seq_B,
                  min_Aa_dG,
                  min_Aa_dG_acrcy,
                  min_Bb_dG,
                  min_Bb_dG_acrcy,
                  target_AOB_dG,
                  target_AOB_dG_acrcy,
                  dna_flag,
                  &w);

        /*
        cout << endl << "Enter any number to continue";
        cin >> pause;
        */
    w.close();
    //system(CLEAN_FUNC);
return 0;
}
