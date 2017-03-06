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

const char* mfe_variants[4] = {"./input/var_1",
                               "./input/var_2",
                               "./input/var_3",
                               "./input/var_4"};

const int   max_loop_size = 10;
const int   absolute_aff_changing = 1;
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

int c_seq_maker (string* seq, string* c_seq)
{
    int l = (*seq).length();
    for (int i = 0; i < l; i++ )
    {
        if ((*seq)[i] == 'A') (*c_seq).push_back('U');
        else if ((*seq)[i] == 'U') (*c_seq).push_back('A');
             else if ((*seq)[i] == 'G') (*c_seq).push_back('C');
                  else if ((*seq)[i] == 'C') (*c_seq).push_back('G');
                       else { cout << "error: incorrect sequence" << endl;
                              return 1;
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
    $NUPACKHOME/bin/complexes -T 23 -material rna $CODDING_HOME/input/complex \n\
    $NUPACKHOME/bin/concentrations $CODDING_HOME/input/complex \n\
"

#define COPY_F "\
    #/bin/bash \n\
    cp ./input/input_seq.fold ./output/input_seq.fold \n\
    cp ./input/input_seq.summary ./output/input_seq.summary \n\
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

int aff_reducer_1(string* input_rna,      //it's like simple total antihairpin (additional)
                  string* new_defined_rna,
                  string* undefined_rna,
                  string* target_struct,
                  int*    number_of_links,
                  float*  new_aff,
                  float   aff,
                  hairpin_nucleotide_links*  HNL)
{
    string last_variant = (*new_defined_rna);
    for (int i = 0; i < *number_of_links; i++)
    {
        if ( (*undefined_rna)[HNL[i].first_nucleotide_number_in_link - 1] == 'N' )
            (*undefined_rna)[HNL[i].second_nucleotide_number_in_link - 1] = 'N';
        else (*undefined_rna)[HNL[i].first_nucleotide_number_in_link - 1] = 'N';
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
    string temp_str;
    if (dG_and_struct_reader (new_defined_rna,
                              &temp_str,
                              HNL,
                              number_of_links)  < 0)
    {
        if (aff_reducer_1(input_rna,
                          new_defined_rna,
                          undefined_rna,
                          target_struct,
                          number_of_links,
                          new_aff,
                          aff,
                          HNL)              )
        return 1;
    }
    return 0; // aff didn't start change significantly
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

float aff_reducer(float   aff,
                  float   target_aff,
                  float   target_aff_accuracy,
                  string* defined_rna,
                  string* defined_rna_struct,
                  string* input_rna,
                  int     number_of_links,
                  hairpin_nucleotide_links*  HNL)
{
    string undefined_rna = (*defined_rna);
    string target_struct;
    string new_defined_rna = (*defined_rna);

    for (int i = 0; i < (*defined_rna).length(); i++)
        target_struct.push_back('.');

    float new_aff;
    string complex_struct;

    //ofstream test_w;
    //test_w.open("./input/last_olig");
    if (aff < 100)
        goto goto_label;

    if ( aff_reducer_1(input_rna,          // aff didn't start change significantly
                       &new_defined_rna,
                       &undefined_rna,
                       &target_struct,
                       &number_of_links,
                       &new_aff,
                       aff,
                       HNL)               == 0 )
    {
        undefined_rna = new_defined_rna;
        target_olig_structure_part_in_complex_reader (&new_defined_rna,
                                                      input_rna,
                                                      &complex_struct);
        aff_reducer_2(input_rna,
                      &complex_struct,
                      &new_defined_rna,
                      &undefined_rna,
                      &target_struct,
                      &new_aff,
                      new_aff          );
    }
    //test_w.close();

    //return new_aff;
    goto_label:
    target_olig_structure_part_in_complex_reader (&new_defined_rna,
                                                  input_rna,
                                                  &complex_struct);
    if ( abs(new_aff - target_aff) < target_aff_accuracy )
    {
        (*defined_rna) = new_defined_rna;
        return new_aff;
    }
    undefined_rna = new_defined_rna;
    careful_aff_reducer( &complex_struct,
                         &undefined_rna,
                         &new_defined_rna,
                         &target_struct,
                         input_rna,
                         target_aff,
                         target_aff_accuracy,
                         &new_aff);
    (*defined_rna) = new_defined_rna;
    return new_aff;
    //string rna_complex = (*input_rna) + "+" + (*defined_rna);
}

int end_program_check (float     dG,
                       float     target_hairpin_dG,
                       float     target_dG_accuracy,
                       float     aff,
                       float     target_aff,
                       float     target_aff_accuracy,
                       string*   defined_rna,
                       ofstream* output_w,
                       string*   defined_rna_struct,
                       string*   input_rna,
                       int       number_of_links,
                       hairpin_nucleotide_links*  HNL)
{
    float new_aff,
          new_dG;
    if ( abs(dG - target_hairpin_dG) < target_dG_accuracy
                                                          ||
         dG > target_hairpin_dG                             )
    {
        if ( aff - target_aff > target_aff_accuracy )
        {
            new_aff = aff_reducer(aff,              //here dG can decrease!
                                  target_aff,
                                  target_aff_accuracy,
                                  defined_rna,
                                  defined_rna_struct,
                                  input_rna,
                                  number_of_links,
                                  HNL                 );

            new_dG = dG_reader(defined_rna);
            (*output_w) << (*defined_rna) << '\t'
                        << new_aff << '\t'
                        << new_dG  << '\t';
            if ( abs(new_aff - target_aff) < target_aff_accuracy )
            {
                (*output_w) << '+' << endl;
            }
            else (*output_w) << '-' << endl;
        }
        else
        {
            (*output_w) << (*defined_rna) << '\t'
                        << aff << '\t'
                        << dG  << '\t';
            if ( abs(aff - target_aff) < target_aff_accuracy )
            {
                (*output_w) << '+' << endl;
            }
            else (*output_w) << '-' << endl;
        }
        return 1;
    }
    return 0;
}


int anti_hairpin  (string* undefined_rna,
                   string* defined_rna,
                   string* input_rna,
                   float   initial_dG,
                   float   initial_cmplx_dG,
                   float   initial_aff,
                   hairpin_nucleotide_links*  HNL,
                   int*    number_of_links,
                   string* target_struct,
                   float*  antihairpin_index_,
                   float*  lost_affinity_index_,
                   float   target_hairpin_dG,
                   float   target_dG_accuracy,
                   float   target_aff,
                   float   target_aff_accuracy,
                   ofstream* output_w,
                   ofstream* stat_w,
                   float break_flag)
                   //ofstream* bigstat_w,
                   //int loop_counter)
{
    int max_dist_HNL_element_number = 0;
    int lack_of_link_groups = *number_of_links; //*bigstat_w << "number_of_links" << '\t' << *number_of_links << endl;
    int i = 0;
    float aff_, dG_;
    while (i - max_dist_HNL_element_number == 0) //number of sequence links (one after another)
    {
        if (lack_of_link_groups <= 1)
        {
            aff_ = aff_reader(input_rna, defined_rna);
            dG_ = dG_and_struct_reader(defined_rna,
                                       undefined_rna, // JUST after return we won't need it
                                       HNL,
                                       number_of_links    );
            (*output_w) << (*defined_rna) << '\t'
                        << aff_ << '\t'
                        << dG_  << '\t';
        if ( (abs( aff_ - target_aff ) < target_dG_accuracy     ||   // HERE CORRECTIONS NEED
             aff_ > target_aff)                                 &&
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
          new_aff;
    hairpin_nucleotide_links* new_HNL;
    new_HNL = new hairpin_nucleotide_links [(int)(rna_length / 2) + 1];
    int new_number_of_links;

    ofstream input_w;
    ifstream output_r;
    //*bigstat_w << "variants" << '\t';
    string complex_seq;
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
                      &new_defined_rna   );
        //*bigstat_w << new_undefined_rna << '\t';
        //*bigstat_w << new_defined_rna << '\t';

        //if (loop_counter == 1 && j == 3)
            //system (COPY_F);

        new_dG = dG_and_struct_reader(&new_defined_rna,
                                      &new_rna_struct,
                                       new_HNL,
                                      &new_number_of_links);//*bigstat_w << new_dG << '\t';

        new_aff = aff_reader(input_rna, &new_defined_rna);
        complex_seq = (*input_rna) + "+" + new_defined_rna;

        antihairpin_index[j]   = (new_dG - initial_dG) / abs(initial_dG); //*bigstat_w << antihairpin_index[j] << '\t';// higher is better
        lost_affinity_index[j] = abs( (initial_cmplx_dG - dG_reader(&complex_seq)) / initial_cmplx_dG );   //*bigstat_w << lost_affinity_index[j] << '\t';// lower is better
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
             new_aff            << endl <<
             new_number_of_links;
        for (int l = 0; l < new_number_of_links; l++)
            input_w << endl << new_HNL[l].first_nucleotide_number_in_link << '\t' <<
                         new_HNL[l].second_nucleotide_number_in_link;
        input_w.close();
    }
    delete [] new_HNL;
    i = max_arr_elem_numb( mutation_comparison_table, 4 ); //*bigstat_w << endl << "i" << '\t' <<  i << endl;

    *antihairpin_index_   = antihairpin_index[i];
    *lost_affinity_index_ = lost_affinity_index[i];

    output_r.open(mfe_variants[i]);
    getline(output_r, new_undefined_rna);
    getline(output_r, *defined_rna);
    getline(output_r, new_rna_struct);

    string temp_str;
    getline(output_r, temp_str);
    decimal_numb_reader(&temp_str, 0, &new_dG);
    getline(output_r, temp_str);
    decimal_numb_reader(&temp_str, 0, &new_aff);

    *stat_w << new_undefined_rna << endl << endl;
    *stat_w << *defined_rna          << '\t'
          << new_dG                 << '\t'
          << new_aff         << '\t'
          << *antihairpin_index_   << '\t'
          << *lost_affinity_index_ << endl;
    *stat_w << new_rna_struct << endl;

    float temp_f;
    getline(output_r, temp_str);
    decimal_numb_reader(&temp_str, 0, &temp_f);
    *number_of_links = (int)temp_f;

    for (int l = 0; l < *number_of_links; l++)
    {
        getline(output_r, temp_str);
        links_reader(&HNL[l], &temp_str);
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
    (*undefined_rna) = new_undefined_rna;

    if (end_program_check(new_dG,
                          target_hairpin_dG,
                          target_dG_accuracy,
                          new_aff,
                          target_aff,
                          target_aff_accuracy,
                          defined_rna,
                          output_w,
                          &new_rna_struct,
                          input_rna,
                          *number_of_links,
                          HNL)           )
        return 1;
    if (break_flag)
    {
        (*output_w) << (*defined_rna) << '\t'
                    << new_aff << '\t'
                    << new_dG  << '\t'
                    << '\t' << '-' << endl;
        return 1;
    }
    return 0;
}

int main_work  (string*   input_rna,
                string*   initial_rna,
                float*    target_hairpin_dG,
                float*    target_dG_accuracy,
                float*    target_aff,
                float*    target_aff_accuracy,
                ofstream* output_w                         )
{
    rna_length = (*initial_rna).length();

    string defined_rna = (*initial_rna),
           undefined_rna = defined_rna;

    string target_struct;
    for (int i = 0; i <  rna_length; i++)
        target_struct.push_back('.');

    hairpin_nucleotide_links* HNL;
    HNL = new hairpin_nucleotide_links [(int)(rna_length / 2) + 1];
    int number_of_links;

    string defined_rna_struct;
    float initial_dG = dG_and_struct_reader(&defined_rna,
                                            &defined_rna_struct,
                                            HNL,
                                            &number_of_links);

    float initial_aff = aff_reader(input_rna, &defined_rna);
    float antihairpin_index = 0,
          lost_affinity_index = 0;

    if ( end_program_check (initial_dG,
                            *target_hairpin_dG,
                            *target_dG_accuracy,
                            initial_aff,
                            *target_aff,
                            *target_aff_accuracy,
                            &defined_rna,
                            output_w,
                            &defined_rna_struct,
                            input_rna,
                            number_of_links,
                            HNL)                 )
        return 0;
    ofstream stat_w;
    stat_w.open(statistic);
    stat_w << (*input_rna) << '\t'
           << (*target_hairpin_dG) << '\t'
           << (*target_dG_accuracy) << '\t'
           << (*target_aff) << '\t'
           << (*target_aff_accuracy) << endl;
    stat_w << "structure"      << '\t'
      << "hairpin_dG"          << '\t'
      << "affinity"            << '\t'
      << "antihairpin_index"   << '\t'
      << "lost_affinity_index" << endl;

    stat_w << defined_rna    << '\t'
      << initial_dG          << '\t'
      << initial_aff         << '\t'
      << antihairpin_index   << '\t'
      << lost_affinity_index << endl;
    stat_w << defined_rna_struct  << endl;

    int break_flag = 0;
    //ofstream bigstat_w;
    //bigstat_w.open("./output/bigstat");
    string complex_seq = (*input_rna) + "+" + defined_rna;
    float initial_cmplx_dG = dG_reader(&complex_seq);
    for ( int i = 0; i < max_loop_size; i++ )
    {
        if (anti_hairpin(&undefined_rna,
                         &defined_rna,
                         input_rna,
                         initial_dG,
                         initial_cmplx_dG,
                         initial_aff,
                         HNL,
                         &number_of_links,
                         &target_struct,
                         &antihairpin_index,
                         &lost_affinity_index,
                         *target_hairpin_dG,
                         *target_dG_accuracy,
                         *target_aff,
                         *target_aff_accuracy,
                         output_w,
                         &stat_w,
                         break_flag)       )
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
                      float*    target_aff,
                      float*    target_aff_accuracy)
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
    decimal_numb_reader( &temp_string, 0, target_aff);
    temp_string.clear();

    string_before_separator_reader (input_str, &temp_string, sep, i);
    decimal_numb_reader( &temp_string, 0, target_aff_accuracy);
}

void initial_rna_maker(string* input_rna, string* c_input_rna)
{
    c_seq_maker( input_rna, c_input_rna );
    string_inverser( c_input_rna );
}


int main()
{
    float target_hairpin_dG,
          target_dG_accuracy,
          target_aff,
          target_aff_accuracy;

    string input_rna,
           initial_rna,
           input_str;

    ifstream r;
    r.open(input_file);

    ofstream w;
    w.open(output_file);
    w    << "final_olig"   << '\t'
         << "affinity, %"  << '\t'
         << "dG, kcal/mol" << endl;

    while ( getline(r, input_str) )
    {
        input_rna_reader(&input_str,
                         &input_rna,
                         &target_hairpin_dG,
                         &target_dG_accuracy,
                         &target_aff,
                         &target_aff_accuracy);
        initial_rna_maker(&input_rna, &initial_rna);
        main_work(&input_rna,
                  &initial_rna,
                  &target_hairpin_dG,
                  &target_dG_accuracy,
                  &target_aff,
                  &target_aff_accuracy,
                  &w);
        input_rna.clear();
        initial_rna.clear();
    }
    r.close();
    w.close();
    return 0;
}
