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

long int fctrl (int n)
{
    if (n == 1) return 1;
    return n*fctrl(n-1);
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


int careful_aff_reducer(string* complex_struct,
                         string* undefined_rna,
                         string* new_defined_rna,
                         string* target_struct,
                         string* input_rna,
                         float   target_aff      )
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
    string new_def_rna = (*new_defined_rna);
    string best_def_rna;
    string old_best_def_rna;
    float d_aff;
    float min_d_aff = target_aff;
    float old_min_d_aff = target_aff;
    for(int mutant_ncltd_numb = 1; ; mutant_ncltd_numb++)
    {
        mut_variants_nmb = C_n_k(cmplx_nmb_of_links, mutant_ncltd_numb); //cout << cmplx_nmb_of_links;
        arr[0] = new int* [mut_variants_nmb];
        for (int i = 0; i < mut_variants_nmb; i++)
            arr[0][i] = new int [mutant_ncltd_numb];
//cout << cmplx_nmb_of_links << ' ' << mutant_ncltd_numb << ' ' << mut_variants_nmb;
//return 1;
        comb_maker (arr[0], cmplx_nmb_of_links, mutant_ncltd_numb, mut_variants_nmb);

        for (int i = 0; i < mut_variants_nmb; i++)
        {
            for (int j = 0; j < mutant_ncltd_numb; j++)
            {
                temp_var = ncltds_in_cmplx_nmbs[ arr[0][i][j] ];
                (*undefined_rna)[temp_var] = except_ncltd_symb ((*undefined_rna)[temp_var]);
            }
            design_maker(target_struct, undefined_rna, &new_def_rna);
            d_aff = abs( aff_reader(input_rna, &new_def_rna) - target_aff );

            if (d_aff < min_d_aff)
            {
                min_d_aff = d_aff;
                best_def_rna = new_def_rna;
            }
            for (int j = 0; j < mutant_ncltd_numb; j++)
            {
                temp_var = ncltds_in_cmplx_nmbs[ arr[0][i][j] ];
                (*undefined_rna)[temp_var] = (*new_defined_rna)[temp_var];
            }
        }
        if (min_d_aff > old_min_d_aff)
        {
            (*new_defined_rna) = old_best_def_rna;
            break;
        }
        old_best_def_rna = best_def_rna;
        old_min_d_aff = min_d_aff;
        min_d_aff = target_aff;
        for (int i = 0; i < mut_variants_nmb; i++)
            delete [] arr[0][i];
    }
    delete [] arr;
    return 1;
}


int main()
{   // DESIGN_FUNC MFE_FUNC CONCENTRATIONS_FUNC
    //system(MFE_FUNC);
    /*string input_rna =       "GCCGGCUCAAUGAGCGUGGGGAUCCCUGGCGUG";
    string new_defined_rna = "UCAUUUUCAUUCAUUUUCGUUCAUUAAGCCGGU";
    string target_struct =   ".................................";
    string complex_struct =  ".................)))))))).)))))).";
    string undefined_rna = new_defined_rna;
    careful_aff_reducer( &complex_struct,
                         &undefined_rna,
                         &new_defined_rna,
                         &target_struct,
                         &input_rna,
                         75.0               );
    */system(CONCENTRATIONS_FUNC);
    return 0;
}
