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

const int   max_loop_size = 15;

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

void result_writer (ofstream* output_w,
                    string*   olig_seq,
                    float     dG,
                    char      plus_or_minus)
{
    (*output_w) << (*olig_seq) << '\t'
                <<     dG      << '\t'
                <<plus_or_minus<< endl;
}

int anti_hairpin  (string* undefined_rna,
                   string* defined_rna,
                   hairpin_nucleotide_links*  HNL,
                   int*    number_of_links,
                   string* target_struct,
                   float   target_hairpin_dG,
                   float   target_dG_accuracy,
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
    float dG_;

    while (i - max_dist_HNL_element_number == 0) //number of sequence links (one after another)
    {
        if (lack_of_link_groups <= 1)
        {
            dG_ = dG_reader(defined_rna);

            (*output_w) << (*defined_rna) << '\t'
                        << dG_  << '\t';

            if ( abs(dG_ - target_hairpin_dG) < target_dG_accuracy ||
                 dG_ > target_hairpin_dG                              )
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
           new_rna_struct,
           old_undefined_rna = (*undefined_rna);

    int chain_flag;

    float new_dG,
          best_new_dG = -1000;

    ofstream input_w;
    ifstream output_r;
    //*bigstat_w << "variants" << '\t';
    int max_dist_last_HNL_element_number = i;

    for (int j = 0; j < 4; j++)
    {
        new_undefined_rna = old_undefined_rna;
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

        new_dG = dG_reader(&new_defined_rna);
        if (new_dG > best_new_dG)
        {
            best_new_dG = new_dG;
            (*defined_rna) = new_defined_rna;
            (*undefined_rna) = new_undefined_rna;
        }
        /*
        dG_and_struct_reader(&new_defined_rna,
                                      &new_rna_struct,
                                       new_HNL,
                                      &new_number_of_links);//*bigstat_w << new_dG << '\t';

        input_w.open(mfe_variants[j]); //*bigstat_w << mutation_comparison_table[j] << endl << '\t';
        input_w << new_undefined_rna << endl <<
             new_defined_rna   << endl <<
             new_rna_struct << endl <<
             new_dG            << endl <<
             new_cmplx_dG            << endl <<
             new_number_of_links;
        for (int l = 0; l < new_number_of_links; l++)
            input_w << endl << new_HNL[l].first_nucleotide_number_in_link << '\t' <<
                         new_HNL[l].second_nucleotide_number_in_link;
        input_w.close();
        */
    }
    //cout << "ENTER any numb to continue" << endl;
    //cin >> pause;
    //delete [] new_HNL;

    dG_and_struct_reader(defined_rna,
                         &new_rna_struct,
                         HNL,
                         number_of_links);

    *stat_w << (*undefined_rna) << endl << endl;
    *stat_w << (*defined_rna)          << '\t'
          << best_new_dG                 << endl;
    *stat_w << new_rna_struct << endl;

    if ( abs(best_new_dG - target_hairpin_dG) < target_dG_accuracy ||
         best_new_dG > target_hairpin_dG                              )
    {
        (*output_w) << (*defined_rna) << '\t'
                    << best_new_dG  << '\t'
                    << '+' << endl;
        return 1;
    }
    if (break_flag)
    {
        (*output_w) << (*defined_rna) << '\t'
                    << best_new_dG  << '\t'
                    << '-' << endl;
        return 1;
    }

    int break_f = 1;
    for (int f = 0; f < new_undefined_rna.length(); f++)
    {
        if ( old_undefined_rna[f] != (*undefined_rna)[f] )
        {
            break_f = 0;
            break;
        }
    }
    if (break_f)
    {
        for (int k = max_dist_last_HNL_element_number; k >= max_dist_HNL_element_number; k--)
        {
            (*undefined_rna)[ HNL[k].first_nucleotide_number_in_link  - 1] =  'N';
            (*undefined_rna)[ HNL[k].second_nucleotide_number_in_link  - 1] =  'N';
        }
    }
    return 0;
}

int main_work  (string*   defined_rna,
                float*    target_hairpin_dG,
                float*    target_dG_accuracy,
                ofstream* output_w,
                int       dna_flag            )
{
    rna_length = (*defined_rna).length();

    string undefined_rna = (*defined_rna);

    string target_struct;
    for (int i = 0; i <  rna_length; i++)
        target_struct.push_back('.');

    hairpin_nucleotide_links* HNL;
    HNL = new hairpin_nucleotide_links [(int)(rna_length / 2) + 1];
    int number_of_links;

    string defined_rna_struct;
    float initial_dG = dG_and_struct_reader(defined_rna,
                                            &defined_rna_struct,
                                            HNL,
                                            &number_of_links);

    int end_prog_flag = 0;
    char plus_or_minus;
    float new_cmplx_dG;

    if (   abs ( initial_dG - (*target_hairpin_dG) ) < (*target_dG_accuracy) ||
           initial_dG >= (*target_hairpin_dG)                                   )
    {
        result_writer(output_w,
                      defined_rna,
                      dG_reader(defined_rna),
                      '+');
        delete [] HNL;
        return 0;
    }

    ofstream stat_w;
    stat_w.open(statistic);
    stat_w << (*defined_rna) << '\t'
           << (*target_hairpin_dG) << '\t'
           << (*target_dG_accuracy) << endl;

    stat_w << "structure"      << '\t'
          << "hairpin_dG"          << endl;

    stat_w << (*defined_rna)    << '\t'
          << initial_dG          << endl
          << defined_rna_struct  << endl;

    int break_flag = 0;
    //ofstream bigstat_w;
    //bigstat_w.open("./output/bigstat");

    for ( int i = 0; i < max_loop_size; i++ )
    {
        if (anti_hairpin(&undefined_rna,
                         defined_rna,
                         HNL,
                         &number_of_links,
                         &target_struct,
                         *target_hairpin_dG,
                         *target_dG_accuracy,
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
                      string*   initial_rna,
                      float*    target_hairpin_dG,
                      float*    target_dG_accuracy,
                      int*      dna_flag)
{
    char sep = ' ';
    int i = string_before_separator_reader (input_str, initial_rna, sep, 0);

    string temp_string;
    i = string_before_separator_reader (input_str, &temp_string, sep, i);
    decimal_numb_reader( &temp_string, 0, target_hairpin_dG);
    temp_string.clear();

    i = string_before_separator_reader (input_str, &temp_string, sep, i);
    decimal_numb_reader( &temp_string, 0, target_dG_accuracy);
    temp_string.clear();

    float temp_var;
    string_before_separator_reader (input_str, &temp_string, sep, i);
    decimal_numb_reader( &temp_string, 0, &temp_var);
    (*dna_flag) = (int) temp_var;
}

int main()
{
    float target_hairpin_dG,
          target_dG_accuracy;

    string initial_rna,
           input_str;

    ifstream r;
    r.open(input_file);

    ofstream w;
    w.open(output_file);
    w    << "final_olig"   << '\t'
         << "hairpin_dG, kcal/mol" << endl;

    int pause;
    int dna_flag;
    while ( getline(r, input_str) )
    {
        input_rna_reader(&input_str,
                         &initial_rna,
                         &target_hairpin_dG,
                         &target_dG_accuracy,
                         &dna_flag);

        main_work(&initial_rna,
                  &target_hairpin_dG,
                  &target_dG_accuracy,
                  &w,
                  dna_flag);
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
