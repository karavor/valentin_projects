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

const char* mfe_variants[4] = {"./input/var_1",
                               "./input/var_2",
                               "./input/var_3",
                               "./input/var_4"};

const int   max_loop_size = 10;
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


struct rna
{
    string seq, seq_struct;
};

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
    $NUPACKHOME/bin/design $CODDING_HOME/input/input_seq \n\
"

#define CONCENTRATIONS_FUNC "\
    #/bin/bash \n\
    export NUPACKHOME=/home/valentin/work_backup/nupack \n\
    export CODDING_HOME=/home/valentin/work_backup/codding \n\
    $NUPACKHOME/bin/complexes -T 23 -material rna $CODDING_HOME/input/complex \n\
    $NUPACKHOME/bin/concentrations $CODDING_HOME/input/complex \n\
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
                            int*    number_of_links,
                            int     struct_flag            )
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
    if (struct_flag)
    {
        getline(r, *seq_struct);
        do
        {
            getline(r, temporal_str);
            links_reader(&HNL[i],&temporal_str);
            i++;
        }
        while (temporal_str[0] != '%');
        *number_of_links = i - 1;
    }
    r.close();
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

int anti_hairpin  (string* undefined_seq,
                   string* defined_seq,
                   string* initial_rna,
                   float   initial_dG,
                   float   initial_complex_dG,
                   hairpin_nucleotide_links*  HNL,
                   int*    number_of_links,
                   string* target_struct,
                   float*  antihairpin_index_,
                   float*  lost_affinity_index_,
                   float   target_hairpin_free_energy_level,
                   float   target_free_enrgy_accuracy,
                   ofstream* stat)
{
    ofstream w;
    int max_dist_HNL_element_number = 0;
    int lack_of_link_groups = *number_of_links;
    int i = 0;
    while (i - max_dist_HNL_element_number == 0) //number of sequence links (one after another)
    {
        if (lack_of_link_groups <= 1)
        {
            w.open(output_file);
            w << *defined_seq << endl <<
              "antihairpin_index = "   << *antihairpin_index_  << endl <<
              "lost_affinity_index = " << *lost_affinity_index_;
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
    }
    string new_undefined_seq,
           new_defined_seq,
           new_defined_seq_struct,
           new_complex;
    int chain_flag;
    float mutation_comparison_table[4],
          antihairpin_index[4],
          lost_affinity_index[4],
          new_dG,
          new_complex_dG;
    ifstream r;
    hairpin_nucleotide_links* new_HNL;
    new_HNL = new hairpin_nucleotide_links [(int)(rna_length / 2) + 1];
    int new_number_of_links;
    for (int j = 0; j < 4; j++)
    {
        new_undefined_seq = *undefined_seq;
        chain_flag = 1;
        for (int k = i - j; k >= max_dist_HNL_element_number; k -= 2)
        {
            if (chain_flag == 1)
            {
                if (new_undefined_seq[ HNL[k].first_nucleotide_number_in_link  - 1] != 'N')
                    new_undefined_seq[ HNL[k].first_nucleotide_number_in_link  - 1] =  'N';
                else
                    new_undefined_seq[ HNL[k].second_nucleotide_number_in_link - 1] =  'N';
            } else
            {
                if (new_undefined_seq[ HNL[k].second_nucleotide_number_in_link - 1] != 'N')
                    new_undefined_seq[ HNL[k].second_nucleotide_number_in_link - 1] =  'N';
                else
                    new_undefined_seq[ HNL[k].first_nucleotide_number_in_link  - 1] =  'N';
            }
            chain_flag *= -1;
        }
        if (j >= 2)
        new_undefined_seq[ HNL[i - (j%2)].second_nucleotide_number_in_link - 1] = 'N';

        w.open(input_design);
        w << *target_struct << endl;
        w << new_undefined_seq;
        w.close();

        system (DESIGN_FUNC);

        r.open(output_design);
        do    {getline(r, new_defined_seq);}
        while (new_defined_seq[0] != 'U' &&
               new_defined_seq[0] != 'C' &&
               new_defined_seq[0] != 'A' &&
               new_defined_seq[0] != 'G'   );
        r.close();

        new_dG = dG_and_struct_reader(&new_defined_seq,
                                      &new_defined_seq_struct,
                                       new_HNL,
                                      &new_number_of_links, 1);

        new_complex = *initial_rna + "+" + new_defined_seq;
        new_complex_dG = dG_and_struct_reader(&new_complex,
                                              &new_complex,
                                              new_HNL,
                                              &new_number_of_links, 0);

        antihairpin_index[j]   = (new_dG - initial_dG) / initial_dG;
        lost_affinity_index[j] = (new_complex_dG - initial_complex_dG) /
                                     initial_complex_dG;
        mutation_comparison_table[j] =                     // how many percents of hairpin destruction
        abs (antihairpin_index[j]/lost_affinity_index[j]); // we have on one percent of affinity loss

        w.open(mfe_variants[j]);
        w << new_undefined_seq << endl <<
             new_defined_seq   << endl <<
             new_defined_seq_struct << endl <<
             new_dG            << endl <<
             new_complex_dG            << endl <<
             new_number_of_links;
        for (int l = 0; l < new_number_of_links; l++)
            w << endl << new_HNL[l].first_nucleotide_number_in_link << '\t' <<
                         new_HNL[l].second_nucleotide_number_in_link;
        w.close();
    }
    delete [] new_HNL;
    i = max_arr_elem_numb( mutation_comparison_table, 4 );

    *antihairpin_index_   = antihairpin_index[i];
    *lost_affinity_index_ = lost_affinity_index[i];

    r.open(mfe_variants[i]);
    getline(r, *undefined_seq);
    getline(r, *defined_seq);
    getline(r, new_defined_seq_struct);

    getline(r, new_complex);
    decimal_numb_reader(&new_complex, 0, &new_dG);
    getline(r, new_complex);
    decimal_numb_reader(&new_complex, 0, &new_complex_dG);

    *stat << *undefined_seq << endl << endl;
    *stat << *defined_seq          << '\t'
          << new_dG                << '\t' << '\t'
          << new_complex_dG        << '\t' << '\t'
          << *antihairpin_index_   << '\t' << '\t'
          << *lost_affinity_index_ << endl;
    *stat << new_defined_seq_struct << endl;

    if ( abs(new_dG -
             target_hairpin_free_energy_level) <
             target_free_enrgy_accuracy   )
    {
        w.open(output_file);
        w << *defined_seq << endl <<
          "antihairpin_index = "   << *antihairpin_index_  << endl <<
          "lost_affinity_index = " << *lost_affinity_index_;
        w.close();
        return 1;
    }
    getline(r, new_complex);
    decimal_numb_reader(&new_complex, 0, &new_dG);
    *number_of_links = (int)new_dG;
    for (int l = 0; l < *number_of_links; l++)
    {
        getline(r, new_complex);
        links_reader(HNL, &new_complex);
    }
    return 0;
}

void main_work (string*   initial_rna,
                float*    target_hairpin_free_energy_level,
                float*    target_free_enrgy_accuracy,
                float*    target_aff,
                float*    target_aff_accuracy)
{
    rna_length = (*initial_rna).length();

    string undefined_rna, defined_rna;
    initial_rna_maker(initial_rna, &defined_rna);
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
                                            &number_of_links, 1);

    string rna_complex = initial_rna + "+" + defined_rna;
    float initial_complex_dG = dG_and_struct_reader(&rna_complex,
                                                    &rna_complex,
                                                    HNL,
                                                    &number_of_links, 0);
    float antihairpin_index = 0,
          lost_affinity_index = 0;

    ofstream w;
    if ( abs(initial_dG -
             *target_hairpin_free_energy_level) <
             target_free_enrgy_accuracy   )
    {
        w.open(output_file);
        w << defined_rna << endl <<
          "antihairpin_index = "   << antihairpin_index  << endl <<
          "lost_affinity_index = " << lost_affinity_index;
        w.close();
        return 0;
    }
    w.open(statistic);
    w << "structure"           << '\t' << '\t' << '\t'
      << "hairpin_dG"          << '\t'
      << "affinity_dG"         << '\t'
      << "antihairpin_index"   << '\t'
      << "lost_affinity_index" << endl;

    w << defined_rna         << '\t' << '\t'
      << initial_dG          << '\t' << '\t'
      << initial_complex_dG  << '\t' << '\t'
      << antihairpin_index   << '\t' << '\t' << '\t'
      << lost_affinity_index << endl;
    w << defined_rna_struct << endl;
    for ( int i = 0; i < max_loop_size; i++ )
    {
        if (anti_hairpin(&undefined_rna,
                         &defined_rna,
                         &initial_rna,
                         initial_dG,
                         initial_complex_dG,
                         HNL,
                         &number_of_links,
                         &target_struct,
                         &antihairpin_index,
                         &lost_affinity_index,
                         target_hairpin_free_energy_level,
                         target_free_enrgy_accuracy,
                         &w) )
            break;
    }
    w.close();
    delete [] HNL;
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
                      float*    target_hairpin_free_energy_level,
                      float*    target_free_enrgy_accuracy,
                      float*    target_aff,
                      float*    target_aff_accuracy)
{
    char sep = ' ';
    int i = string_before_separator_reader (input_str, input_rna, sep, 0);

    string temp_string;
    i = string_before_separator_reader (input_str, &temp_string, sep, i);
    decimal_numb_reader( &temp_string, 0, target_hairpin_free_energy_level);
    temp_string.clear();

    i = string_before_separator_reader (input_str, &temp_string, sep, i);
    decimal_numb_reader( &temp_string, 0, target_free_enrgy_accuracy);
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
    float target_hairpin_free_energy_level,
          target_free_enrgy_accuracy,
          target_aff,
          target_aff_accuracy;

    string input_rna,
           initial_rna;
           input_str;

    ifstream r;
    r.open(input_file);
    while ( getline(r, input_str) )
    {
        input_rna_reader(&input_str,
                         &input_rna,
                         &target_hairpin_free_energy_level,
                         &target_free_enrgy_accuracy,
                         &target_aff,
                         &target_aff_accuracy);
        initial_rna_maker(&input_rna, &initial_rna);
        main_work(&initial_rna,
                  &target_hairpin_free_energy_level,
                  &target_free_enrgy_accuracy,
                  &target_aff,
                  &target_aff_accuracy);
        input_rna.clear();
    }
    r.close();
    return 0;
}
