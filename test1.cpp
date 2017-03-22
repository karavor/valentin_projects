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
const int GEN_TERM  = 1;//ok, terminate
const int GEN_EMPTY = 2; //ok, print EMPTY SET and continue
const int GEN_ERROR = 3;
int pause;

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
    $NUPACKHOME/bin/complexes -T 23 -material dna -quiet -mfe $PWD/input/complex \n\
    $NUPACKHOME/bin/concentrations -quiet $PWD/input/complex \n\
"

struct hairpin_nucleotide_links
{
    float first_nucleotide_number_in_link;
    float second_nucleotide_number_in_link;
    float distance_between_nucleotides;
};

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

char random_ncltd (int dna_flag, int rnd_seed)
{
    srand (rnd_seed);

    switch(rand() % 4)
    {
        case 0:
            return 'A';
        case 1:
            if (dna_flag)
                return 'T';
            return 'U';
        case 2:
            return 'G';
        case 3:
            return 'C';
    }
}

void olig_generator (int olig_length,
                     string* olig_seq,
                     int dna_flag)
{
    for (int i = 0; i < olig_length; i++)

        (*olig_seq).push_back( random_ncltd(dna_flag, time(NULL) + i) );
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


float olig_generator_with_pair (int     initial_olig_length,
                                string* olig_seq,
                                string* c_olig_seq,
                                float   target_cmplx_dG,
                                int     dna_flag)
{
    int olig_length = initial_olig_length;
    string cmplx;
    float cmplx_dG = 0;
    olig_generator( initial_olig_length, olig_seq, dna_flag);
    c_seq_maker (olig_seq, c_olig_seq, dna_flag);
    cmplx = (*olig_seq) + "+" + (*c_olig_seq);
    cmplx_dG = dG_reader(&cmplx);

    char new_ncltd,
         c_new_ncltd;
    while (cmplx_dG > target_cmplx_dG)
    {
        new_ncltd = random_ncltd(dna_flag, time(NULL) + olig_length);
        (*olig_seq).push_back( new_ncltd );
        olig_length++;


        c_new_ncltd = complement_ncltd(new_ncltd, dna_flag);
        (*c_olig_seq).insert( (*c_olig_seq).begin(), 1, c_new_ncltd );

        cmplx.insert( cmplx.begin() + olig_length - 1, 1, new_ncltd   );
        cmplx.insert( cmplx.begin() + olig_length + 1, 1, c_new_ncltd );

        cmplx_dG = dG_reader(&cmplx);
    }
    return cmplx_dG;
}

int simple_total_antihairpin ( string* new_defined_rna,
                               string* undefined_rna,
                               string* target_struct,
                               int*    number_of_links,
                               hairpin_nucleotide_links*  HNL,
                               int     dna_flag)
{
    for (int i = 0; i < *number_of_links; i++)
    {
        if ( (*undefined_rna)[HNL[i].first_nucleotide_number_in_link - 1] == 'N' )
            (*undefined_rna)[HNL[i].second_nucleotide_number_in_link - 1] = 'N';
        else (*undefined_rna)[HNL[i].first_nucleotide_number_in_link - 1] = 'N';

        design_maker(target_struct,
                     undefined_rna,
                     new_defined_rna,
                     dna_flag);

        if ( dG_reader(new_defined_rna) == 0 )
            return 1;
    }

    string temp_str;
    if (dG_and_struct_reader (new_defined_rna,
                              &temp_str,
                              HNL,
                              number_of_links)  < 0)
    {
        simple_total_antihairpin(new_defined_rna,
                          undefined_rna,
                          target_struct,
                          number_of_links,
                          HNL,
                          dna_flag);
    }
    return 0;
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

int main()
{
    // DESIGN_FUNC MFE_FUNC CONCENTRATIONS_FUNC
    system(CONCENTRATIONS_FUNC);

    //string olig;
    //olig_generator(12, &olig, 1);
/*
    string olig = "ACCCTATTCAGG";
    string olig_struct;
    string undef_olig = olig;
    string target_struct;

    for (int i = 0; i < olig.length(); i++)
        target_struct.push_back('.');

    int links_nmb;

    hairpin_nucleotide_links HNL [6];

    dG_and_struct_reader(&olig, &olig_struct, HNL, &links_nmb);

    simple_total_antihairpin ( &olig,
                               &undef_olig,
                               &target_struct,
                               &links_nmb,
                               HNL,
                               1);
/*
    string olig = "AAGATCGTTGGGCAAGGGGAAAGTAC";
    string target_olig = "GTACTTTCCCCTTGCCCAACGATCGC";

    float cmplx_dG =
    cmplx_dG_reducer(-30,
                     3,
                     &olig,          // here is resulting rna
                     &target_olig);
 /*
    string olig_seq,
           c_olig_seq;

    float cmplx_dG =
    olig_generator_with_pair (15,
                              &olig_seq,
                              &c_olig_seq,
                              -45,
                              1);

    ofstream w;
    w.open("./input/except_olig");
    w << olig;
    w.close();
*/
    return 0;
}
