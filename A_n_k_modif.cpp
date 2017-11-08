#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <cmath>
#include <vector>
#include <iomanip>

using namespace std;

const char* test_output    = "./test_output";

int int_pow (int n, int power)
{
    int result = 1;

    for (int i = 0; i < power; i++)
        result *= n;

    return result;
}

struct A_n_k_elmnt_with_nmb
{
    unsigned short int nmb;
    vector<unsigned short int> elmnt;
};

bool A_n_k_elmnt_test ( A_n_k_elmnt_with_nmb* A_n_k_elmnt )
{
    //vector <unsigned short int> frbdn_vect_0;
    //frbdn_vect_0.push_back(1);

    vector <unsigned short int> frbdn_vect_1;
    frbdn_vect_1.push_back(1);
    frbdn_vect_1.push_back(0);

    if ( ( (*A_n_k_elmnt).elmnt == frbdn_vect_1 ) )
        //||
         //( (*A_n_k_elmnt).elmnt == frbdn_vect_0 )) //for example
    {
        return 0;
    }

    return 1;
}

bool forbidden_A_n_k_elmnt_nmb_test( unsigned short int n,
                                     unsigned short int current_A_n_k_elmnt_nmb,
                                     vector< vector<unsigned short int> >* forbidden_A_n_k_elmnts_nmbs_stat )
{
    for (int i = 0; i < (*forbidden_A_n_k_elmnts_nmbs_stat).size(); i++)
    {
        for (int j = 0; j < (*forbidden_A_n_k_elmnts_nmbs_stat)[i].size(); j++)
        {
            if ( ( current_A_n_k_elmnt_nmb -
                   (*forbidden_A_n_k_elmnts_nmbs_stat)[i][j] ) % int_pow(n, i + 1) == 0 )
            {
                return 0;
            }
        }
    }
    return 1;
}

template<typename _DT_>
void print_vect_as_string (vector<_DT_>* vect,
                           ofstream* write    )
{
    for (int i = 0; i < (*vect).size(); i++)
    {
        (*write) << (*vect)[i];
    }
}

void A_n_k_maker_stat_writer ( ofstream* write,
                               int iter_nmb,
                               vector< A_n_k_elmnt_with_nmb >* A_n_k,
                               vector< vector<unsigned short int> >*
                                     forbidden_A_n_k_elmnts_nmbs_stat)
{
    (*write) << endl << iter_nmb << endl << endl;

    char sep = '\t';

    for (int i = 0; i < (*A_n_k).size(); i++)
    {
        (*write) << sep << (*A_n_k)[i].nmb << sep;

        print_vect_as_string <unsigned short int> ( &((*A_n_k)[i].elmnt),
                                                     write                );
        (*write) << endl;
    }
    (*write) << endl;

    for (int j = 0; j < (*forbidden_A_n_k_elmnts_nmbs_stat)[iter_nmb].size(); j++)
    {
        (*write) << sep << (*forbidden_A_n_k_elmnts_nmbs_stat)[iter_nmb][j] << endl;
    }
}

void dead_end_A_n_k_elmnts_nmbs_add (int  n,
                                     int  i,
                                     int  A_n_k_insert_nmb_nmb,
                                     int* last_erased_elmnt_nmb,
                                     int  fix_arr_size,
                                     int  iter_nmb,
                                     vector< vector<unsigned short int> >*
                                     forbidden_A_n_k_elmnts_nmbs_stat)
{
    if ( (i == 0) && ( A_n_k_insert_nmb_nmb > 0) )
    {
        (*last_erased_elmnt_nmb) = -1;
    }

    if ( (i == (fix_arr_size - 1))
        &&
         ( A_n_k_insert_nmb_nmb < (int_pow(n, iter_nmb) - 1) ) )
    {
        for (int j = ( A_n_k_insert_nmb_nmb + 1 ) * n;
                 j < int_pow(n, iter_nmb)               * n; j++)
        {
            (*forbidden_A_n_k_elmnts_nmbs_stat)
            [(*forbidden_A_n_k_elmnts_nmbs_stat).size() - 1].
            push_back( j );
        }
    }

    if ( ( A_n_k_insert_nmb_nmb - 1 ) > (*last_erased_elmnt_nmb) )
    {
        for (int j = ((*last_erased_elmnt_nmb) + 1) * n;
                 j < A_n_k_insert_nmb_nmb           * n; j++)
        {
            (*forbidden_A_n_k_elmnts_nmbs_stat)
            [(*forbidden_A_n_k_elmnts_nmbs_stat).size() - 1].
            push_back( j );
        }
    }
}

void A_n_k_arr_init (int n,
                     vector< A_n_k_elmnt_with_nmb >* A_n_k,
                     vector< vector<unsigned short int> >*
                     forbidden_A_n_k_elmnts_nmbs_stat)
{
    A_n_k_elmnt_with_nmb current_A_n_k_elmnt;
    vector<unsigned short int> init_vect;

    (*forbidden_A_n_k_elmnts_nmbs_stat).push_back(init_vect);

    for (unsigned short int i = 0; i < n; i++) // arr initialisation
    {
        current_A_n_k_elmnt.elmnt.push_back(i);
        current_A_n_k_elmnt.nmb = i;

        if ( A_n_k_elmnt_test( &current_A_n_k_elmnt ) )
        {
            (*A_n_k).push_back(current_A_n_k_elmnt);
        }
        else (*forbidden_A_n_k_elmnts_nmbs_stat)
             [(*forbidden_A_n_k_elmnts_nmbs_stat).size() - 1].
             push_back( current_A_n_k_elmnt.nmb );

        current_A_n_k_elmnt.elmnt.clear();
    }
}

bool modified_A_n_k_maker ( unsigned short int n,
                            unsigned short int k,
                            vector< A_n_k_elmnt_with_nmb >* A_n_k )
{
    if ( n < k )
        return 1;

    vector< vector<unsigned short int> > forbidden_A_n_k_elmnts_nmbs_stat;

    A_n_k_arr_init (n,
                    A_n_k,
                    &forbidden_A_n_k_elmnts_nmbs_stat);

    ofstream write;
    write.open(test_output);

    A_n_k_maker_stat_writer(&write,
                            0,
                            A_n_k,
                            &forbidden_A_n_k_elmnts_nmbs_stat);

    int fix_arr_size;
    int insert_nmb;
    int last_erased_elmnt_nmb;
    A_n_k_elmnt_with_nmb current_A_n_k_elmnt;
    vector<unsigned short int> init_vect;

    for (int iter_nmb = 1; iter_nmb < k; iter_nmb++)
    {
        insert_nmb = 0;
        current_A_n_k_elmnt.nmb = 0;
        last_erased_elmnt_nmb = 0;

        fix_arr_size = (*A_n_k).size();

        forbidden_A_n_k_elmnts_nmbs_stat.push_back(init_vect);

        for (int i = 0; i < fix_arr_size; i++)
        {
            current_A_n_k_elmnt.elmnt = ((*A_n_k)[insert_nmb]).elmnt;

            current_A_n_k_elmnt.nmb = n * ((*A_n_k)[insert_nmb]).nmb;

            dead_end_A_n_k_elmnts_nmbs_add (n,
                                            i,
                                            ((*A_n_k)[insert_nmb]).nmb,
                                            &last_erased_elmnt_nmb,
                                            fix_arr_size,
                                            iter_nmb,
                                            &forbidden_A_n_k_elmnts_nmbs_stat);

            last_erased_elmnt_nmb = ((*A_n_k)[insert_nmb]).nmb;

            (*A_n_k).erase( (*A_n_k).begin() + insert_nmb );

            for (int j = 0; j < n; j++)
            {
                if ( forbidden_A_n_k_elmnt_nmb_test( n,
                                                     current_A_n_k_elmnt.nmb,
                                                     &forbidden_A_n_k_elmnts_nmbs_stat ) )
                {
                    current_A_n_k_elmnt.elmnt.push_back(j);

                    if ( A_n_k_elmnt_test( &current_A_n_k_elmnt ) )
                    {
                        (*A_n_k).insert( (*A_n_k).begin() + insert_nmb,
                                  current_A_n_k_elmnt);
                        insert_nmb ++;
                    }
                    else forbidden_A_n_k_elmnts_nmbs_stat
                         [forbidden_A_n_k_elmnts_nmbs_stat.size() - 1].
                         push_back( current_A_n_k_elmnt.nmb );

                    current_A_n_k_elmnt.elmnt.
                    erase( current_A_n_k_elmnt.elmnt.end() - 1 );
                }

                current_A_n_k_elmnt.nmb ++;
            }
        }
        A_n_k_maker_stat_writer(&write,
                                iter_nmb,
                                A_n_k,
                                &forbidden_A_n_k_elmnts_nmbs_stat);
    }
    write.close();
    return 0;
}

int main()
{
    unsigned short int n = 3,
                       k = 3;

    vector< A_n_k_elmnt_with_nmb > A_n_k;

    modified_A_n_k_maker ( n,
                           k,
                           &A_n_k );

    return 0;
}
