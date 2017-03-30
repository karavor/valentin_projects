#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <cmath>

using namespace std;

int pause;

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

char random_undef_ncltd (int rnd_seed)
{
    srand (rnd_seed);

    switch(rand() % 4)
    {
        case 0:
            return 'D';
        case 1:
            return 'H';
        case 2:
            return 'B';
        case 3:
            return 'V';
    }
}

template <typename T>

void print_arr (T* arr, int length)
{
    for (int i = 0; i < length; i++)
        cout << arr[i] << ' ';
    cout << endl;
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

void n_pow_k_comb_print (int n, int k)
{
    int comb_nmb = (int)pow(n, k);
    int n_pow_k_arr [k];

    for (int i = 0; i < k; i++)
        n_pow_k_arr[i] = 0;

    for (int i = 0; i < comb_nmb; i++)
    {
        print_arr <int> (n_pow_k_arr, k);
        n_pow_k_arr[0]++;

        for (int l = 0; l < (k - 1); l++)
        {
            if ( left_decade_pop( &(n_pow_k_arr[l + 1]),
                                  &(n_pow_k_arr[l]),
                                  (n - 1)               ) )
                break;
        }
    }
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

int main()
{

    string str = "AGCT",
           except_str;
    int n = 3,
        k = str.length();

    int comb_nmb = (int)pow(n, k);
    int n_pow_k_arr [k];

    for (int i = 0; i < k; i++)
        n_pow_k_arr[i] = 0;

    for (int i = 0; i < comb_nmb; i++)
    {
        for (int j = 0; j < k; j++)
        {
            except_str.push_back( except_defined_ncltd (str[j], n_pow_k_arr[j], 1) );
        }
        cout << except_str << endl;
        except_str.clear();
        n_pow_k_arr[0]++;

        for (int l = 0; l < (k - 1); l++)
        {
            if ( left_decade_pop( &(n_pow_k_arr[l + 1]),
                                  &(n_pow_k_arr[l]),
                                  (n - 1)               ) )
                break;
        }
    }
    return 0;
}
