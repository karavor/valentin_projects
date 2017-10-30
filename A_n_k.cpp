#include <iostream>
#include <string>
#include <vector>

using namespace std;

bool A_n_k_maker ( unsigned short int n,
                   unsigned short int k,
                   vector< vector<unsigned short int> >* A_n_k )
{
    if ( n < k )
        return 1;

    vector<unsigned short int> current_A_n_k_elmnt;

    for (int i = 0; i < n; i++) // arr initialisation
    {
        current_A_n_k_elmnt.push_back(i);
        (*A_n_k).push_back(current_A_n_k_elmnt);

        current_A_n_k_elmnt.clear();
    }

    int f;
    int fix_arr_size;

    for (int l = 0; l < k - 1; l++)
    {
        f = 0;

        fix_arr_size = (*A_n_k).size();

        for (int i = 0; i < fix_arr_size; i++)
        {
            current_A_n_k_elmnt = (*A_n_k)[f];

            (*A_n_k).erase( (*A_n_k).begin() + f );

            for (int j = 0; j < n; j++)
            {
                current_A_n_k_elmnt.push_back(j);
                (*A_n_k).insert( (*A_n_k).begin() + f + j, current_A_n_k_elmnt );
                current_A_n_k_elmnt.erase( current_A_n_k_elmnt.end() - 1 );
            }
            f += n;
        }
    }
    return 0;
}


int main()
{
    vector< vector<unsigned short int> > A_n_k;
    if ( A_n_k_maker (3, 2, &A_n_k) )
    {
        cout << endl << "error in A_n_k_maker: n < k" << endl;
        return 0;
    }

    for (int i = 0; i < A_n_k.size(); i++)
    {
        cout << endl;

        for (int j = 0; j < A_n_k[i].size(); j++)
            cout << A_n_k[i][j] << ' ';
    }
    cout << endl;

    return 0;
}
