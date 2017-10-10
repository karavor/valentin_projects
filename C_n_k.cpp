#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <cmath>
#include <vector>

using namespace std;

template<class _T_>

bool C_n_k_generator (vector<_T_>* elements_arr,
                      vector< vector<_T_> >* result_arr,
                      int k)
{
    int n = (*elements_arr).size();

    if (n < k)
        return 1;

    vector<_T_> current_res;

    for (int j = 0; j <= n - k; j++)
    {
        for (int l = j; l < j + k; l++)
            current_res.push_back( (*elements_arr)[l] );

        (*result_arr).push_back(current_res);

        for (int i = j + k; i < n; i++)
        {
            current_res.erase(current_res.begin() + k - 1);
            current_res.push_back( (*elements_arr)[i] );
            (*result_arr).push_back(current_res);
        }
        current_res.clear();
    }
    return 0;
}

bool C_n_k_generator_caller ()
{
    vector<int> elements_arr;
    elements_arr.push_back(1);
    elements_arr.push_back(2);
    elements_arr.push_back(3);
    elements_arr.push_back(4);
    elements_arr.push_back(5);

    vector< vector<int> > result_arr;
    int k = 4;

    if ( C_n_k_generator<int> (&elements_arr, &result_arr, k) )
    {
        cout << endl << "fail: n < k" << endl;
        return 1;
    }

    for (int i = 0; i < result_arr.size(); i++)
    {
        for (int j = 0; j < k; j++)
        {
            cout << result_arr[i][j] << ' ';
        }
        cout << endl;
    }
    return 0;
}

int main()
{
    C_n_k_generator_caller();
    return 0;
}
