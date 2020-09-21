//
// Created by CLAKERS on 2020/8/4.
//

#ifndef MOV_SAMPLEGREEDY_H
#define MOV_SAMPLEGREEDY_H
#include "onepass.h"
#include "random"

default_random_engine e(1234);
S_gamma samplegreedy(double B, double p, dataset *data, long long int &oracle_times)
{
    bernoulli_distribution u(p);

    int node_num = data->index.size();
    int i_star = -1;
    double i_star_value = -999999999;
    for (int iter = 0; iter < node_num; iter++)
    {
        oracle_times++;

        double value = f_u(iter, node_num, data);
        if (data->rate_costs[iter] <= B)
        {
            if (value > i_star_value)
            {
                i_star_value = value;
                i_star = iter;
            }
        }
    }
    S_gamma S;
    vector<int> F;
    double R = B;
    for (int iter = 0; iter < node_num; iter++)
    {
        oracle_times++;

        double value = f_u(iter, node_num, data);
        if (value > 0 && data->rate_costs[iter] <= B)
        {
            F.push_back(iter);
        }
    }
    vector<int> A_available(node_num, 1);
    while (!F.empty())
    {
        //argmax i
        int node_i = -1;
        double max_marginal_cost = -999999999;
        for (vector<int>::iterator k = F.begin(); k != F.end(); k++)
        {
            oracle_times++;

            double temp_marginal_cost = S.marginal(*k) / data->rate_costs[*k];
            if (temp_marginal_cost > max_marginal_cost)
            {
                max_marginal_cost = temp_marginal_cost;
                node_i = *k;
            }
        }
        double random = u(e);
        //cout<<"random: "<<random<<endl;
        if (random == 1)
        {
            S.Set.push_back(node_i);
            S.S_cost += data->rate_costs[node_i];
            S.S_revenue += max_marginal_cost * data->rate_costs[node_i];
            R -= data->rate_costs[node_i];
        }
        A_available[node_i] = 0;
        F.clear();
        for (int iter = 0; iter < node_num; iter++)
        {
            if (A_available[iter] == 1)
            {
                oracle_times++;

                double value = S.marginal(iter);
                if (value > 0 && data->rate_costs[iter] <= R)
                {
                    F.push_back(iter);
                }
            }
        }
    }
    if (i_star_value > S.S_revenue)
    {
        S.Set.clear();
        S.Set.push_back(i_star);
        S.S_revenue = i_star_value;
        S.S_cost = data->rate_costs[i_star];
    }
    return S;
}

void iter_sample(double B, double default_p, dataset *data, double &zep_rev, long long int &zep_ora)
{

    cout << "SampleGreedy ---------start--------- " << endl;

    long long int oracle_times = 0;
#ifndef ZEP
    cout << "SampleGreedy & B: " << B << endl;
#endif
    S_gamma max_S = samplegreedy(B, default_p, data, oracle_times);

#ifndef ZEP
    cout << "S*:" << endl;
    cout << "  revenue: " << max_S.S_revenue << " cost: " << max_S.S_cost << " size: " << max_S.Set.size() << endl;
    for (int i = 0; i < max_S.Set.size(); i++)
        cout << max_S.Set[i] << " ";
    cout << endl;
    cout << "oracle times: " << oracle_times << endl;
    cout << "SampleGreedy ---------end--------- " << endl;
#endif

    zep_ora = oracle_times;
    zep_rev = max_S.S_revenue;
}

#endif //MOV_SAMPLEGREEDY_H
