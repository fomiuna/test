//
// Created by CLAKERS on 2020/8/5.
//

#ifndef MOV_MULTIPASS_H
#define MOV_MULTIPASS_H
#include "onepass.h"
void MultiPass(dataset *data, double B, double eps, double &zep_rev, long long int &zep_ora)
{

    cout << "MultiPass ---------start--------- " << endl;

    long long int oracle_times = 0;
    int node_num = data->index.size();

    double tmp_rev;
    long long int tmp_ora;
    S_gamma S_star = OnePass(data, B, eps, oracle_times, tmp_rev, tmp_ora);
    double gamma_max = 8 * (1 + eps) * S_star.S_revenue / B;
    double alpha = 1.0 / 7.0;
    double gamma_min = alpha * S_star.S_revenue / B;

    vector<S_gamma> Si;
    Si.push_back(S_gamma());
    Si.push_back(S_gamma());
    vector<int> S1_S2_selected(node_num, 0);
    vector<int> L;
    for (double gamma = gamma_max; gamma > (gamma_min / (1 + eps)); gamma /= (1 + eps))
    {
        for (int e = 0; e < node_num; e++)
        {
            //build L
            L.clear();
            if (S1_S2_selected[e] == 0)
            {
                for (int iter = 0; iter < 2; iter++)
                {
                    if (Si[iter].S_cost + data->rate_costs[e] <= B)
                        L.push_back(iter);
                }
            }
            if (!L.empty())
            {
                double max_marginal = -999999999;
                int max_S = -1;
                for (auto i = L.begin(); i != L.end(); i++)
                {
                    oracle_times++;

                    double temp_marginal = Si[*i].marginal(e);
                    if (temp_marginal > max_marginal)
                    {
                        max_marginal = temp_marginal;
                        max_S = *i;
                    }
                }
                if ((max_marginal / data->rate_costs[e]) >= gamma)
                {
                    Si[max_S].Set.push_back(e);
                    Si[max_S].S_cost += data->rate_costs[e];
                    Si[max_S].S_revenue += max_marginal;
                    S1_S2_selected[e] = 1;
                }
            }
        }
    }
    double max_S_value = S_star.S_revenue;
    for (int iter = 0; iter < 2; iter++)
    {
        if (Si[iter].S_revenue > max_S_value)
            S_star.copy(Si[iter]);
    }
    for (int e = 0; e < node_num; e++)
    {
        //cout<<e<<endl;
        for (int iter = 0; iter < 2; iter++)
        {
            //j=0
            //e is not in Sij anyway because Si0 is empty;
            //and the cost of Si0=c(e),f(Si0+e)=f(e)
            if (data->rate_costs[e] <= B)
            {
                oracle_times++;
                double f_value = f_u(e, node_num, data);
                if (f_value > S_star.S_revenue)
                {
                    S_star.Set.clear();
                    S_star.Set.push_back(e);
                    S_star.S_revenue = f_value;
                    S_star.S_cost = data->rate_costs[e];
                }
            }
            //j>=1
            vector<int> Sij_selected(node_num, 0);
            S_gamma Sij;
            for (int j = 0; j < Si[iter].Set.size(); j++)
            {
                //build Sij
                int uij = Si[iter].Set[j]; //Si,j-1 and uij = now Sij
                Sij_selected[uij] = 1;
                Sij.Set.push_back(uij);
                Sij.S_cost += data->rate_costs[uij];
                Sij.S_revenue += Sij.marginal(uij); //just for running quickly, so it is not included in oracle times;

                //e is not in Sij
                if (Sij_selected[e] == 0)
                {
                    if (Sij.S_cost + data->rate_costs[e] <= B)
                    {
                        oracle_times++;
                        double f_Sij_and_e = Sij.marginal(e) + Sij.S_revenue;
                        if (f_Sij_and_e > S_star.S_revenue)
                        {
                            S_star.copy(Sij);
                            S_star.S_revenue = f_Sij_and_e;
                            S_star.S_cost += data->rate_costs[e];
                            S_star.Set.push_back(e);
                        }
                    }
                }
            }
        }
    }

#ifndef ZEP
    cout << "MultiPass & Budget: " << B << endl;
    cout << "S*:" << endl;
    cout << "  revenue: " << S_star.S_revenue << " cost: " << S_star.S_cost << " size: " << S_star.Set.size() << endl;
    for (int i = 0; i < S_star.Set.size(); i++)
        cout << S_star.Set[i] << " ";
    cout << endl;
    cout << "oracle times: " << oracle_times << endl;
    cout << "MultiPass ---------end--------- " << endl
         << endl;
#endif

    zep_ora = oracle_times;
    zep_rev = S_star.S_revenue;
    // zep_cos = S_star.S_cost;
}

#endif //MOV_MULTIPASS_H
