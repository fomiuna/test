//
// Created by CLAKERS on 2020/8/2.
//

#ifndef MOV_FANTOM_H
#define MOV_FANTOM_H
#include "onepass.h"
#include <list>
S_gamma GDT(double rho, list<int> &N, double B, dataset *data, long long int &oracle_times)
{
    int node_num = data->index.size();
    //copy omega to N for Greedy algorithm
    double max_marginal = -999999999;
    list<int>::iterator max_node;
    S_gamma S;

    //find z
    int z = -1;
    double max_z = -999999999;
    for (list<int>::iterator p = N.begin(); p != N.end(); p++)
    {
        oracle_times++;

        double temp_f = f_u(*p, data->index.size(), data);
        if (data->rate_costs[*p] <= B)
        {
            if (temp_f > max_z)
            {
                max_z = temp_f;
                z = (*p);
            }
        }
    }
    while (!N.empty())
    {
        max_marginal = -999999999;
        //only visit node in N, which is avaiable now
        for (list<int>::iterator p = N.begin(); p != N.end(); p++)
        {
            //if (S.selected[*p] == 0)//useless，because N always maintain right
            //{

            oracle_times++;

            double temp_marginal = S.marginal(*p);
            if (temp_marginal > max_marginal)
            {
                max_marginal = temp_marginal;
                max_node = p;
            }
            //}
        }
        if (((max_marginal / data->rate_costs[*max_node]) > rho) && (S.S_cost + data->rate_costs[*max_node] <= B))
        {
            S.Set.push_back(*max_node);
            S.S_cost += data->rate_costs[*max_node];
            //S.selected[*max_node]=1;
            S.S_revenue += max_marginal;

            N.erase(max_node++);
        }
        else
        {
            N.erase(max_node++);
        }
    }
    if (max_z > S.S_revenue)
    {
        S.Set.clear();
        S.Set.push_back(z);
        S.S_revenue = max_z;
        S.S_cost = data->rate_costs[z];
        //fill(S.selected.begin(), S.selected.end(), 0);
        //S.selected[z]=1;
    }
    return S;
}
S_gamma USM(const S_gamma &N, dataset *data, long long int &oracle_times)
{
    int node_num = data->index.size();
    S_gamma X;
    S_gamma Y;
    Y.Set.assign(N.Set.begin(), N.Set.end());
    //Y.selected.assign(N.selected.begin(),N.selected.end());
    Y.S_cost = N.S_cost;
    Y.S_revenue = N.S_revenue;
    for (auto u = N.Set.begin(); u != N.Set.end(); u++)
    {
        //calculate ai
        double f_Xi_1 = X.S_revenue;
        //X.selected[*u]=1;
        X.Set.push_back(*u);

        oracle_times++;

        double f_Xi_1_and_u = X.f_S();
        X.Set.pop_back();
        //X.selected[*u]=0;
        double ai = f_Xi_1_and_u - f_Xi_1;

        //calculate bi
        double f_Yi_1 = Y.S_revenue;
        //Y.selected[*u]=0;

        oracle_times++;

        double f_Yi_1_sub_u = Y.S_sub_u(*u);
        //Y.selected[*u]=1;
        double bi = f_Yi_1_sub_u - f_Yi_1;
        if (ai > bi)
        {
            //X.selected[*u]=1;
            X.Set.push_back(*u);
            X.S_revenue = f_Xi_1_and_u;
            X.S_cost += data->rate_costs[*u];
        }
        else
        {
            //Y.selected[*u]=0;
            //delete u_i
            for (vector<int>::iterator p = Y.Set.begin(); p != Y.Set.end();)
            {
                if ((*p) == (*u))
                {
                    Y.Set.erase(p++);
                    break;
                }
                else
                {
                    p++;
                }
            }
            Y.S_revenue = f_Yi_1_sub_u;
            Y.S_cost -= data->rate_costs[*u];
        }
    }
    return X;
}
S_gamma IGDT(double rho, int p, double B, dataset *data, long long int &oracle_times)
{
    int node_num = data->index.size();
    vector<S_gamma> U;
    vector<int> node_available(node_num, 1);
    list<int> omega;
    for (int i = 0; i < node_num; i++) //omega has all nodes at first time
        omega.push_back(i);
    //int available_num=node_num;
    for (int i = 0; i < p + 1; i++)
    {
        S_gamma Si = GDT(rho, omega, B, data, oracle_times);
        S_gamma Si_ = USM(Si, data, oracle_times);
        U.push_back(Si);
        U.push_back(Si_);
        for (int j = 0; j < Si.Set.size(); j++) //mark omega=omega-S_i
        {
            node_available[Si.Set[j]] = 0;
        }
        omega.clear();
        for (int iter = 0; iter < node_num; iter++) //rebuild omega
        {
            if (node_available[iter] == 1)
                omega.push_back(iter);
        }
    }
    double max_value = -999999999;
    S_gamma max_S;
    for (int i = 0; i < U.size(); i++)
    {
        if (U[i].S_revenue > max_value)
        {
            max_S = U[i];
            max_value = U[i].S_revenue;
        }
        /*
        cout<<"U_"<<i<<": "<<endl;
        cout<<"  revenue: "<<U[i].S_revenue<<" cost: "<<U[i].S_cost<<" size: "<<U[i].Set.size()<<endl;
        cout<<"  all nodes: "<<endl;
        for(int j=0;j<U[i].Set.size();j++)
            cout<<U[i].Set[j]<<" ";
        cout<<endl;//*/
    }
    return max_S;
}
void fantom(int p, double B, double eps, dataset *data, double &zep_rev, long long int &zep_ora)
{
    cout << "fantom ---------start--------- " << endl;

    long long int oracle_times = 0;

    int node_num = data->index.size();
    double M = -999999999;
    for (int i = 0; i < node_num; i++)
    {
        oracle_times++;

        double value = f_u(i, node_num, data);
        if (value > M)
        {
            M = value;
        }
    }
    double gamma = 2 * p * M / ((p + 1) * (2 * p + 1));
    int m = 0;
    vector<double> R;
    while (1)
    {
        if (pow(1 + eps, m) < node_num)
            R.push_back(pow(1 + eps, m) * gamma);
        else
        {
            R.push_back(gamma * node_num);
            break;
        }
        m++;
    }
    //for rho\in R
    vector<S_gamma> U;
    for (int i = 0; i < R.size(); i++)
    {
        //cout<<"rho: "<<R[i]<<endl;
        S_gamma S = IGDT(R[i], p, B, data, oracle_times);
        U.push_back(S);
    }

    double max_value = -999999999;
    S_gamma max_S;
    for (int i = 0; i < U.size(); i++)
    {
        if (U[i].S_revenue > max_value)
        {
            max_S = U[i];
            max_value = U[i].S_revenue;
        }
    }

#ifndef ZEP
    cout << "fantom & Budget: " << B << endl;
    cout << "S*:" << endl;
    cout << "  revenue: " << max_S.S_revenue << " cost: " << max_S.S_cost << " size: " << max_S.Set.size() << endl;
    for (int i = 0; i < max_S.Set.size(); i++)
        cout << max_S.Set[i] << " ";
    cout << endl;
    cout << "oracle times: " << oracle_times << endl;
    cout << "fantom ---------end--------- " << endl
         << endl;
#endif

    zep_ora = oracle_times;
    zep_rev = max_S.S_revenue;
}

#endif //MOV_FANTOM_H
