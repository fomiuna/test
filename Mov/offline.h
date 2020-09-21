//
// Created by amax on 2020/7/31.
//

#ifndef MOVREC_OFFLINE_H
#define MOVREC_OFFLINE_H

#include <iostream>
#include <fstream>
#include "random"
#include "set"
#include "utility_functions.h"

using namespace std;

double f_S(const vector<int> &Set, dataset *data)
{
    double sum_value = 0.0;
    float M1 = 0, M2 = 0;
    for (int it = 0; it < Set.size(); it++)
    {
        M1 += sim_mat_sum[Set[it]];
        for (int its = 0; its < Set.size(); its++)
            M2 += sim_mat[Set[it]][Set[its]];
    }
    sum_value = M1 - M2;
    return sum_value;
}
double marginal_gain(int e, const vector<int> &Set, dataset *data)
{
    float M1 = 0, M2 = 0;
    M1 += sim_mat_sum[e];
    for (int it = 0; it < Set.size(); it++)
    {
        M2 += sim_mat[e][Set[it]];
        M2 += sim_mat[Set[it]][e];
    }
    M2 += sim_mat[e][e];
    return M1 - M2;
}
void Offline(dataset *data, double B, double &zep_rev, long long int &zep_ora)
{

    cout << "Offline ---------start--------- " << endl;

    long long int oracle_times = 0;
    int node_num = data->index.size();

    //test value
    /*
    vector<int> t1={1754, 50, 169 ,273 ,330 ,508, 743, 770 ,860 ,980 ,1072, 1563, 1518};
    cout<<f_S(t1,data)<<endl;
    double sum_cost=0.0;
    for(int i=0;i<t1.size();i++)
        sum_cost+=data->rate_costs[t1[i]];
    cout<<sum_cost<<endl;

    vector<int> t2={50, 159 ,174, 316, 333, 449 ,669 ,750 ,860, 1754};
    cout<<f_S(t2,data)<<endl;

    vector<int> t3={860, 316, 330 ,50 ,1398, 743, 1092 ,1754 ,174, 1563, 1761};
    cout<<f_S(t3,data)<<endl;
//*/

    double S1_cost = 0.0;
    double S2_cost = 0.0;
    double S3_cost = 0.0;
    double S1_revenue = 0.0;
    double S2_revenue = 0.0;
    double S3_revenue = 0.0;

    vector<int> S1, S2, S3;
    //use for fast calculate f_S, 0->not selected,1->selected
    // vector<int> S1_selected(node_num,0);
    // vector<int> S2_selected(node_num,0);
    // vector<int> S3_selected(node_num,0);

    vector<int> S1_S2_selected(node_num, 0);
    vector<int> C1, C2;
    while (1)
    {
        C1.clear();
        C2.clear();
        for (int i = 0; i < node_num; i++)
        {
            if (S1_S2_selected[i] == 0)
            {
                if (S1_cost + data->rate_costs[i] <= B)
                {
                    C1.push_back(i);
                }
                if (S2_cost + data->rate_costs[i] <= B)
                {
                    C2.push_back(i);
                }
            }
        }
        if (!C1.empty() && !C2.empty()) //L={1,2}
        {
            //C1
            int a1 = -1;
            double max_marginal_cost1 = -999999999; //max f(e|S1)/c(e)
            int b1 = -1;
            double max_marginal1 = -999999999; //max f(e|S1)
            //double f_S_value1=f_S(S1_selected,S1);

            //arg max e\in C_i
            for (int i = 0; i < C1.size(); i++)
            {
                //S1_selected[C1[i]]=1;
                //S1.push_back(C1[i]);

                oracle_times++;

                double marginal = marginal_gain(C1[i], S1, data);          //f(e|S1)
                double marginal_cost = marginal / data->rate_costs[C1[i]]; //f(e|S1)/c(e)

                if (marginal_cost > max_marginal_cost1)
                {
                    max_marginal_cost1 = marginal_cost;
                    a1 = C1[i];
                }
                if (marginal > max_marginal1)
                {
                    max_marginal1 = marginal;
                    b1 = C1[i];
                }
                //S1_selected[C1[i]]=0;
                //S1.pop_back();
            }
            double f_S_union_e1 = max_marginal1 + S1_revenue; //f(S1+b1)
            if (f_S_union_e1 > S3_revenue)
            {
                S3.assign(S1.begin(), S1.end());
                //S3_selected.assign(S1_selected.begin(),S1_selected.end());
                S3.push_back(b1);
                //S3_selected[b1]=1;

                //cout<<"S3 select node: "<<b1<<" from b1"<<endl;

                S3_revenue = f_S_union_e1;
                S3_cost = S1_cost + data->rate_costs[b1];
            }
            //C2
            int a2 = -1;
            double max_marginal_cost2 = -999999999;
            int b2 = -1;
            double max_marginal2 = -999999999;

            //arg max e\in C_i
            for (int i = 0; i < C2.size(); i++)
            {
                //S2_selected[C2[i]]=1;
                //S2.push_back(C2[i]);

                oracle_times++;

                double marginal = marginal_gain(C2[i], S2, data);          //f(e|S2)
                double marginal_cost = marginal / data->rate_costs[C2[i]]; //f(e|S2)/c(e)
                if (marginal_cost > max_marginal_cost2)
                {
                    max_marginal_cost2 = marginal_cost;
                    a2 = C2[i];
                }
                if (marginal > max_marginal2)
                {
                    max_marginal2 = marginal;
                    b2 = C2[i];
                }
                //S2_selected[C2[i]]=0;
                //S2.pop_back();
            }
            double f_S_union_e2 = max_marginal2 + S2_revenue; //f(S2+b2)
            if (f_S_union_e2 > S3_revenue)
            {
                S3.assign(S2.begin(), S2.end());
                // S3_selected.assign(S2_selected.begin(),S2_selected.end());
                S3.push_back(b2);
                // S3_selected[b2]=1;
                //cout<<"S3 select node: "<<b2<<" from b2"<<endl;

                S3_revenue = f_S_union_e2;
                S3_cost = S2_cost + data->rate_costs[b2];
            }

            if (max_marginal_cost1 > max_marginal_cost2) //j=1,S1
            {
                double a1_marginal = max_marginal_cost1 * data->rate_costs[a1];
                if (a1_marginal > 0)
                {
                    S1.push_back(a1);
                    // S1_selected[a1]=1;
                    S1_S2_selected[a1] = 1;
                    S1_revenue += a1_marginal;
                    S1_cost += data->rate_costs[a1];

                    //cout<<"S1 select node: "<<a1<<endl;
                }
                else
                {
                    C1.clear();
                    C2.clear();
                }
            }
            else
            { //j=2,S2
                double a2_marginal = max_marginal_cost2 * data->rate_costs[a2];
                if (a2_marginal > 0)
                {
                    S2.push_back(a2);
                    //  S2_selected[a2]=1;
                    S1_S2_selected[a2] = 1;
                    S2_revenue += a2_marginal;
                    S2_cost += data->rate_costs[a2];

                    //cout<<"S2 select node: "<<a2<<endl;
                }
                else
                {
                    C1.clear();
                    C2.clear();
                }
            }
        }
        else if (!C1.empty() && C2.empty()) //L={1}
        {
            int a1 = -1;
            double max_marginal_cost1 = -999999999;
            int b1 = -1;
            double max_marginal1 = -999999999;
            //double f_S_value1=f_S(S1_selected,S1);

            //arg max e\in C_i
            for (int i = 0; i < C1.size(); i++)
            {
                // S1_selected[C1[i]]=1;
                //S1.push_back(C1[i]);

                // double marginal=f_S(S1_selected)-S1_revenue;//f(e|S1)

                oracle_times++;

                double marginal = marginal_gain(C1[i], S1, data);          //f(e|S1)
                double marginal_cost = marginal / data->rate_costs[C1[i]]; //f(e|S1)/c(e)
                if (marginal_cost > max_marginal_cost1)
                {
                    max_marginal_cost1 = marginal_cost;
                    a1 = C1[i];
                }
                if (marginal > max_marginal1)
                {
                    max_marginal1 = marginal;
                    b1 = C1[i];
                }
                //S1_selected[C1[i]]=0;
                //S1.pop_back();
            }
            double f_S_union_e1 = max_marginal1 + S1_revenue; //f(S1+b1)
            if (f_S_union_e1 > S3_revenue)
            {
                S3.assign(S1.begin(), S1.end());
                //S3_selected.assign(S1_selected.begin(),S1_selected.end());
                S3.push_back(b1);
                // S3_selected[b1]=1;

                S3_revenue = f_S_union_e1;
                S3_cost = S1_cost + data->rate_costs[b1];

                //cout<<"S3 select node: "<<b1<<" from b1"<<endl;
            }
            double a1_marginal = max_marginal_cost1 * data->rate_costs[a1];
            if (a1_marginal > 0)
            {
                S1.push_back(a1);
                //S1_selected[a1]=1;
                S1_S2_selected[a1] = 1;
                S1_revenue += a1_marginal;
                S1_cost += data->rate_costs[a1];

                //cout<<"S1 select node: "<<a1<<endl;
            }
            else
            {
                C1.clear();
                C2.clear();
            }
        }
        else if (C1.empty() && !C2.empty()) //L={2}
        {
            int a2 = -1;
            double max_marginal_cost2 = -999999999;
            int b2 = -1;
            double max_marginal2 = -999999999;

            //arg max e\in C_i
            for (int i = 0; i < C2.size(); i++)
            {
                // S2_selected[C2[i]]=1;
                //S2.push_back(C2[i]);

                oracle_times++;

                double marginal = marginal_gain(C2[i], S2, data);          //f(e|S2)
                double marginal_cost = marginal / data->rate_costs[C2[i]]; //f(e|S2)/c(e)
                if (marginal_cost > max_marginal_cost2)
                {
                    max_marginal_cost2 = marginal_cost;
                    a2 = C2[i];
                }
                if (marginal > max_marginal2)
                {
                    max_marginal2 = marginal;
                    b2 = C2[i];
                }
                //S2_selected[C2[i]]=0;
                //S2.pop_back();
            }
            double f_S_union_e2 = max_marginal2 + S2_revenue; //f(S2+b2)
            if (f_S_union_e2 > S3_revenue)
            {
                S3.assign(S2.begin(), S2.end());
                //S3_selected.assign(S2_selected.begin(),S2_selected.end());
                S3.push_back(b2);
                //S3_selected[b2]=1;

                S3_revenue = f_S_union_e2;
                S3_cost = S2_cost + data->rate_costs[b2];
                //cout<<"S3 select node: "<<b2<<" from b2"<<endl;
            }
            double a2_marginal = max_marginal_cost2 * data->rate_costs[a2];
            if (a2_marginal > 0)
            {
                S2.push_back(a2);
                //S2_selected[a2]=1;
                S1_S2_selected[a2] = 1;
                S2_revenue += a2_marginal;
                S2_cost += data->rate_costs[a2];

                //cout<<"S3 select node: "<<a2<<endl;
            }
            else
            {
                C1.clear();
                C2.clear();
            }
        }

        //cout<<"C1 size: "<<C1.size()<<endl;
        //cout<<"C2 size: "<<C2.size()<<endl;

        if (C1.empty() && C2.empty())
            break;
    }
#ifndef ZEP
    cout << "Offline & Budget: " << B << endl;

    cout << "S1:" << endl;
    cout << "  revenue: " << S1_revenue << " cost: " << S1_cost << " size: " << S1.size() << endl;
    cout << "  all nodes: " << endl;
    for (int i = 0; i < S1.size(); i++)
    {
        cout << S1[i] << " ";
    }
    cout << endl;
    cout << "S2:" << endl;
    cout << "  revenue: " << S2_revenue << " cost: " << S2_cost << " size: " << S2.size() << endl;
    cout << "  all nodes: " << endl;
    for (int i = 0; i < S2.size(); i++)
    {
        cout << S2[i] << " ";
    }
    cout << endl;
    cout << "S3:" << endl;
    cout << "  revenue: " << S3_revenue << " cost: " << S3_cost << " size: " << S3.size() << endl;
    cout << "  all nodes: " << endl;
    for (int i = 0; i < S3.size(); i++)
    {
        cout << S3[i] << " ";
    }
    cout << endl;
#endif

    double best_revenue = max(S1_revenue, S2_revenue);
    best_revenue = max(best_revenue, S3_revenue);

#ifndef ZEP
    cout << "best revenue: " << best_revenue << " oracle times: " << oracle_times << endl;
    cout << "Offline ---------end--------- " << endl
         << endl;
#endif

    zep_ora = oracle_times;
    zep_rev = best_revenue;
}

#endif //MOVREC_OFFLINE_H
