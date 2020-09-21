
#ifndef STREAMING_ALGORITHM_ONEPASS_H
#define STREAMING_ALGORITHM_ONEPASS_H
#include <iostream>
#include <fstream>
#include <cstring>
#include <cstdlib>
#include <sys/types.h>
#include <sstream>
#include <climits>
#include <string>
#include <vector>
#include <cmath>
#include <set>
#include "offline.h"

using namespace std;

//int node_num;

class S_gamma
{
public:
    S_gamma()
    {
        //sim_mat = simmat;
        //node_num = nodenum;
        S_cost = 0.0;
        S_revenue = 0.0;
        //selected.resize(node_num,0);
        /*
        sim_mat_sum = vector<float>(nodenum,0);
        for(int it=0;it < nodenum;it++){
            for(int itv=0;itv < nodenum;itv++)
                sim_mat_sum[it] += sim_mat[it][itv];
        }*/
    }
    vector<int> Set;
    double S_cost;
    double S_revenue;
    //vector<int> selected;
    //int node_num;
    //vector<vector<float>> sim_mat;
    //vector<float> sim_mat_sum;
    double f_S()
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
    double marginal(int e)
    {
        //Set.push_back(node);
        //selected[node]=1;
        //double f_S_and_u=f_S();
        //Set.pop_back();
        //selected[node]=0;
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
    /*
    double Sij_marginal(int e,int max_node_index)
    {
        //Set.push_back(node);
        //selected[node]=1;
        //double f_S_and_u=f_S();
        //Set.pop_back();
        //selected[node]=0;
        float M1=0,M2=0;
        M1 += sim_mat_sum[e];
        for(int it=0;it <= max_node_index;it++)
        {
            M2 += sim_mat[e][Set[it]];
            M2 += sim_mat[Set[it]][e];
        }
        M2 += sim_mat[e][e];
        return M1-M2;
    }*/
    double S_sub_u(int e)
    {
        double sum_value = 0.0;
        float M1 = 0, M2 = 0;
        for (int it = 0; it < Set.size(); it++)
        {
            if (Set[it] == e)
                continue;
            M1 += sim_mat_sum[Set[it]];
            for (int its = 0; its < Set.size(); its++)
            {
                if (Set[its] == e)
                    continue;
                M2 += sim_mat[Set[it]][Set[its]];
            }
        }
        sum_value = M1 - M2;
        return sum_value;
    }
    void copy(const S_gamma &temp)
    {
        S_revenue = temp.S_revenue;
        S_cost = temp.S_cost;
        Set.assign(temp.Set.begin(), temp.Set.end());
        //S* need not maintain selected and sim_mat to calculate f(S)
        //selected.assign(temp.selected.begin(),temp.selected.end());
    }
};
class S_pair
{
public:
    vector<S_gamma> spair;
    double gamma;
    S_pair() {}
    S_pair(double g)
    {
        gamma = g;
        spair.push_back(S_gamma());
        spair.push_back(S_gamma());
    }
};
//double f_u(int node,int node_num,dataset *data);
//S_gamma OnePass(dataset* data, double B,double eps,long long int &oracle_times);

double f_u(int node, int node_num, dataset *data)
{
    double sum_value = 0.0;
    for (int i = 0; i < node_num; i++)
    {
        sum_value += sim_mat[node][i];
    }
    return sum_value;
}

S_gamma OnePass(dataset *data, double B, double eps, long long int &oracle_times, double &zep_rev, long long int &zep_ora)
{

    cout << "OnePass ---------start--------- " << endl;

    oracle_times = 0;
    S_gamma S_best;
    vector<S_pair> S_array;
    double max_gamma = 0.0;        //max gamma in S_array, help to update S_array
    int min_gamma_index_in_C = -1; //-1 means S_array is empty, else means min gamma index which is still available now
    int max_gamma_index_in_C = -1;
    for (int u = 0; u < data->index.size(); u++)
    {
        //generate C
        if (data->rate_costs[u] > B)
            continue;

        oracle_times++;

        double fu_value = f_u(u, data->index.size(), data);
        if (fu_value > S_best.S_revenue)
        {
            S_best.Set.clear();
            S_best.Set.push_back(u);
            S_best.S_revenue = fu_value;
            S_best.S_cost = data->rate_costs[u];
            //S* need not maintain selected to calculate f(S)

            //fill(S_best.selected.begin(), S_best.selected.end(), 0);
            //S_best.selected[u]=1;
        }
        if (fu_value <= 0)
            continue;
        double left_temp = S_best.S_revenue / (4.0 * B * (1 + eps));
        double right_temp = fu_value / data->rate_costs[u];
        if (left_temp > right_temp)
            continue;

        double left = ceil(log(left_temp) / log(1 + eps));
        double right = floor(log(right_temp) / log(1 + eps));

        //vector<double> C;
        if (min_gamma_index_in_C == -1) //first visit S_array
        {
            for (int t = left; t <= right; t++)
            {
                double gamma_temp = pow(1 + eps, t);
                S_array.push_back(S_pair(gamma_temp));
            }
            min_gamma_index_in_C = 0;                  //the index of min gamma which is still available now
            max_gamma = pow(1 + eps, right);           //now max gamma
            max_gamma_index_in_C = S_array.size() - 1; //max gamma index in C now
        }
        else //not first visit
        {
            double now_min_gamma_in_C = pow(1 + eps, left);  //min gamma in C now
            double now_max_gamma_in_C = pow(1 + eps, right); //max gamma in C now
            if (now_min_gamma_in_C > max_gamma)              //all old S_arrary shouble be removed
            {
                min_gamma_index_in_C = S_array.size(); //the index of min gamma which is still available now
                for (int t = left; t <= right; t++)    //then add new S gamma pair
                {
                    double gamma_temp = pow(1 + eps, t);
                    S_array.push_back(S_pair(gamma_temp));
                }
                max_gamma_index_in_C = S_array.size() - 1;     //max gamma index in C now
                max_gamma = S_array[S_array.size() - 1].gamma; //the last element always is the max_gamma anyway
            }
            else //else find where is the min gamma index now, is equivalent to remove all S gamma < left
            {
                bool need_update = true; //judge S_array shouble be updated or not
                for (int iter = min_gamma_index_in_C; iter < S_array.size(); iter++)
                {
                    if (S_array[iter].gamma < now_min_gamma_in_C)
                        min_gamma_index_in_C++;
                    if (S_array[iter].gamma >= now_max_gamma_in_C)
                    {
                        max_gamma_index_in_C = iter;
                        need_update = false;
                        break;
                    }
                }
                if (need_update)
                {
                    //finally, go through all gamma now, put new S gamma pair if needed
                    for (int t = left; t <= right; t++)
                    {
                        double gamma_temp = pow(1 + eps, t);
                        if (gamma_temp > max_gamma)
                            S_array.push_back(S_pair(gamma_temp));
                    }
                    max_gamma = S_array[S_array.size() - 1].gamma; //the last element always is the max_gamma anyway
                    max_gamma_index_in_C = S_array.size() - 1;
                }
            }
        }
        //foreach gamma\in C
        for (int iter = min_gamma_index_in_C; iter <= max_gamma_index_in_C; iter++)
        {
            oracle_times++;
            oracle_times++;

            double marginal1 = S_array[iter].spair[0].marginal(u);
            double marginal2 = S_array[iter].spair[1].marginal(u);
            if (marginal1 > marginal2)
            {
                if ((marginal1 / data->rate_costs[u] >= S_array[iter].gamma) && (S_array[iter].spair[0].S_cost + data->rate_costs[u] <= B))
                {
                    S_array[iter].spair[0].S_cost += data->rate_costs[u];
                    S_array[iter].spair[0].S_revenue += marginal1;
                    //S_array[iter].spair[0].selected[u]=1;
                    S_array[iter].spair[0].Set.push_back(u);

                    if (S_array[iter].spair[0].S_revenue > S_best.S_revenue)
                    {
                        S_best.copy(S_array[iter].spair[0]);
                    }
                }
            }
            else
            {
                if ((marginal2 / data->rate_costs[u] >= S_array[iter].gamma) && (S_array[iter].spair[1].S_cost + data->rate_costs[u] <= B))
                {
                    S_array[iter].spair[1].S_cost += data->rate_costs[u];
                    S_array[iter].spair[1].S_revenue += marginal2;
                    //S_array[iter].spair[1].selected[u]=1;
                    S_array[iter].spair[1].Set.push_back(u);

                    if (S_array[iter].spair[1].S_revenue > S_best.S_revenue)
                    {
                        S_best.copy(S_array[iter].spair[1]);
                    }
                }
            }
        }
    }
#ifndef ZEP
    cout << "OnePass & Budget: " << B << endl;
    cout << "S*:" << endl;
    cout << "  revenue: " << S_best.S_revenue << " cost: " << S_best.S_cost << " size: " << S_best.Set.size() << endl;
    for (int i = 0; i < S_best.Set.size(); i++)
        cout << S_best.Set[i] << " ";
    cout << endl;
    cout << "oracle times: " << oracle_times << endl;

    cout << "OnePass ---------end--------- " << endl
         << endl;
#endif

    zep_ora = oracle_times;
    zep_rev = S_best.S_revenue;

    return S_best;
}

#endif //STREAMING_ALGORITHM_ONEPASS_H
