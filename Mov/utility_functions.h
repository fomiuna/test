
#ifndef STREAMING_ALGORITHM_UTILITY_FUNCTIONS_H
#define STREAMING_ALGORITHM_UTILITY_FUNCTIONS_H
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

#define ZEP

using namespace std;

vector<float> sim_mat_sum;
vector<vector<float>> sim_mat;
class dataset
{
public:
    vector<int> index;
    set<int> v;
    vector<vector<float>> matrix;
    vector<vector<string>> genres;
    vector<string> titles;
    float lambda;
    vector<int> year;
    vector<float> rating;
    vector<float> year_costs;
    vector<float> rate_costs;
    //vector<vector<float>> sim_mat;
    dataset(string);
    void cost_normalize();
    void readfile(string filename);
    vector<float> feature_extractor(string features, char delim);
    vector<string> genre_extractor(string genres, char delim);
    void sim_mat_generator();
};
dataset::dataset(string filename)
{
    readfile(filename);
    cost_normalize();
    sim_mat_generator();

    int node_num = index.size();
    sim_mat_sum = vector<float>(node_num, 0.0);
    for (int it = 0; it < node_num; it++)
    {
        for (int itv = 0; itv < node_num; itv++)
            sim_mat_sum[it] += sim_mat[it][itv];
    }
}

void dataset::cost_normalize()
{
    float rate = rating.size() / 10.0;
    float rate_coef, year_coef;
    float rate_sum, year_sum;
    float rate_base = 10.0, year_base = 1985;
    for (int i = 0; i < rating.size(); i++)
    {
        rate_sum += rating.at(i);
        year_sum += abs(year.at(i) - year_base);
    }
    rate_sum = rating.size() * rate_base - rate_sum;
    //year_sum = rating.size()*year_base;
    rate_coef = rate / rate_sum;
    year_coef = rate / year_sum;

    for (int i = 0; i < rating.size(); i++)
    {
        rate_costs.push_back(fabs(rating.at(i) - rate_base) * rate_coef);
        year_costs.push_back(fabs(year.at(i) - year_base) * year_coef);
    }
}

void dataset::readfile(string filename)
{
    //dataset* data = new dataset();
    ifstream infofile;
    infofile.open(filename.c_str(), ios::in);
    int count = 0;
    char delim = '|';
    if (infofile.is_open())
    {
        while (!infofile.eof())
        {
            std::string line;
            getline(infofile, line);
            if (line.empty())
                continue;
            std::string::size_type pos = line.find_first_of(delim);
            v.insert(count++);
            int prevpos = 0;
            int length = line.length();
            // id
            string str = line.substr(prevpos, pos - prevpos);
            index.push_back(atoi(str.c_str()));
            // titles
            prevpos = line.find_first_not_of(delim, pos);
            pos = line.find_first_of(delim, prevpos);
            str = line.substr(prevpos, pos - prevpos);
            titles.push_back(str);
            // genres
            prevpos = line.find_first_not_of(delim, pos);
            pos = line.find_first_of(delim, prevpos);
            str = line.substr(prevpos, pos - prevpos);
            genres.push_back(genre_extractor(str, ' '));
            // lambda
            lambda = 0.2;
            //features
            prevpos = line.find_first_not_of(delim, pos);
            pos = line.find_first_of(delim, prevpos);
            str = line.substr(prevpos, pos - prevpos);
            matrix.push_back(feature_extractor(str, ' '));

            // year
            prevpos = line.find_first_not_of(delim, pos);
            pos = line.find_first_of(delim, prevpos);
            str = line.substr(prevpos, pos - prevpos);
            year.push_back(atoi(str.c_str()));

            // rating
            prevpos = line.find_first_not_of(delim, pos);
            pos = line.find_first_of(delim, prevpos);
            str = line.substr(prevpos, pos - prevpos);
            rating.push_back(atof(str.c_str()));
        }
    }
}

void dataset::sim_mat_generator()
{ //计算Mij
    sim_mat = vector<vector<float>>(index.size(), vector<float>(index.size(), 0));
    for (int i = 0; i < index.size(); i++)
    {
        for (int j = i + 1; j < index.size(); j++)
        {
            float distance = 0;
            float sim_dis = 0;
            for (int k = 0; k < matrix[0].size(); k++)
            {
                distance += pow((matrix[i][k] - matrix[j][k]), 2);
            }
            distance = sqrt(distance);
            sim_dis = exp(-1 * lambda * distance);
            sim_mat[i][j] = sim_dis;
            sim_mat[j][i] = sim_dis;
        }
    }
}

vector<float> dataset::feature_extractor(string features, char delim)
{
    vector<float> str;
    std::string::size_type pos = features.find_first_of(delim);
    int prevpos = 0;
    while (string::npos != pos || string::npos != prevpos)
    {
        str.push_back(atof(features.substr(prevpos, pos - prevpos).c_str()));
        prevpos = features.find_first_not_of(delim, pos);
        pos = features.find_first_of(delim, prevpos);
    }
    return str;
}

vector<string> dataset::genre_extractor(string genres, char delim)
{
    vector<string> str;
    std::string::size_type pos = genres.find_first_of(delim);
    int prevpos = 0;
    while (string::npos != pos || string::npos != prevpos)
    {
        str.push_back(genres.substr(prevpos, pos - prevpos).c_str());
        prevpos = genres.find_first_not_of(delim, pos);
        pos = genres.find_first_of(delim, prevpos);
    }
    return str;
}

#endif //STREAMING_ALGORITHM_UTILITY_FUNCTIONS_H
