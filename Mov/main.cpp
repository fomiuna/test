#include "offline.h"
#include "onepass.h"
#include "fantom.h"
#include "samplegreedy.h"
#include "multipass.h"
#include <stdlib.h>
#include <iomanip>

int main(int argc, char *argv[])
{
    double offline_rev, onepass_rev, multi_rev, sample_rev, fantom_rev;
    long long int offline_ora, onepass_ora, multi_ora, sample_ora, fantom_ora;

    double B = atof(argv[1]);
    double eps = atof(argv[2]);
    // cout << "eps: " << eps << endl;

    string filename = "movies_1793.txt";
    dataset *da = new dataset(filename);

    double default_p = sqrt(2) - 1;
    Offline(da, B, offline_rev, offline_ora);
    long long int onepass_oracle = 0;
    OnePass(da, B, eps, onepass_oracle, onepass_rev, onepass_ora);
    MultiPass(da, B, eps, multi_rev, multi_ora);
    iter_sample(B, default_p, da, sample_rev, sample_ora);
    fantom(1, B, eps, da, fantom_rev, fantom_ora);

    cout << "\n\n";
    cout << "eps = " << eps << "   budget = " << B << endl;
    cout << "------------------------------------------" << endl;
    cout << setiosflags(ios::left) << setw(14) << "Methods" << resetiosflags(ios::left)
         << setiosflags(ios::right) << setw(11) << "Revenue" << setw(16) << "Oracles"
         << resetiosflags(ios::right) << endl;
    cout << "------------------------------------------" << endl;
    string methods[] = {"Offline", "OnePass", "MultiPass", "SampleGreedy", "Fantom"};
    double rev[] = {offline_rev, onepass_rev, multi_rev, sample_rev, fantom_rev};
    long long int ora[] = {offline_ora, onepass_ora, multi_ora, sample_ora, fantom_ora};
    for (int i = 0; i < 5; ++i)
    {
        cout << setiosflags(ios::left) << setw(14) << methods[i] << resetiosflags(ios::left)
             << setiosflags(ios::right) << setw(11) << rev[i] << setw(16) << ora[i]
             << resetiosflags(ios::right) << endl;
    }
    cout << "------------------------------------------" << endl;
    cout << "\n\n";

    return 0;
}
