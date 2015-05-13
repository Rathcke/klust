#include <algorithm>
#include <cstdlib>
#include <vector>
#include <iomanip>

#include "Seq.h"
#include "Utils.h"
#include "Distance.h"

namespace Utils {

using namespace std;

/**
 * Return random int in interval [a,b].
 */
int get_rand(int a, int b) {
    random_device rd;
    mt19937 mt(rd());
    uniform_int_distribution<int> dist(a, b);

    return dist(mt);
}

char get_rand_base_not(char c) {
    char not_a[] = "CGT";
    char not_c[] = "AGT";
    char not_g[] = "ACT";
    char not_t[] = "ACG";

    int r = get_rand(0,2);

    switch (c) {
        case 'A' : return not_a[r];
        case 'C' : return not_c[r];
        case 'G' : return not_g[r];
        case 'T' : return not_t[r];
    }
    return 'A';
}

void permute(vector<Seq>& seqs, int count, double ratio, ofstream& fs_cts) {

    for (auto& s : seqs) {
        // write original sequence
        fs_cts << '>' << s.desc << '\n'
               << s.to_string() << '\n';

        for (int i = 0; i < count; ++i) {

            int rand_indices_count = s.length * ratio;

            vector<int> rand_indices;   // indices in sequence to be changed

            int j = 0;
            while (j < rand_indices_count) {
                int random_index = get_rand(0, s.length-1);
                if (find(rand_indices.begin(),rand_indices.end(), random_index)
                        == rand_indices.end()) {
                    rand_indices.push_back(random_index);
                    ++j;
                }
            }

            string s_perm = s.to_string();

            for (int r : rand_indices)
                s_perm[r] = get_rand_base_not(s_perm[r]);

            fs_cts << '>' << s.desc << '\n'
                   << s_perm << '\n';

        }
    }
}

void permute_chunks(vector<Seq>& seqs, int count, double ratio, ofstream& fs_cts) {


    int chunk_size = 5;

    for (auto& s : seqs) {
        // write original sequence
        fs_cts << '>' << s.desc << '\n'
               << s.to_string() << '\n';

        for (int i = 0; i < count; ++i) {

            int rand_indices_count = s.length * ratio;

            vector<int> rand_indices;   // indices in sequence to be changed

            int j = 0;
            while (j < rand_indices_count) {
                int random_index = get_rand(0, s.length-1);
                if (find(rand_indices.begin(),rand_indices.end(), random_index)
                        == rand_indices.end()) {
                    rand_indices.push_back(random_index);
                    ++j;
                }
            }

            string s_perm = s.to_string();

            for (int r : rand_indices)
                for (int j = 0; j < chunk_size; ++j)
                    s_perm[r] = get_rand_base_not(s_perm[min(r+j, (int) s.length-1)]);

            fs_cts << '>' << s.desc << '\n'
                   << s_perm << '\n';

        }
    }
}


void print_matrix(vector<Seq>& seqs, ostream& fs_mat, Distance& dist) {

    double matrix[seqs.size()][seqs.size()];

    for (unsigned int i = 0; i < seqs.size(); ++i) {
        for (unsigned int j = 0; j < seqs.size(); ++j) {

            if (i == j)
                matrix[i][j] = 1;

            matrix[i][j] = dist.distance(seqs[i], seqs[j]);

        }
    }
    for (unsigned int i = 0; i < seqs.size(); ++i) {
        for (unsigned int j = 0; j < seqs.size(); ++j) {

            fs_mat << setw(6) << setprecision(3) << matrix[i][j];

        }
        fs_mat << '\n';
    }
}

} // namespace Util
