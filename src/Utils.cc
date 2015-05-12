#include <algorithm>
#include <cstdlib>
#include <vector>
#include <iomanip>

#include "Seq.h"
#include "Utils.h"
#include "Distance.h"

namespace Utils {

using namespace std;

char get_rand_base(char c) {
    char not_a[] = "CGT";
    char not_c[] = "AGT";
    char not_g[] = "ACT";
    char not_t[] = "ACG";

    switch (c) {
        case 'A' : return not_a[rand() % 3];
        case 'C' : return not_c[rand() % 3];
        case 'G' : return not_g[rand() % 3];
        case 'T' : return not_t[rand() % 3];
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
                int random_index = rand() % s.length;
                if (find(rand_indices.begin(),rand_indices.end(), random_index)
                        == rand_indices.end()) {
                    rand_indices.push_back(random_index);
                    ++j;
                }
            }

            string s_perm = s.to_string();

            for (int r : rand_indices)
                s_perm[r] = get_rand_base(s_perm[r]);

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
