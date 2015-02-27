/**
 * A memoized version of the naive implementation of the Levenshtein distance.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int levenshtein(const char *s, const char *t) {
    int s_len = strlen(s), t_len = strlen(t);

    int *d = malloc(sizeof(int) * s_len * t_len);  // distance matrix

    int i,j;
    for (i = 0; i <= s_len; i++)                   // set all entries to -1
        for (j = 0; j <= t_len; j++)
            d[sizeof(int) * i + sizeof(int) * j] = -1;


    int dist(int ls, int lt) {
        // see if already calculated
        if(d[sizeof(int) * ls + sizeof(int) * lt] >= 0)
            return d[sizeof(int) * ls + sizeof(int) * lt];

        int ret; // variable for holding distance return value

        // if one of the lengths are equal to 0, return the other length
        if (!ls)
            ret = lt;
        else if (!lt)
            ret = ls;
        // if last character of the two strings are identical, no added distance
        else if (s[ls] == t[lt])
            ret = dist(ls-1, lt-1);
        else {
            // tree of recursive calls
            int a = dist(ls-1, lt)   + 1;
            int b = dist(ls,   lt-1) + 1;
            int c = dist(ls-1, lt-1) + 1;

            // min{a,b,c}
            if(a < b && a < c)
                ret = a;
            else if(b < a && b < c)
                ret = b;
            else
                ret = c;
        }
        return d[sizeof(int) * ls + sizeof(int) * lt] = ret;
    }

    int distance = dist(s_len, t_len);
    free(d);
    return distance;
}

int main(int argc, char *argv[])
{

    char *a = "caaccugguugauccugccaguagucauaugcuugucugcuaaugggagacacugcugucuaaguauaaacaccuuuauacacgugaaacugcgaauggcucauuaaaacaguuauaguuccuuugagaacacugcuacuuggauaaccguaguaauucuagagcuaauacaugccuaaacgcccgacucacggagggugguauuuauuggauaaaaaaccaucgugcccucucgggggcacuuugagaugauucacaauaacuagucggaucgcauagcuuugugcuggcgacuccucauucaaauuucugcccuaucaacuuucgaugguagaguauuggucuaccauggugucgacgggugacggggaauuaggguucgauuccggagagggagccugagaaacggcuaccacauccaaggaaggcagcaggcgcgcaaauuacccaaucccgauucggggagguagugacaaaaaauaacaauagggggcccuuugggucuucuaauuggaaugagaacaauuuaaaucccuuaucgaggauccauuggagggcaagucuggugccagcagccgcgguaauuccagcuccaauagcguauauuaaaguuguugcaguuaaaacgcucguagucggaccucggagcacccccacgggucgggcuuuugccuaugugcucgugggcgguuugucuccuuuugucgggacgcggugcgcggggacuuuacugucucugugcguuauaccugcgggccgaccguuuacugugaagaaaguaaaguguucaaagcaggccauugccuugaauaugugagcauggaauaauagaauaggacuugggcucuauuuuguugguuuccagugaccaaguaaugauuaauagggacgguugggggcauucguauuucauugucagaggugaaauucuuggauugauggaagacgaacaacugcgaaagcaucugccauggauguuuucauugaucaagaacgaaaguuaggggaucgaagacgaucagauaccgucguagucuuaaccauaaacgaugccgacuggggauugacggggguuuuuugaacgacucugucagcacccugagggaaaccaaagucuuuggguucuggggggaguauggucgcaaggcugaaacuuaaaggaauugacggaagggcaccaccaggaguggagccugcggcuuaauuugacucaacacgggaaaacuuaccagguccggacagaagaaugauugacagacugaaagcucuuucuugauuuuuugguugguggugcauggccguucuuaguugguggagugauuugucugguuaauuccguuaacgaacgagaccucguccugcuaaauaggugcgcgcaugcuaguacugcguucuuaccuucuuagagggacuaugcgcgucuagcguauggaagauugaggcaauaacaggucugugaugcccuuagauguucugggccgcacgcgcgcuacacugaugcauucaacgaguuucacuugauagcccugggucggaaggccuggguaaucuuuugaaagugcaucgugcuggggacagaucauugcaauuauugaucucaaacgaggaauuccuuguaggcgcaggucaucagccugcgccgaauacgucccugcccuuuguacacaccgcccgucgcuccuaccgauugaaugguccgaugaaauguccggacggcggcucucaggacggacuugcgccagaaugcuugucguggaaagcucauuaaaucuuaucauuuagaggaaggagaagucguaacaagguuuccguaggugaaccugcagaaggauc";
    char *b = "ggcucaggaugaacgcuggcggcaugcuuaacacaugcaagucguacgcaugcaauuuggcuugccagauugcgaugaguggcggacgggugaguaacacguaagaaccuaccuuuuggagagggauaaccauuggaaacgauggcuaauaccucguauugcugagaagugaaagaugaaaaucgccaauagaugggcuugcggcugauuagcuuguuggugagguaauggcuuaccaaggcaaugaucaguagcuggucugagaggaugaucagccacacugggacugagacacggcccagacuccuacgggaggcagcagugaggaauuuuccgcaaugggcgacagccugacggagcaaugccgcgugaaggaugaaggccuauggguuguaaacuucuuuucucagagaagaagcauugacgguaucugaggaauaagcaucggcuaacucugugccagcagccgcgguaagacagaggaugcaagcguuauccggaaugauugggcguaaagcgucuguagguggcuuaaaaagucuccugucagagaucagggcuuaacccugggccggcaggagaaacucuuaggcuagaguuugguaggggcagagggaauucccgguggagcggugaaaugcguagagaucgggaggaacaccaaaggcgaaagcacucugcugggccauaacugacacugagagacgaaagcgaggggagcaaaagggauuagauaccccuguaguccucgccguaaacgauggauacuagauguuggguagguuaaaucacucaguaucguagcuaacgcgugaaguaucccgccuggggaguaugcucgcaagagugaaacucaaaggaauugacgggggcccgcacaagcgguggagcaugugguuuaauucgaugcaacgcgaagaaccuuaccaggacuugacaugccacuuuuucccugaaaggggaaguuccagaguggacacagguggugcauggcugucgucagcucgugucuugagauguuggguuaagucccgcaacgagcgcaacccuuguuuugaauugccaguaaugggaaauucaaaagacugccggugacaagccggaggaaggugaggaugacgucaagucagcaugccccuuacguccugggcgacacacgugcuacaauggccgggacaaagagaugcaaacccgcgagggcuagccaaccucaaaaacccggucucaguucggauugcaggcugcaacucgccugcaugaagucggaaucgcuaguaaucgcaggucagccauacugcggugaauacguucccgggccuuguacacaccgcccgucacaccaugggagcuggcuaugcccaaagucguuaccccaaccuuuuaggagggggacgccuaaggcagagcuagugacuagggugaagucguaacaag";


    //char *a = "aaaaaaaaaaaaaaaaaaaaaaaaaa";
    //char *b = "abcdefghijklmnopqrstuvwxyz";

    int dist = levenshtein(a, b);

    printf("Levenshtein distance between \"%s\" and \"%s\": %d\n", a, b, dist);

    return 0;
}
