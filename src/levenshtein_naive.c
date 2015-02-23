#include <stdio.h>
#include <string.h>

int levenshtein(char *s, int ls, char *t, int lt) {
    // if one of the lengths are equal to 0, return the other length
    if(!ls) return lt;
    if(!lt) return ls;

    // if last character of the two string are identical, no added distance
    if(s[ls] == t[lt])
        return levenshtein(s, ls-1, t, lt-1);

    // tree of recursive calls; no reuse of already calculated distances
    int a = levenshtein(s, ls-1, t, lt)   + 1;
    int b = levenshtein(s, ls,   t, lt-1) + 1;
    int c = levenshtein(s, ls-1, t, lt-1) + 1;

    // min{a,b,c}
    if(a < b && a < c)
        return a;
    else if(b < a && b < c)
        return b;
    else return c;
}

int main(int argc, char *argv[])
{
    char *a = "aaaaaaaaaaaaaa";
    char *b = "abcdefghijklmn";

    int dist = levenshtein(a, strlen(a), b, strlen(b));

    printf("Levenshtein distance between \"%s\" and \"%s\": %d\n", a, b, dist);

    return 0;
}
