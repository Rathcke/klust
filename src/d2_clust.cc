#include <iostream>
#include <cstring>
#include <algorithm>
#include <cmath>
using namespace std;

int gram_pos(string s) {
	int slen = s.length();
	int cost = 0;
	for (int i = slen - 1; i >= 0; i--) {
		switch (s[i]) {
			case 'a':
				break;
			case 'c':
				cost += pow(4,slen-i-1);
				break;
			case 'g':
				cost += 2*pow(4,slen-i-1);
				break;
			case 't':
			case 'u':
				cost += 3*pow(4,slen-i-1);
				break;
			default:
				cout << "Unknown char passed to gram_pos";
		}
	}
	return cost;
}

int d2_clust(string s, string t, int k) {
	int gramlen = pow(4,k);
	int slen = s.length();
	int tlen = t.length();
	int dist = 0;
	int v0[gramlen];
	int v1[gramlen];
	int i;
	for (i = 0; i < gramlen; i++) {
		v0[i] = 0;
		v1[i] = 0;
	}
	for (i = 0; i <= slen-k; i++) {
		v0[gram_pos(s.substr(i,k))]++;
	}
	for (i = 0; i <= tlen-k; i++) {
		v1[gram_pos(t.substr(i,k))]++;
	}
	for (i = 0; i < gramlen; i++) {
		dist += pow(v0[i]-v1[i],2);
	}
	return sqrt(dist);
}

int main() {
	
}