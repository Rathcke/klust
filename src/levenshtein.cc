#include <iostream>
#include <cstring>
#include <algorithm>
using namespace std;

int leven(string s, string t) {
	int slen = s.length();
	int tlen = t.length();
	if (slen == 0) {
		return tlen;
	}
	if (tlen == 0) {
		return slen;
	}
	int col[slen + 1];
	int pcol[slen + 1];
	for (int i = 0; i < slen+1; i++) {
		pcol[i] = i;
	}
	for (int i = 0; i < tlen; i++) {
		col[0] = i+1;
		for (int j = 0; j < slen; j++) {
			int cost = !(s[i] == t[j]);
			col[j+1] = min(col[j] + 1, min(pcol[j+1] + 1, pcol[j] + cost));
		}
		for (int j = 0; j < slen + 1; j++) {
			pcol[j] = col[j];
		}
	}
	return col[slen];
}

int main() {
	cout << leven("actgac", "acgga") << '\n';
}