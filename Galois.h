#include <cmath>
#include <stdlib.h>

int GFOrder(unsigned char * a){
	int highest = 0;
	for (int i = 0; i<15; i++) {
		highest = (a[i] == 1 ? i : highest);
	}
	return (highest);
}

unsigned char GFMul(unsigned char a, unsigned char b){
	unsigned char c = 0;
	unsigned char bits[16];
	unsigned char mask = 1;

	// initialize the bit sequence
	for (int i = 0; i<15; i++) {
		bits[i] = 0;
	}
	// multiply in gf
	for (int i = 0; i<8; i++) {
		for (int j = 0; j<8; j++) {
			bits[i + j] = bits[i + j] ^ (((a >> i) & mask) & ((b >> j) & mask));
		}
	}
	// find modulo x8 + x4 + x3 + x + 1 =
	int order = GFOrder(bits) - 8;
	while (order >= 0) {
		bits[order + 8] = bits[order + 8] ^ 1;
		bits[order + 4] = bits[order + 4] ^ 1;
		bits[order + 3] = bits[order + 3] ^ 1;
		bits[order + 1] = bits[order + 1] ^ 1;
		bits[order] = bits[order] ^ 1;
		order = GFOrder(bits) - 8;
	}
	// initialize the bit sequence
	for (int i = 0; i<15; i++) {
		//cout << (int) bits[i] << " " ;
		c = c + bits[i] * pow(2, i);
	}
	return (c);
}

void GFInitInverse(unsigned char * inverses) {
	unsigned char c;
	for (int i = 1; i<256; i++) {
		for (int j = 1; j<256; j++) {
			c = GFMul(i, j);
			if (c == 1) {
				inverses[i] = j;
			}
		}
	}

}

unsigned char GFPower(unsigned char a, int exp){
	unsigned char c = 1;
	for (int i = 1; i <= exp; i++){
		c = GFMul(c, a);
	}
	return (c);
}

unsigned char GFAdd(unsigned char a, unsigned char b){
	unsigned char c;
	c = a ^ b;
	return (c);
}
//------------- Galois Field Operations END ---------------------------