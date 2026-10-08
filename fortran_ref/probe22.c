#include <stdio.h>
#include <complex.h>
#include <string.h>
#include <stdint.h>
static void bits(double v){ uint64_t u; memcpy(&u,&v,8); printf("%016llX\n",(unsigned long long)u); }
int main(void){
    int k;
    for(k=1;k<=8;k++){
        double a=1.0+ k*0.13, b=2.0- k*0.07, c=3.0+ k*0.29, d=0.5- k*0.11;
        double _Complex z = (a+b*I);
        double _Complex w = (c+d*I);
        double _Complex m = z*w;
        bits(creal(m)); bits(cimag(m));
        double _Complex q = z/w;
        bits(creal(q)); bits(cimag(q));
    }
    return 0;
}
