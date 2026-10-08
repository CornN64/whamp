#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdint.h>
static void bits(double v){ uint64_t u; memcpy(&u,&v,8); printf("%016llX\n",(unsigned long long)u); }
int main(void){
    int i;
    for(i=1;i<=10;i++){
        double x = 0.1*i + 0.37;
        bits(exp(x));
        bits(log(x));
        bits(sin(x));
        bits(cos(x));
    }
    return 0;
}
