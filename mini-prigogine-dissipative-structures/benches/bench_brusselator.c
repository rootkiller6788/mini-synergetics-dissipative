#include "pds_brusselator.h"
#include <stdio.h>
#include <time.h>
int main(void){clock_t s=clock();PDSBrusselatorParams p={1,3,1,1,1,1,1,1};PDSBrusselatorState i={0.5,1};for(int k=0;k<50;k++){PDSBrusselatorResult*r=pds_brusselator_simulate(&i,&p,10,0.01,20);pds_brusselator_result_free(r);}printf("Time: %.3f s
",(double)(clock()-s)/CLOCKS_PER_SEC);return 0;}
