#include "pds_core.h"
#include "pds_thermo.h"
#include "pds_brusselator.h"
#include <stdio.h>
int main(void){puts("=== Prigogine Dissipative Structures Demo ===");puts("1. Creating dissipative system (Brusselator)...");PDSDissipativeSystem*sys=pds_system_create("demo",PDS_BRUSSELATOR,2,1,2);printf("System created: %s, species=%d
",sys->name,sys->n_species);puts("2. Simulating Brusselator beyond Hopf bifurcation...");PDSBrusselatorParams par={1.0,3.0,1,1,1,1,1,1};PDSBrusselatorState init={0.5,1.0};PDSBrusselatorResult*r=pds_brusselator_simulate(&init,&par,20,0.01,20);printf("Limit cycle detected: %s
",r->has_limit_cycle?"YES":"no");if(r->has_limit_cycle)printf("Period: %.2f, Amplitude: %.2f
",r->oscillation_period,r->oscillation_amplitude);pds_brusselator_result_free(r);pds_system_free(sys);puts("Demo complete.");return 0;}
