#include "pds_brusselator.h"
#include <stdio.h>
#include <stdlib.h>

int main(void) {
    PDSBrusselatorParams params = {1.0, 3.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
    PDSBrusselatorState initial = {0.5, 1.0};
    double bc = pds_brusselator_bifurcation_point(&params);
    printf("Brusselator: a=%.1f, b=%.1f, b_c=%.2f", params.a, params.b, bc);
    putchar(0x0a);
    puts("b > b_c: expecting dissipative structure (limit cycle)");
    PDSBrusselatorResult* result = pds_brusselator_simulate(
        &initial, &params, 30.0, 0.01, 50);
    pds_brusselator_print_result(result);
    pds_brusselator_result_free(result);
    return 0;
}
