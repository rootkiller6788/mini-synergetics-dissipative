#include "pds_core.h"
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <stdlib.h>

int main(void) {
    PDSDissipativeSystem* sys = pds_system_create("test", PDS_BRUSSELATOR, 2, 2, 2);
    assert(sys != NULL);
    pds_system_set_parameter(sys, 0, 1.0);
    pds_system_set_parameter(sys, 1, 2.5);
    assert(fabs(pds_system_get_parameter(sys, 0) - 1.0) < 1e-10);
    double x0[] = {0.5, 1.0};
    pds_system_set_state(sys, x0);
    double x_out[2];
    pds_system_get_state(sys, x_out);
    assert(fabs(x_out[0] - 0.5) < 1e-10);
    for (int i = 0; i < 10; i++) pds_system_record_state(sys);
    assert(sys->history_count == 10);
    PDSVector* v1 = pds_vector_create(3);
    pds_vector_set(v1, 0, 1.0); pds_vector_set(v1, 1, 2.0); pds_vector_set(v1, 2, 3.0);
    double norm = pds_vector_norm(v1);
    assert(fabs(norm - sqrt(14.0)) < 1e-10);
    PDSMatrix* A = pds_matrix_create(2, 2);
    pds_matrix_set(A, 0, 0, 1.0); pds_matrix_set(A, 0, 1, 2.0);
    pds_matrix_set(A, 1, 0, 3.0); pds_matrix_set(A, 1, 1, 4.0);
    assert(fabs(pds_matrix_determinant_2x2(A) - (-2.0)) < 1e-10);
    assert(pds_sign(5.0) == 1.0);
    assert(pds_sign(-3.0) == -1.0);
    assert(pds_is_near(0.1, 0.1000001, 0.001));
    assert(pds_vector_create(-1) == NULL);
    pds_vector_free(v1);
    pds_matrix_free(A);
    pds_system_free(sys);
    puts("ALL CORE TESTS PASSED");
    return 0;
}
