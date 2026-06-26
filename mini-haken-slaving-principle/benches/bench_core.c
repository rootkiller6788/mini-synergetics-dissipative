/* Benchmark: eigenvalue computation and adiabatic elimination performance */
#include "haken_core.h"
#include "haken_mode_decomp.h"
#include "haken_adiabatic.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int main(void) {
    printf("=== Haken Slaving Principle — Benchmarks ===\n\n");
    clock_t start, end;

    /* Benchmark 1: Eigenvalue computation scaling */
    printf("Benchmark 1: QR eigenvalue computation\n");
    int sizes[] = {4, 8, 16, 32};
    for (int si = 0; si < 4; si++) {
        int n = sizes[si];
        Haken_System* sys = haken_system_create(n, NULL, NULL);
        for (int i = 0; i < n; i++) {
            sys->linear_matrix[i * n + i] = -1.0 - (double)i;
            if (i > 0) sys->linear_matrix[i * n + (i - 1)] = 0.1;
            if (i < n-1) sys->linear_matrix[(i + 1) * n + i] = 0.1;
        }
        start = clock();
        haken_compute_eigenmodes(sys);
        end = clock();
        double elapsed = (double)(end - start) / CLOCKS_PER_SEC;
        printf("  n=%2d: %.4f seconds\n", n, elapsed);
        haken_system_free(sys);
    }

    /* Benchmark 2: Adiabatic elimination */
    printf("\nBenchmark 2: Adiabatic elimination (n=8, 2 critical, 6 stable)\n");
    Haken_System* sys = haken_system_create(8, NULL, NULL);
    for (int i = 0; i < 8; i++) {
        sys->linear_matrix[i * 8 + i] = (i < 2) ? -0.001 * (i+1) : -1.0 - i;
    }
    haken_compute_eigenmodes(sys);
    haken_classify_modes(sys, 0.1, 0.5);

    start = clock();
    Haken_SlavingResult result;
    Haken_AdiabaticConfig config = {HAKEN_ADIABATIC_LINEAR, 1e-8, 50, true, 0.1};
    for (int trial = 0; trial < 100; trial++) {
        haken_adiabatic_eliminate(sys, &config, &result);
    }
    end = clock();
    printf("  100 iterations: %.4f seconds (avg %.4f ms/iter)\n",
           (double)(end - start) / CLOCKS_PER_SEC,
           (double)(end - start) / CLOCKS_PER_SEC * 10.0);

    haken_system_free(sys);

    printf("\n=== Benchmarks Complete ===\n");
    return 0;
}
