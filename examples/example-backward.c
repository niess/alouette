/* Standard library includes. */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/* ALOUETTE: a TAUOLA BMC wrapper. */
#include "alouette.h"

/* Tau longitudinal polariser. */
static void polarise(int pid, const double momentum[3], double * polarisation)
{
        const double polar = 1.;
        double nrm = momentum[0] * momentum[0] + momentum[1] * momentum[1] +
            momentum[2] * momentum[2];
        if (nrm <= 0.) {
                polarisation[0] = 0.;
                polarisation[1] = 0.;
                polarisation[2] = 0.;
                return;
        }
        nrm = polar / sqrt(nrm);
        polarisation[0] = nrm * momentum[0];
        polarisation[1] = nrm * momentum[1];
        polarisation[2] = nrm * momentum[2];
}

int main()
{
        /* Randomise a few tau decays. */
        const int mode = 0;
        const int pid = 16;
        const double momentum[3] = { 0., 0., 1. };

        int i;
        for (i = 0; i < 3; i++) {
                struct alouette_products products;
                if (alouette_undecay(
                    mode, pid, momentum, polarise, 0., &products) !=
                    ALOUETTE_RETURN_SUCCESS)
                        continue;

                printf("# Event %d (%.5E) :\n", i + 1, products.weight);
                int j;
                for (j = 0; j < products.size; j++) {
                        printf("    %4d %12.5lE %12.5lE %12.5lE\n",
                            products.pid[j], products.P[j][0], products.P[j][1],
                            products.P[j][2]);
                }
        }

        exit(EXIT_SUCCESS);
}
