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
        /* Initialise the ALOUETTE library. */
        enum alouette_return rc;
        if ((rc = alouette_initialise(NULL, NULL)) != ALOUETTE_RETURN_SUCCESS) {
                fprintf(stderr, "alouette: %s\n", alouette_strerror(rc));
                exit(EXIT_FAILURE);
        };

        /* Randomise a few tau decays. */
        int i;
        for (i = 0; i < 3; i++) {
                double weight, momentum[3] = { 0., 0., 1. };
                if (alouette_undecay(16, momentum, polarise, 0., &weight) !=
                    ALOUETTE_RETURN_SUCCESS)
                        continue;

                int pid;
                printf("# Event %d (%.5E) :\n", i + 1, weight);
                while (alouette_product(&pid, momentum) ==
                    ALOUETTE_RETURN_SUCCESS) {
                        printf("    %4d %12.5lE %12.5lE %12.5lE\n", pid,
                            momentum[0], momentum[1], momentum[2]);
                }
        }

        exit(EXIT_SUCCESS);
}
