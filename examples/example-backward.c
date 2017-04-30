/* Standard library includes. */
#include <math.h>
#include <stdio.h>

/* TAUOLA C wrapper. */
#include "alouette.h"

/* Tau longitudinal polariser. */
static void polarise(
    int pid, const double momentum[3], double * polarisation)
{
        const double polar = 1.;
        double nrm = momentum[0] * momentum[0] +
            momentum[1] * momentum[1] + momentum[2] * momentum[2];
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
        /* Initialise TAUOLA through its wrapper. */
        alouette_initialise(1, NULL);

        /* Randomise a few tau decays. */
        int i;
        for (i = 0; i < 3; i++) {
                double weight, momentum[3] = { 0., 0., 1. };
                if (alouette_undecay(16, momentum, polarise, &weight) != 0)
                        continue;

                int pid;
                printf("# Event %d (%.5E) :\n", i + 1, weight);
                while (alouette_product(&pid, momentum)) {
                        printf("    %4d %12.5lE %12.5lE %12.5lE\n", pid,
                            momentum[0], momentum[1], momentum[2]);
                }
        }

        /* Clean the wrapper. */
        alouette_finalise();

        return 0;
}
