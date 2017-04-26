/* Standard library includes. */
#include <stdio.h>

/* TAUOLA C wrapper. */
#include "alouette.h"

int main()
{
        /* Initialise TAUOLA through its wrapper. */
        alouette_initialise(1, NULL);

        /* Randomise a few tau decays. */
        const double polarisation[3] = { 1., 0., 0. };
        int i;
        for (i = 0; i < 3; i++) {
                double momentum[3] = { 0., 0., 1. };
                if (!alouette_decay(15, momentum, polarisation)) continue;

                int pid;
                printf("# Event %d :\n", i + 1);
                while (alouette_product(&pid, momentum)) {
                        printf("    %4d %12.5lE %12.5lE %12.5lE\n", pid,
                            momentum[0], momentum[1], momentum[2]);
                }
        }

        /* Clean the wrapper. */
        alouette_finalise();

        return 0;
}
