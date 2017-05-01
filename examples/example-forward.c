/* Standard library includes. */
#include <stdio.h>
#include <stdlib.h>

/* ALOUETTE: a TAUOLA BMC wrapper. */
#include "alouette.h"

int main()
{
        /* Initialise the ALOUETTE library. */
        enum alouette_return rc;
        if ((rc = alouette_initialise(1, NULL)) != ALOUETTE_RETURN_SUCCESS) {
                fprintf(stderr, "alouette: %s\n", alouette_strerror(rc));
                exit(EXIT_FAILURE);
        };

        /* Randomise a few tau decays. */
        const double polarisation[3] = { 1., 0., 0. };
        int i;
        for (i = 0; i < 3; i++) {
                double momentum[3] = { 0., 0., 1. };
                if (alouette_decay(15, momentum, polarisation) !=
                    ALOUETTE_RETURN_SUCCESS)
                        continue;

                int pid;
                printf("# Event %d :\n", i + 1);
                while (alouette_product(&pid, momentum) ==
                    ALOUETTE_RETURN_SUCCESS) {
                        printf("    %4d %12.5lE %12.5lE %12.5lE\n", pid,
                            momentum[0], momentum[1], momentum[2]);
                }
        }

        /* Finalise the library. */
        if ((rc = alouette_finalise()) != ALOUETTE_RETURN_SUCCESS) {
                fprintf(stderr, "alouette: %s\n", alouette_strerror(rc));
                exit(EXIT_FAILURE);
        };

        return 0;
}
