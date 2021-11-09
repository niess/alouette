/* Standard library includes. */
#include <stdio.h>
#include <stdlib.h>

/* ALOUETTE: a TAUOLA BMC wrapper. */
#include "alouette.h"

int main()
{
        /* Randomise a few tau decays. */
        const double momentum[3] = { 0., 0., 1. };
        const double polarisation[3] = { 1., 0., 0. };

        int i;
        for (i = 0; i < 3; i++) {
                struct alouette_products products;
                if (alouette_decay(15, momentum, polarisation, &products) !=
                    ALOUETTE_RETURN_SUCCESS)
                        continue;

                printf("# Event %d :\n", i + 1);
                int j;
                for (j = 0; j < products.size; j++) {
                        printf("    %4d %12.5lE %12.5lE %12.5lE\n",
                            products.pid[j], products.P[j][0], products.P[j][1],
                            products.P[j][2]);
                }
        }

        return 0;
}
