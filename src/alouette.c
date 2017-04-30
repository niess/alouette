/*
 * Copyright (C) 2017 Université Clermont Auvergne, CNRS/IN2P3, LPC
 * Author: Valentin NIESS (niess@in2p3.fr)
 *
 * This software is a C library for the backward decay of taus wrapping TAUOLA.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>
 */

/* Standard library includes. */
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* The ALOUETTE API. */
#include "alouette.h"

/* Stack for TAUOLA's decay products. */
static struct {
#define STACK_MAX_DEPTH 16
        int length, index;
        int pid[STACK_MAX_DEPTH];
        double polarimetric[4];
        double p[4 * STACK_MAX_DEPTH];
} stack = { 0, -1 };

/* TAUOLA's common block for particle masses. */
extern struct {
        float amtau, amnuta, amell, amnue, ammu, amnumu, ampiz, ampi, amro,
            gamro, ama1, gama1, amk, amkz, amkst, gamkst;
} parmas_;

/* TAUOLA's main decay routine. */
extern void dekay_(int * state, double polarimetric[4]);

/* Initialisation of TAUOLA's engine. */
static void tauola_initialise(char * tauola_log, int seed)
{
        if (tauola_log != NULL) {
                /* Redirect the output messages from TAUOLA. */
                extern void openout_(int *, char *, int);
                int iout = 17;
                openout_(&iout, tauola_log, strlen(tauola_log));
        }

        /* Set the seed for the random engine. */
        extern void rmarin_(int * ijklin, int * ntotin, int * ntot2n);
        int ntotin = 0, ntot2n = 0;
        rmarin_(&seed, &ntotin, &ntot2n);

        /* Configure the common block for selecting the decay modes. */
        extern struct {
                int jak1, jak2, jakp, jakm, ktom;
        } jaki_;
        jaki_.jak1 = 0;
        jaki_.jak2 = 0;

        /* Unknown common block. */
        extern struct {
                double xk0dec;
                int itdkrc;
        } taurad_;
        taurad_.itdkrc = 1;
        taurad_.xk0dec = 0.001;

        /* Configure the common block with tau PDG ID. */
        extern struct {
                int idff;
        } idfc_;
        idfc_.idff = 15;

        /* Initialise masses. */
        extern void inimas_();
        inimas_();

        /* Initialise decay parameters? */
        extern void initdk_();
        initdk_();

        /* Initialise Physics? */
        extern void iniphy_(float * i);
        float iniphy_parameter = 0.1;
        iniphy_(&iniphy_parameter);

        /* Initialise the decay generator. */
        int state = -1;
        double polarimetric[4];
        dekay_(&state, polarimetric);
}

/* Finalise TAUOLA's engine. */
static void tauola_finalise()
{
        /* Dump a summary info. */
        int state = 100;
        double polarimetric[4];
        dekay_(&state, polarimetric);

        /* Close any file redirection. */
        extern void closout_();
        closout_();
}

/* Low level decay routine with TAUOLA. */
static void tauola_decay(int pid, int pull)
{
        /* Configure the positions of taus in the LUND common block */
        extern struct {
                int npa, npb;
        } taupos_;
        taupos_.npa = 1;
        taupos_.npb = 1;

        int type = (pid > 0) ? 1 : 2;
        type += 10 * (pull != 0);
        dekay_(&type, stack.polarimetric);
}

/* Callback for TAUOLA, used for retrieving decay products. */
void filhep_(int * n, int * status, int * pid, int * mother_first,
    int * mother_last, int * daughter_first, int * daughter_last, float p4[4],
    float * p_inv_mass, int * photos_flag)
{
        if (*status != 1) return;
        const int aid = abs(*pid);
        if ((aid == 24) || (aid > 9999)) return;
        stack.pid[stack.length] = *pid;
        double * const P = stack.p + 4 * stack.length;
        P[0] = (double)p4[0];
        P[1] = (double)p4[1];
        P[2] = (double)p4[2];
        P[3] = (double)p4[3];
        stack.length++;
}

/* Callback for TAUOLA, not needed since all decays are in the center of mass
 * frame.
 */
void tralo4_(float * kto, float p[4], float q[4], float * ams) {}

/* Status flag for the wrapper's initialisation. */
static int initialised = 0;

/* Initialise TAUOLA and its wrapper. */
void alouette_initialise(int mute, int * seed_p)
{
        if (initialised) return;

        /* Set the seed for the random engine. */
        int seed;
        if (seed_p == NULL) {
                seed = 0;
                FILE * stream = fopen("/dev/urandom", "rb");
                if (stream != NULL) {
                        if (fread(&seed, sizeof(seed), 1, stream) <= 0) {
                                /* TODO: manage the failure in reading from
                                 * /dev/urandom.
                                 */
                        };
                        fclose(stream);
                }
        } else {
                seed = *seed_p;
        }

        /* Initialise TAUOLA's library. */
        char * tauola_log = (mute != 0) ? "/dev/null" : NULL;
        tauola_initialise(tauola_log, seed);
        initialised = 1;
}

/* Finalise TAUOLA and its wrapper. */
void alouette_finalise(void)
{
        if (!initialised) return;
        tauola_finalise();
        initialised = 0;
}

/* Uniform PRN over [0,1]. */
static double uniform01(void)
{
        /* TODO: uniformise the PRN streams. */
        return rand() / (double)(RAND_MAX);
}

/* Decay a tau with TAUOLA. */
int alouette_decay(
    int pid, const double momentum[3], const double * polarisation)
{
        /* Reset the stack. */
        stack.length = 0;
        stack.index = -1;
        if (abs(pid) != 15) return stack.index;

        /* Decay a tau in its rest frame. */
        for (;;) {
                tauola_decay(pid, 0);
                if (polarisation == NULL) break;
                const double w =
                    0.5 * (1. + stack.polarimetric[0] * polarisation[0] +
                              stack.polarimetric[1] * polarisation[1] +
                              stack.polarimetric[2] * polarisation[2]);
                if (uniform01() <= w) break;
                stack.length = 0;
        }
        tauola_decay(pid, 1);

        /* Check the result. */
        if (isinf(stack.p[0])) {
                stack.length = 0;
                return stack.index;
        }

        /* Boost the daughters to the lab frame. */
        const double tau[3] = { momentum[0] / parmas_.amtau,
                momentum[1] / parmas_.amtau, momentum[2] / parmas_.amtau };
        const double gamma =
            sqrt(1. + tau[0] * tau[0] + tau[1] * tau[1] + tau[2] * tau[2]);
        stack.index = 0;
        if (gamma <= 1. + FLT_EPSILON) return stack.index;

        int i;
        double * P;
        for (i = 0, P = stack.p; i < stack.length; i++, P += 4) {
                /* Apply the boost. */
                const double ptau =
                    P[0] * tau[0] + P[1] * tau[1] + P[2] * tau[2];
                const double tmp = ptau / (gamma + 1.) + P[3];
                P[0] += tmp * tau[0];
                P[1] += tmp * tau[1];
                P[2] += tmp * tau[2];
                P[3] = gamma * P[3] + ptau;
        }

        return stack.index;
}

/* Backward decay from a tau neutrino to a tau. */
int alouette_undecay(int pid, const double momentum[3],
    polarisation_cb * polarisation, double * weight)
{
        /* Reset the stack. */
        stack.length = 0;
        stack.index = -1;
        if (abs(pid) != 16) return stack.index;

        /* Decay an unpolarised tau in its rest frame. */
        tauola_decay(pid, 0);
        tauola_decay(pid, 1);

        /* Check the result. */
        if (isinf(stack.p[0])) {
                stack.length = 0;
                return stack.index;
        }

        /* Get the daughter's index. */
        int i = 0;
        double * pi;
        for (i = 0, pi = stack.p; i < stack.length; i++, pi += 4)
                if (stack.pid[i] == pid) break;

        /* Compute the parameters of the frame transform. */
        const double energy = sqrt(momentum[0] * momentum[0] +
            momentum[1] * momentum[1] + momentum[2] * momentum[2]);
        const double d = energy * pi[3] + momentum[0] * pi[0] +
            momentum[1] * pi[1] + momentum[2] * pi[2];
        const double ee = energy + pi[3];
        const double gamma = ee * ee / d - 1.;
        const double t0 = ee / d;
        const double tau[3] = { t0 * (momentum[0] - pi[0]),
                t0 * (momentum[1] - pi[1]), t0 * (momentum[2] - pi[2]) };

        /* Boost the daughters to the lab frame. */
        double Et = 0., Pt[3] = { 0., 0., 0. };
        int j;
        double * pj;
        for (j = 0, pj = stack.p; j < stack.length; j++, pj += 4) {
                if (j == i) {
                        /* Update with the provided daughter data. */
                        Et += energy;
                        Pt[0] += momentum[0];
                        Pt[1] += momentum[1];
                        Pt[2] += momentum[2];
                        continue;
                }

                /* Apply the boost and update the sum. */
                const double ptau =
                    pj[0] * tau[0] + pj[1] * tau[1] + pj[2] * tau[2];
                const double tmp = ptau / (gamma + 1.) + pj[3];
                pj[0] += tmp * tau[0];
                pj[1] += tmp * tau[1];
                pj[2] += tmp * tau[2];
                pj[3] = gamma * pj[3] + ptau;
                Pt[0] += pj[0];
                Pt[1] += pj[1];
                Pt[2] += pj[2];
                Et += pj[3];
        }

        /* Set the backward Monte-Carlo weight. */
        *weight = parmas_.amtau * parmas_.amtau * parmas_.amtau *
            fabs(t0 * t0 * gamma / energy);

        /* Weight for the spin polarisation of the tau mother. */
        const int tau_pid = pid > 0 ? 15 : -15;
        if (polarisation != NULL) {
                double pol[3];
                polarisation(tau_pid, Pt, pol);
                *weight *= 1. + stack.polarimetric[0] * pol[0] +
                    stack.polarimetric[1] * pol[1] +
                    stack.polarimetric[2] * pol[2];
        }

        /* Replace the daughter particle with the tau mother. In order to
         * ensure the conservation of the energy-momentum, the mother's 4
         * momentum has been computed from the boosted products.
         */
        stack.pid[i] = tau_pid;
        pi[0] = Pt[0];
        pi[1] = Pt[1];
        pi[2] = Pt[2];
        pi[3] = Et;

        /* Set the stack index and return. */
        stack.index = 0;
        return stack.index;
}

/* Iterator over the tau decay products. */
int alouette_product(int * pid, double momentum[3])
{
        if ((stack.length <= 0) || (stack.index < 0)) return 0;
        if (stack.index < stack.length) {
                double * p = stack.p + 4 * stack.index;
                *pid = stack.pid[stack.index];
                momentum[0] = p[0];
                momentum[1] = p[1];
                momentum[2] = p[2];
                return ++stack.index;
        }
        stack.length = 0;
        return 0;
}

/* Getter for the polarimetric vector of the last decay. */
void alouette_polarimetric(double polarimetric[3])
{
        if (stack.length <= 0) return;
        polarimetric[0] = stack.polarimetric[0];
        polarimetric[1] = stack.polarimetric[1];
        polarimetric[2] = stack.polarimetric[2];
}
