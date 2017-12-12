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
#include <setjmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* The ALOUETTE API. */
#include "alouette.h"

#ifndef M_PI
/* Define pi, if unknown. */
#define M_PI 3.14159265358979323846
#endif

/* Stack for TAUOLA's decay products. */
static struct {
#define STACK_MAX_DEPTH 16
        int length, index;
        int pid[STACK_MAX_DEPTH];
        double polarimetric[4];
        double p[4 * STACK_MAX_DEPTH];
} stack = { 0, -1 };

/* Jump buffer for calling back from a TAUOLA error. */
static jmp_buf jump_buffer;

/* Jump back to the calling context from within TAUOLA, instead of issuing
 * a hard stop.
 */
void softstp_() { longjmp(jump_buffer, 1); }

/* TAUOLA's common block for particle masses. */
extern struct {
        float amtau, amnuta, amell, amnue, ammu, amnumu, ampiz, ampi, amro,
            gamro, ama1, gama1, amk, amkz, amkst, gamkst;
} parmas_;

/* TAUOLA's main decay routine. */
extern void dekay_(int * state, double polarimetric[4]);

/* Initialisation of TAUOLA's engine. */
static enum alouette_return tauola_initialise(char * tauola_log, int seed)
{
        /* Register a rally in case of a TAUOLA error. */
        if (setjmp(jump_buffer) != 0) return ALOUETTE_RETURN_TAULOA_ERROR;

        if (tauola_log != NULL) {
                /* Redirect the output messages from TAUOLA. */
                extern void openout_(
                    int *, char *, int); /* TODO: check success. */
                int iout = 17;
                openout_(&iout, tauola_log, strlen(tauola_log));
        }

        /* Set the seed for the random engine. */
        extern void rmarin_(int * ijklin, int * ntotin, int * ntot2n);
        int ntotin = 0, ntot2n = 0;
        seed = seed % 900000001;
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

        return ALOUETTE_RETURN_SUCCESS;
}

/* Finalise TAUOLA's engine. */
static enum alouette_return tauola_finalise()
{
        /* Register a rally in case of a TAUOLA error. */
        if (setjmp(jump_buffer) != 0) return ALOUETTE_RETURN_TAULOA_ERROR;

        /* Dump a summary info. */
        int state = 100;
        double polarimetric[4];
        dekay_(&state, polarimetric);

        /* Close any file redirection. */
        extern void closout_();
        closout_();

        return ALOUETTE_RETURN_SUCCESS;
}

/* Low level decay routine with TAUOLA. */
static enum alouette_return tauola_decay(int pid, int pull)
{
        /* Register a rally in case of a TAUOLA error. */
        if (setjmp(jump_buffer) != 0) return ALOUETTE_RETURN_TAULOA_ERROR;

        /* Configure the positions of taus in the LUND common block */
        extern struct {
                int npa, npb;
        } taupos_;
        taupos_.npa = 1;
        taupos_.npb = 1;

        /* Call TAUOLA's decay routine. */
        int type = (pid > 0) ? 1 : 2;
        type += 10 * (pull != 0);
        dekay_(&type, stack.polarimetric);

        return ALOUETTE_RETURN_SUCCESS;
}

/* Callback for TAUOLA, used for retrieving decay products. */
void filhep_(int * n, int * status, int * pid, int * mother_first,
    int * mother_last, int * daughter_first, int * daughter_last, float p4[4],
    float * p_inv_mass, int * photos_flag)
{
        if (*status != 1) return;
        const int aid = abs(*pid);
        if ((aid == 24) || (aid == 313) || (aid == 323) || (aid > 9999)) return;
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
enum alouette_return alouette_initialise(int mute, int * seed_p)
{
        if (initialised) return ALOUETTE_RETURN_SUCCESS;

        /* Set the seed for the random engine. */
        int seed;
        if (seed_p == NULL) {
                seed = 0;
                FILE * stream = fopen("/dev/urandom", "rb");
                if (stream != NULL) {
                        if (fread(&seed, sizeof(seed), 1, stream) <= 0) {
                                fclose(stream);
                                return ALOUETTE_RETURN_IO_ERROR;
                        }
                        fclose(stream);
                } else
                        return ALOUETTE_RETURN_PATH_ERROR;
        } else {
                seed = *seed_p;
        }

        /* Initialise TAUOLA's library. */
        char * tauola_log = (mute != 0) ? "/dev/null" : NULL;
        enum alouette_return rc;
        if ((rc = tauola_initialise(tauola_log, seed)) !=
            ALOUETTE_RETURN_SUCCESS)
                return rc;
        initialised = 1;
        return ALOUETTE_RETURN_SUCCESS;
}

/* Finalise TAUOLA and its wrapper. */
enum alouette_return alouette_finalise(void)
{
        if (!initialised) return ALOUETTE_RETURN_SUCCESS;
        enum alouette_return rc;
        if ((rc = tauola_finalise()) != ALOUETTE_RETURN_SUCCESS) return rc;
        initialised = 0;
        return ALOUETTE_RETURN_SUCCESS;
}

/* Get a return code as a string. */
const char * alouette_strerror(enum alouette_return rc)
{
        static const char * msg[ALOUETTE_N_RETURNS] = { "Operation succeeded",
                "A value is out of range", "A floating point error occured",
                "Couldn't read or write file", "No such file or directory",
                "A Tauola error occured" };

        if ((rc < 0) || (rc >= ALOUETTE_N_RETURNS))
                return NULL;
        else
                return msg[rc];
}

/* Uniform PRN over [0,1]. */
static double uniform01(void)
{
        /* TODO: uniformise the PRN streams. */
        return rand() / (double)(RAND_MAX);
}

/* Decay a tau with TAUOLA. */
enum alouette_return alouette_decay(
    int pid, const double momentum[3], const double * polarisation)
{
        enum alouette_return rc;

        /* Reset the stack. */
        stack.length = 0;
        stack.index = -1;
        if (abs(pid) != 15) return ALOUETTE_RETURN_DOMAIN_ERROR;

        /* Decay a tau in its rest frame. */
        for (;;) {
                if ((rc = tauola_decay(pid, 0)) != ALOUETTE_RETURN_SUCCESS)
                        return rc;
                if (polarisation == NULL) break;
                const double w =
                    0.5 * (1. + stack.polarimetric[0] * polarisation[0] +
                              stack.polarimetric[1] * polarisation[1] +
                              stack.polarimetric[2] * polarisation[2]);
                if (uniform01() <= w) break;
                stack.length = 0;
        }
        if ((rc = tauola_decay(pid, 1)) != ALOUETTE_RETURN_SUCCESS) return rc;

        /* Check the result. */
        if (isinf(stack.p[0])) {
                stack.length = 0;
                return ALOUETTE_RETURN_FLOATING_ERROR;
        }

        /* Boost the daughters to the lab frame. */
        const double tau[3] = { momentum[0] / parmas_.amtau,
                momentum[1] / parmas_.amtau, momentum[2] / parmas_.amtau };
        const double gamma =
            sqrt(1. + tau[0] * tau[0] + tau[1] * tau[1] + tau[2] * tau[2]);
        stack.index = 0;
        rc = ALOUETTE_RETURN_SUCCESS;
        if (gamma <= 1. + FLT_EPSILON) return rc;

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

        return rc;
}

/* Multiply a matrix and a vector. */
static void mv_multiply(const double * M, double * v)
{
        const double v0 = M[0] * v[0] + M[1] * v[1] + M[2] * v[2];
        const double v1 = M[3] * v[0] + M[4] * v[1] + M[5] * v[2];
        const double v2 = M[6] * v[0] + M[7] * v[1] + M[8] * v[2];
        v[0] = v0;
        v[1] = v1;
        v[2] = v2;
}

/* Rotate a direction randomly but constraining the polar angle. */
static enum alouette_return rotate_direction(
    double cos_theta, double * direction)
{
        /* Check the numerical sine. */
        const double stsq = 1. - cos_theta * cos_theta;
        if (stsq <= 0.) return ALOUETTE_RETURN_FLOATING_ERROR;
        const double st = sqrt(stsq);

        /* select the co-vectors for the local basis. */
        double u0x = 0., u0y = 0., u0z = 0.;
        const double a0 = fabs(direction[0]);
        const double a1 = fabs(direction[1]);
        const double a2 = fabs(direction[2]);
        if (a0 > a1) {
                if (a0 > a2) {
                        const double nrm =
                            1. / sqrt(direction[0] * direction[0] +
                                     direction[2] * direction[2]);
                        u0x = -direction[2] * nrm, u0z = direction[0] * nrm;
                } else {
                        const double nrm =
                            1. / sqrt(direction[1] * direction[1] +
                                     direction[2] * direction[2]);
                        u0y = direction[2] * nrm, u0z = -direction[1] * nrm;
                }
        } else {
                if (a1 > a2) {
                        const double nrm =
                            1. / sqrt(direction[0] * direction[0] +
                                     direction[1] * direction[1]);
                        u0x = direction[1] * nrm, u0y = -direction[0] * nrm;
                } else {
                        const double nrm =
                            1. / sqrt(direction[1] * direction[1] +
                                     direction[2] * direction[2]);
                        u0y = direction[2] * nrm, u0z = -direction[1] * nrm;
                }
        }
        const double u1x = u0y * direction[2] - u0z * direction[1];
        const double u1y = u0z * direction[0] - u0x * direction[2];
        const double u1z = u0x * direction[1] - u0y * direction[0];

        /* Apply the rotation. */
        const double phi = M_PI * (1. - 2. * uniform01());
        const double cp = cos(phi);
        const double sp = sin(phi);
        direction[0] = cos_theta * direction[0] + st * (cp * u0x + sp * u1x);
        direction[1] = cos_theta * direction[1] + st * (cp * u0y + sp * u1y);
        direction[2] = cos_theta * direction[2] + st * (cp * u0z + sp * u1z);

        return ALOUETTE_RETURN_SUCCESS;
}

/* Build a rotation matrix from vi to vf. */
static enum alouette_return build_rotation(
        double * vi, double * vf, double norm_f, double * R)
{
        /* Build the rotation axis. */
        double n[3] = { vi[1] * vf[2] - vi[2] * vf[1],
            vi[2] * vf[0] - vi[0] * vf[2], vi[0] * vf[1] - vi[1] * vf[0] };
        double nrm = n[0] * n[0] + n[1] * n[1] + n[2] * n[2];
        if (fabs(nrm) <= FLT_EPSILON) return ALOUETTE_RETURN_FLOATING_ERROR;
        nrm = 1. / sqrt(nrm);
        n[0] *= nrm;
        n[1] *= nrm;
        n[2] *= nrm;

        /* Compute the rotation angle. */
        const double theta =
            -acos((vi[0] * vf[0] + vi[1] * vf[1] + vi[2] * vf[2]) / norm_f);

        /* Fill the rotation matrix. */
        const double c = cos(theta);
        const double s = sin(theta);
        const double c1 = 1. - c;
        R[0] = n[0] * n[0] * c1 + c;
        R[1] = n[0] * n[1] * c1 - n[2] * s;
        R[2] = n[0] * n[2] * c1 + n[1] * s;
        R[3] = n[0] * n[1] * c1 + n[2] * s;
        R[4] = n[1] * n[1] * c1 + c;
        R[5] = n[1] * n[2] * c1 - n[0] * s;
        R[6] = n[0] * n[2] * c1 - n[1] * s;
        R[7] = n[1] * n[2] * c1 + n[0] * s;
        R[8] = n[2] * n[2] * c1 + c;

        return ALOUETTE_RETURN_SUCCESS;
}

/* Backward decay from a tau neutrino to a tau. */
enum alouette_return alouette_undecay(int pid, const double momentum[3],
    alouette_polarisation_cb * polarisation, double bias, double * weight)
{
        enum alouette_return rc;

        /* Reset the stack. */
        stack.length = 0;
        stack.index = -1;
        if ((abs(pid) != 16) || (bias <= -1.))
            return ALOUETTE_RETURN_DOMAIN_ERROR;

        /* Decay an unpolarised tau in its rest frame. */
        if ((rc = tauola_decay(pid, 0)) != ALOUETTE_RETURN_SUCCESS) return rc;
        if ((rc = tauola_decay(pid, 1)) != ALOUETTE_RETURN_SUCCESS) return rc;

        /* Check the result. */
        if (isinf(stack.p[0])) {
                stack.length = 0;
                return ALOUETTE_RETURN_FLOATING_ERROR;
        }

        /* Get the daughter's index. */
        int i = 0;
        double * pi;
        for (i = 0, pi = stack.p; i < stack.length; i++, pi += 4)
                if (stack.pid[i] == pid) break;

        /* Re-draw the daughter's direction in the mother's rest frame using
         * a biased distribution, towards the daughter's direction in the
         * laboratory frame.
         */
        *weight = 1.;
        const double energy = sqrt(momentum[0] * momentum[0] +
            momentum[1] * momentum[1] + momentum[2] * momentum[2]);
        double * R = NULL;
        double R_storage[9];
        if (bias != 0.) {
                /* Draw the new direction. */
                double u[3] = { momentum[0] / energy, momentum[1] / energy,
                    momentum[2] / energy };
                const double r = uniform01();
                const double cos_theta = 2. * pow(r, 1. / (bias + 1.)) - 1.;
                if ((rotate_direction(cos_theta, u) ==
                        ALOUETTE_RETURN_SUCCESS) &&
                    (build_rotation(u, pi, pi[3], R_storage) ==
                        ALOUETTE_RETURN_SUCCESS)) {
                        /* Update the generated state and the BMC weight. */
                        R = R_storage;
                        mv_multiply(R, pi);
                        mv_multiply(R, stack.polarimetric);
                        *weight = 0.5 * (cos_theta + 1.) / ((bias + 1.) * r);
                }
        }

        /* Compute the parameters of the frame transform. */
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

                /* Apply the bias rotation. */
                if (R != NULL) mv_multiply(R, pj);

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
        *weight *= parmas_.amtau * parmas_.amtau * parmas_.amtau *
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
        return ALOUETTE_RETURN_SUCCESS;
}

/* Iterator over the tau decay products. */
enum alouette_return alouette_product(int * pid, double momentum[3])
{
        if ((stack.length <= 0) || (stack.index < 0)) return 0;
        if (stack.index < stack.length) {
                double * p = stack.p + 4 * stack.index;
                *pid = stack.pid[stack.index++];
                momentum[0] = p[0];
                momentum[1] = p[1];
                momentum[2] = p[2];
                return ALOUETTE_RETURN_SUCCESS;
        }
        stack.length = 0;
        return ALOUETTE_RETURN_DOMAIN_ERROR;
}

/* Getter for the polarimetric vector of the last decay. */
enum alouette_return alouette_polarimetric(double polarimetric[3])
{
        if (stack.length <= 0) return ALOUETTE_RETURN_DOMAIN_ERROR;
        polarimetric[0] = stack.polarimetric[0];
        polarimetric[1] = stack.polarimetric[1];
        polarimetric[2] = stack.polarimetric[2];
        return ALOUETTE_RETURN_SUCCESS;
}
