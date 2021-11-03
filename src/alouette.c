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
#include <stdarg.h>
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
#define STACK_MAX_SIZE 16
        int size, index;
        int pid[STACK_MAX_SIZE];
        double polarimetric[4];
        double p[4 * STACK_MAX_SIZE];
} _stack = { 0, -1 };

/* Stack for message(s) */
static struct {
#define MESSAGE_MAX_SIZE 1024
        int size;
        enum alouette_return code;
        char data[];
} _message = {0, ALOUETTE_RETURN_SUCCESS, {0x0}};

/* Some TAUOLA common blocks and routine(s) */
extern struct {
    int jak1, jak2, jakp, jakm, ktom;
} tauola_jaki;

extern struct {
    float amtau, amnuta, amel, amnue, ammu, amnumu, ampiz, ampi, amro, gamro,
          ama1, gama1, amk, amkz, amkst, gamkst;
} tauola_parmas;

extern void tauola_decay(int * state, double polarimeter[4]);

/* Jump buffer for calling back from a TAUOLA error. */
static jmp_buf _jump_buffer;

/* Jump back to the calling context from within TAUOLA, instead of issuing
 * a hard stop.
 */
void tauola_stop(void)
{
        longjmp(_jump_buffer, 1);
}

/* Redirect TAUOLA printing to the message stack. */
void tauola_print(const char * msg)
{
        const char * prefix = (_message.size > 0) ? "\n" : "";
        _message.size += snprintf(_message.data + _message.size,
            MESSAGE_MAX_SIZE - _message.size - 1, "%s%s", prefix, msg);
}

/* Callback for TAUOLA, used for retrieving decay products. */
void tauola_filhep(int * n, int * status, int * pid, int * mother_first,
    int * mother_last, int * daughter_first, int * daughter_last, float p4[4],
    float * p_inv_mass, int * photos_flag)
{
        if ((*status != 1) || (_stack.size >= STACK_MAX_SIZE))
                return;

        const int aid = abs(*pid);
        if ((aid == 24) || (aid == 313) || (aid == 323) || (aid > 9999)) return;
        _stack.pid[_stack.size] = *pid;
        double * const P = _stack.p + 4 * _stack.size;
        P[0] = (double)p4[0];
        P[1] = (double)p4[1];
        P[2] = (double)p4[2];
        P[3] = (double)p4[3];
        _stack.size++;
}

/* Utility function for reseting the message stack. */
static void message_reset(void)
{
        _message.size = 0;
        _message.code = ALOUETTE_RETURN_SUCCESS;
        _message.data[0] = 0x0;
}

/* Utility function for dumping an error to the message stack. */
static enum alouette_return message_error(
    enum alouette_return code, const char * fmt, ...)
{
        _message.code = code;

        if (fmt != NULL) {
                va_list args;
                va_start(args, fmt);
                _message.size += vsnprintf(_message.data,
                    MESSAGE_MAX_SIZE - 1 - _message.size, fmt, args);
                va_end(args);
        }

        return code;
}

/* Data structure for the built-in MT pseudo random engine. */
#define MT_PERIOD 624
static struct {
        unsigned long seed;
        int index;
        unsigned long data[MT_PERIOD];
} _random_stream;

/* Get the random seed for the built-in PRNG. */
unsigned long alouette_seed_get(void)
{
        return _random_stream.seed;
}

/* Get a random seed from /dev/urandom */
static enum alouette_return random_get_seed(unsigned long * seed)
{
        const char * urandom = "/dev/urandom";
        FILE * fp = fopen(urandom, "rb");
        if (fp == NULL) {
                return message_error(ALOUETTE_RETURN_PATH_ERROR,
                    "could not open %s", urandom);
        }
        if (fread(seed, sizeof(long), 1, fp) <= 0) {
                fclose(fp);
                return message_error(ALOUETTE_RETURN_IO_ERROR,
                    "could not read from %s", urandom);
        }
        fclose(fp);
        return ALOUETTE_RETURN_SUCCESS;
}

/* Initialise the built-in PRNG. */
static enum alouette_return random_initialise(unsigned long * seed_ptr)
{
        if (seed_ptr == NULL) {
                /* Get a seed from /dev/urandom */
                enum alouette_return rc;
                if ((rc = random_get_seed(&_random_stream.seed)) !=
                    ALOUETTE_RETURN_SUCCESS)
                        return rc;
        } else {
                _random_stream.seed = *seed_ptr;
        }

        /* Set the Mersenne Twister initial state. */
        _random_stream.data[0] = _random_stream.seed & 0xffffffffUL;
        int j;
        for (j = 1; j < MT_PERIOD; j++) {
                _random_stream.data[j] = (1812433253UL *
                        (_random_stream.data[j - 1] ^
                            (_random_stream.data[j - 1] >> 30)) +
                    j);
                _random_stream.data[j] &= 0xffffffffUL;
        }
        _random_stream.index = MT_PERIOD;

        return ALOUETTE_RETURN_SUCCESS;
}

/* Set the random seed for the built-in PRNG. */
enum alouette_return alouette_seed_set(unsigned long * seed)
{
        message_reset();
        return random_initialise(seed);
}

/* Uniform pseudo random distribution over (0,1) from a Mersenne Twister */
static float random_uniform01(void)
{
        /* Check the buffer */
        if (_random_stream.index < MT_PERIOD - 1) {
                _random_stream.index++;
        } else {
                /* Update the MT state */
                const int M = 397;
                const unsigned long UPPER_MASK = 0x80000000UL;
                const unsigned long LOWER_MASK = 0x7fffffffUL;
                static unsigned long mag01[2] = { 0x0UL, 0x9908b0dfUL };
                unsigned long y;
                int kk;
                for (kk = 0; kk < MT_PERIOD - M; kk++) {
                        y = (_random_stream.data[kk] & UPPER_MASK) |
                            (_random_stream.data[kk + 1] & LOWER_MASK);
                        _random_stream.data[kk] =
                            _random_stream.data[kk + M] ^ (y >> 1) ^
                            mag01[y & 0x1UL];
                }
                for (; kk < MT_PERIOD - 1; kk++) {
                        y = (_random_stream.data[kk] & UPPER_MASK) |
                            (_random_stream.data[kk + 1] & LOWER_MASK);
                        _random_stream.data[kk] =
                            _random_stream.data[kk + (M - MT_PERIOD)] ^
                            (y >> 1) ^ mag01[y & 0x1UL];
                }
                y = (_random_stream.data[MT_PERIOD - 1] & UPPER_MASK) |
                    (_random_stream.data[0] & LOWER_MASK);
                _random_stream.data[MT_PERIOD - 1] =
                    _random_stream.data[M - 1] ^ (y >> 1) ^
                    mag01[y & 0x1UL];
                _random_stream.index = 0;
        }

        /* Tempering */
        unsigned long y = _random_stream.data[_random_stream.index];
        y ^= (y >> 11);
        y ^= (y << 7) & 0x9d2c5680UL;
        y ^= (y << 15) & 0xefc60000UL;
        y ^= (y >> 18);

        /* Convert to a 32 bits floating point in (0,1) */
        const double u = (((double)y) + 0.5) * (1.0 / 4294967296.0);

        return (float)u;
}

/* Random stream wrapper for TAUOLA */
void tauola_random(float * r, int * n)
{
        int i;
        for (i = 0; i < *n; i++) {
                r[i] = random_uniform01();
        }
}

/* Status flag for the wrapper's initialisation. */
static int _initialised = 0;

/* Initialise TAUOLA and its wrapper. */
enum alouette_return alouette_initialise(
    unsigned long * seed, double * xk0dec)
{
        message_reset();
        if (_initialised) return ALOUETTE_RETURN_SUCCESS;

        /* Rally point in case of a TAUOLA error. */
        if (setjmp(_jump_buffer) != 0) {
                return message_error(ALOUETTE_RETURN_TAULOA_ERROR, NULL);
        }

        /* TAUOLA common blocks and routines. */
        extern struct { int idff; } tauola_idfc;
        extern struct { double xk0dec; int itdkrc; } tauola_taurad;

        /* Parse settings and configure. */
        if (xk0dec == NULL) {
                /* Cut parameter for radiative corrections in leptonic
                 * decays. Set default value according to Jezabek et al.,
                 * CPC 70 (1992) 69-76.
                 */
                tauola_taurad.itdkrc = 1;
                tauola_taurad.xk0dec = 0.001;
        } else if (*xk0dec > 0.) {
                tauola_taurad.itdkrc = 1;
                tauola_taurad.xk0dec = *xk0dec;
        } else {
                tauola_taurad.itdkrc = 0;
                tauola_taurad.xk0dec = 0.;
        }

        /* PDG identifier for tau. */
        tauola_idfc.idff = 15;

        /* Configure the positions of taus in TAUOLA's stack. */
        extern struct {
                int npa, npb;
        } tauola_taupos;
        tauola_taupos.npa = 1;
        tauola_taupos.npb = 1;

        /* Initialise the builtin random engine with an arbitrary seed for
         * TAUOLA warmup.
         */
        unsigned long tmp = 1357894;
        enum alouette_return rc;
        if ((rc = random_initialise(&tmp)) != ALOUETTE_RETURN_SUCCESS)
                return rc;

        /* Initialise the decay routine. */
        double polarimeter[4];
        int state = -1;
        tauola_jaki.jak1 = 0;
        tauola_jaki.jak2 = 0;
        tauola_decay(&state, polarimeter);

        /* Clear the particle stack. */
        memset(&_stack, 0x0, sizeof _stack);

        /* Initialise the random engine with the the user supplied seed. */
        if ((rc = random_initialise(seed)) != ALOUETTE_RETURN_SUCCESS)
                return rc;

        /* Flag as initialised and return. */
        _initialised = 1;

        return ALOUETTE_RETURN_SUCCESS;
}

/* Get the last (error) message(s). */
const char * alouette_message(void)
{
        static const char * msg[ALOUETTE_N_RETURNS] = { "Operation succeeded",
                "A value is out of range", "A floating point error occured",
                "Couldn't read or write file", "No such file or directory",
                "A Tauola error occured" };

        if (_message.code == ALOUETTE_RETURN_SUCCESS) {
                if (_message.data[0] == 0x0) {
                        return NULL;
                } else {
                        return _message.data;
                }
        } else if (_message.code == ALOUETTE_RETURN_TAULOA_ERROR) {
                        return _message.data;
        } else {
                return msg[_message.code];
        }
}

/* Low level decay routine with TAUOLA. */
static enum alouette_return decay(int pid, int pull)
{
        /* Register a rally in case of a TAUOLA error. */
        if (setjmp(_jump_buffer) != 0) {
                return message_error(ALOUETTE_RETURN_TAULOA_ERROR, NULL);
        }

        /* XXX Allow to select the decay mode. */

        /* Call TAUOLA's decay routine. */
        int type = (pid > 0) ? 1 : 2;
        type += 10 * (pull != 0);
        tauola_decay(&type, _stack.polarimetric);

        return ALOUETTE_RETURN_SUCCESS;
}

/* Decay a tau with TAUOLA. */
enum alouette_return alouette_decay(
    int pid, const double momentum[3], const double * polarisation)
{
        enum alouette_return rc;

        /* Reset the stacks. */
        message_reset();
        _stack.size = 0;
        _stack.index = -1;
        if (abs(pid) != 15) {
                return message_error(ALOUETTE_RETURN_DOMAIN_ERROR,
                    "bad pid for mother particle (%d)", pid);
        }

        /* Decay a tau in its rest frame. */
        for (;;) {
                /* XXX Rotate instead of re-decaying */
                if ((rc = decay(pid, 0)) != ALOUETTE_RETURN_SUCCESS)
                        return rc;
                if (polarisation == NULL) break;
                const double w =
                    0.5 * (1. + _stack.polarimetric[0] * polarisation[0] +
                              _stack.polarimetric[1] * polarisation[1] +
                              _stack.polarimetric[2] * polarisation[2]);
                if (random_uniform01() <= w) break;
        }
        if ((rc = decay(pid, 1)) != ALOUETTE_RETURN_SUCCESS) return rc;

        /* Check the result. */
        if (isinf(_stack.p[0])) {
                _stack.size = 0;
                return message_error(ALOUETTE_RETURN_FLOATING_ERROR,
                        "floating point exception (inf)");
        }

        /* Boost the daughters to the lab frame. */
        const double amtau = tauola_parmas.amtau;
        const double tau[3] = { momentum[0] / amtau,
                momentum[1] / amtau, momentum[2] / amtau };
        const double gamma =
            sqrt(1. + tau[0] * tau[0] + tau[1] * tau[1] + tau[2] * tau[2]);
        _stack.index = 0;
        rc = ALOUETTE_RETURN_SUCCESS;
        if (gamma <= 1. + FLT_EPSILON) return rc;

        int i;
        double * P;
        for (i = 0, P = _stack.p; i < _stack.size; i++, P += 4) {
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
        if (stsq <= 0.) {
                return message_error(ALOUETTE_RETURN_FLOATING_ERROR,
                        "floating point exception (%g <= 0.)", stsq);
        }
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
        const double phi = M_PI * (1. - 2. * random_uniform01());
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
        if (fabs(nrm) <= FLT_EPSILON) {
                return message_error(ALOUETTE_RETURN_FLOATING_ERROR,
                    "floating point exception (%g <= %g)", fabs(nrm),
                    FLT_EPSILON);
        }
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

        /* Reset the stacks. */
        message_reset();
        _stack.size = 0;
        _stack.index = -1;
        if (abs(pid) != 16) {
                return message_error(ALOUETTE_RETURN_DOMAIN_ERROR,
                        "bad pid value for daugther particle (%d)", pid);
        } else if (bias <= -1.) {
                return message_error(ALOUETTE_RETURN_DOMAIN_ERROR,
                        "bad bias value (%g <= -1)", bias);
        }

        /* Decay an unpolarised tau in its rest frame. */
        if ((rc = decay(pid, 0)) != ALOUETTE_RETURN_SUCCESS) return rc;
        if ((rc = decay(pid, 1)) != ALOUETTE_RETURN_SUCCESS) return rc;

        /* Check the result. */
        if (isinf(_stack.p[0])) {
                _stack.size = 0;
                return message_error(ALOUETTE_RETURN_FLOATING_ERROR,
                        "floating point exception (inf)");
        }

        /* Get the daughter's index. */
        int i = 0;
        double * pi;
        for (i = 0, pi = _stack.p; i < _stack.size; i++, pi += 4)
                if (_stack.pid[i] == pid) break;

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
                const double r = random_uniform01();
                const double cos_theta = 2. * pow(r, 1. / (bias + 1.)) - 1.;
                if ((rotate_direction(cos_theta, u) ==
                        ALOUETTE_RETURN_SUCCESS) &&
                    (build_rotation(u, pi, pi[3], R_storage) ==
                        ALOUETTE_RETURN_SUCCESS)) {
                        /* Update the generated state and the BMC weight. */
                        R = R_storage;
                        mv_multiply(R, pi);
                        mv_multiply(R, _stack.polarimetric);
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
        for (j = 0, pj = _stack.p; j < _stack.size; j++, pj += 4) {
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
        const double amtau = tauola_parmas.amtau;
        *weight *= amtau * amtau * amtau * fabs(t0 * t0 * gamma / energy);

        /* Weight for the spin polarisation of the tau mother. */
        const int tau_pid = pid > 0 ? 15 : -15;
        if (polarisation != NULL) {
                double pol[3];
                polarisation(tau_pid, Pt, pol);
                *weight *= 1. + _stack.polarimetric[0] * pol[0] +
                    _stack.polarimetric[1] * pol[1] +
                    _stack.polarimetric[2] * pol[2];
        }

        /* Replace the daughter particle with the tau mother. In order to
         * ensure the conservation of the energy-momentum, the mother's 4
         * momentum has been computed from the boosted products.
         */
        _stack.pid[i] = tau_pid;
        pi[0] = Pt[0];
        pi[1] = Pt[1];
        pi[2] = Pt[2];
        pi[3] = Et;

        /* Set the stack index and return. */
        _stack.index = 0;
        return ALOUETTE_RETURN_SUCCESS;
}

/* Iterator over the tau decay products. */
enum alouette_return alouette_product(int * pid, double momentum[3])
{
        if ((_stack.size <= 0) || (_stack.index < 0)) return 0;
        if (_stack.index < _stack.size) {
                double * p = _stack.p + 4 * _stack.index;
                *pid = _stack.pid[_stack.index++];
                momentum[0] = p[0];
                momentum[1] = p[1];
                momentum[2] = p[2];
                return ALOUETTE_RETURN_SUCCESS;
        }
        _stack.size = 0;
        return ALOUETTE_RETURN_DOMAIN_ERROR;
}

/* Getter for the polarimetric vector of the last decay. */
enum alouette_return alouette_polarimetric(double polarimetric[3])
{
        if (_stack.size <= 0) return ALOUETTE_RETURN_DOMAIN_ERROR;
        polarimetric[0] = _stack.polarimetric[0];
        polarimetric[1] = _stack.polarimetric[1];
        polarimetric[2] = _stack.polarimetric[2];
        return ALOUETTE_RETURN_SUCCESS;
}
