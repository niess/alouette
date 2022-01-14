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

/* The Alouette API. */
#include "alouette.h"

#ifndef M_PI
/* Define pi, if unknown. */
#define M_PI 3.14159265358979323846
#endif

/*  Extra debugging options, on linux. */
#ifdef ALOUETTE_DEBUG
#ifndef __USE_GNU
#define __USE_GNU
#endif
#include <fenv.h>
#endif

/* Some TAUOLA common blocks and routine(s) */
extern struct {
    int jak1, jak2, jakp, jakm, ktom;
} tauola_jaki;

struct tauola_taukle {
    float bra1, brk0, brk0b, brks;
};

extern struct tauola_taukle tauola_taukle;

extern struct {
    float amtau, amnuta, amel, amnue, ammu, amnumu, ampiz, ampi, amro, gamro,
          ama1, gama1, amk, amkz, amkst, gamkst;
} tauola_parmas;

extern void tauola_decay(int * state, double polarimeter[4]);

/* Jump buffer for calling back from a TAUOLA error.
 *
 * Note: setting this to static results in undefined behaviour.
 */
jmp_buf alouette_context;

/* Jump back to the calling context from within TAUOLA, instead of issuing
 * a hard stop.
 */
void tauola_stop(void)
{
        longjmp(alouette_context, 1);
}

/* Stack for message(s) */
static struct {
#define MESSAGE_MAX_SIZE 1024
        int size;
        enum alouette_return code;
        char data[MESSAGE_MAX_SIZE];
} _message = {0, ALOUETTE_RETURN_SUCCESS, {0x0}};

/* Redirect TAUOLA printing to the message stack. */
void tauola_print(const char * msg)
{
        if (_message.size >= MESSAGE_MAX_SIZE - 2) return;

        const char * prefix = (_message.size > 0) ? "\n" : "";
        _message.size += snprintf(_message.data + _message.size,
            MESSAGE_MAX_SIZE - _message.size - 1, "%s%s", prefix, msg);
}

/* Pointer to current decay products. */
static struct alouette_products * _products = NULL;

/* Callback for TAUOLA, used for retrieving decay products. */
void tauola_filhep(int * n, int * status, int * pid, int * mother_first,
    int * mother_last, int * daughter_first, int * daughter_last, float p4[4],
    float * p_inv_mass, int * photos_flag)
{
        const int i = _products->size;
        if ((*status != 1) || (i >= ALOUETTE_MAX_SIZE))
                return;

        const int aid = abs(*pid);
        if ((aid == 24) || (aid == 313) || (aid == 323) || (aid > 9999)) return;
        _products->pid[i] = *pid;
        int j;
        for (j = 0; j < 4; j++) {
                _products->P[i][j] = (double)p4[j];
        }
        _products->size++;
}

/* Utility function for reseting a products container. */
static void products_reset(struct alouette_products * products)
{
        products->size = 0;
        products->weight = 0.;
        memset(products->polarimeter, 0x0, sizeof products->polarimeter);
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
                    MESSAGE_MAX_SIZE - 1, fmt, args);
                va_end(args);
        }

        return code;
}

/* Data structure for the built-in MT pseudo random engine. */
#define MT_PERIOD 624
struct random_stream {
        int initialised;
        unsigned long seed;
        int index;
        unsigned long data[MT_PERIOD];
};

static struct random_stream _random_stream = {0};

/* Get the random seed of the built-in PRNG. */
unsigned long alouette_random_seed(void)
{
        if (!_random_stream.initialised) {
                alouette_random_set(NULL);
        }

        return _random_stream.seed;
}

/* Get a random seed from the OS, e.g. /dev/urandom */
static unsigned long random_get_seed(void)
{
        FILE * fp = fopen("/dev/urandom", "rb");
        if (fp != NULL) {
                unsigned long seed;
                size_t n = fread(&seed, sizeof(long), 1, fp);
                fclose(fp);
                if (n == 1) return seed;
        }

        return 0; /* XXX Fallback method? */
}

/* Set the random seed for the built-in PRNG. */
void alouette_random_set(unsigned long * seed)
{
        if (seed == NULL) {
                /* Get a seed from the OS. */
                _random_stream.seed = random_get_seed();
        } else {
                _random_stream.seed = *seed;
        }
        _random_stream.initialised = 1;

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
}

/* Uniform pseudo random distribution over (0,1) from a Mersenne Twister */
static float random_uniform01(void)
{
        /* Initialise the PRNG, if not already done. */
        if (!_random_stream.initialised) {
                alouette_random_set(NULL);
        }

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

/* Export the built-in PRNG. */
typedef float random_cb(void);
random_cb * alouette_random = &random_uniform01;

/* Random stream wrapper for TAUOLA */
void tauola_random(float * r, int * n)
{
        int i;
        for (i = 0; i < *n; i++) {
                r[i] = alouette_random();
        }
}

/* Data relative to decay channels. */
#define TAUOLA_MAX_CHANNELS 30
#define N_CHANNELS 30 /* Effective number of channels */
#define N_DAUGHTERS 13
static struct {
        int n;
        double total_width;
        double total_weight[N_DAUGHTERS];
        int mode[N_CHANNELS];
        int subchannel[N_CHANNELS];
        double branching_ratio[N_CHANNELS];
        double multiplicity[N_CHANNELS][N_DAUGHTERS];
        double weight[N_CHANNELS][N_DAUGHTERS];
        struct tauola_taukle default_taukle;
} _channels = {0};

/* Compute the array index for a given daughter pid. */
static int daughter_index(int pid)
{
        const int pids[N_DAUGHTERS] = {11, -12, 13, -14, 16, 111, 211, -211,
            221, 310, 130, 321, -321};

        int i;
        for (i = 0; i < N_DAUGHTERS; i++) {
                if (pid == pids[i]) return i;
        }
        return -1;
}

/* Converse (anti)particle. */
static int daughter_converse(int pid)
{
        if ((pid == 111) || (pid == 221) || (pid == 310) || (pid == 130)) {
                return pid;
        } else {
                return -pid;
        }
}

static enum alouette_return channel_parse(int * mode_ptr, int * sub_ptr)
{
        int mode, sub;
        if (*mode_ptr >= 100) {
                mode = *mode_ptr / 100;
                sub = *mode_ptr - 100 * mode;
        } else {
                mode = *mode_ptr;
                sub = 0;
        }
        *sub_ptr = sub;

        if ((mode < 0) || (mode > _channels.n)) {
            goto error;
        }

        if (sub) {
                if ((mode == 5) || (mode == 16) || (mode == 19) ||
                    (mode == 22)) {
                        if ((sub != 1) && (sub != 2)) {
                                goto error;
                        }
                } else if ((mode == 7) || (mode == 15)) {
                        if ((sub < 1) || (sub > 3)) {
                                goto error;
                        }
                } else {
                        goto error;
                }
        }

        *mode_ptr = mode;

        return ALOUETTE_RETURN_SUCCESS;

error:
        return message_error(ALOUETTE_RETURN_VALUE_ERROR,
            "bad decay mode (%d)", *mode_ptr);
}

/* Low level routine for configuring decay channels. */
static void channel_configure(int mode, int sub)
{
        memcpy(&tauola_taukle, &_channels.default_taukle,
            sizeof(struct tauola_taukle));

        if (sub) {
                if (mode == 5) {
                        if (sub == 1) {
                                tauola_taukle.bra1 = 1.;
                        } else {
                                tauola_taukle.bra1 = 0.;
                        }
                } else if (mode == 7) {
                        if (sub == 1) {
                                tauola_taukle.brks = 1.;
                                tauola_taukle.brk0 = 1.;
                                tauola_taukle.brk0b = 1.;
                        } else if (sub == 2) {
                                tauola_taukle.brks = 1.;
                                tauola_taukle.brk0 = 0.;
                                tauola_taukle.brk0b = 0.;
                        } else {
                                tauola_taukle.brks = 0.;
                        }
                } else if (mode == 15) {
                        if (sub == 1) {
                                tauola_taukle.brk0 = 1.;
                                tauola_taukle.brk0b = 1.;
                        } else if (sub == 2) {
                                const double p0 =
                                    _channels.default_taukle.brk0 *
                                    (1. - _channels.default_taukle.brk0b);
                                const double p1 =
                                    (1. - _channels.default_taukle.brk0) *
                                    _channels.default_taukle.brk0b;
                                if (alouette_random() * (p0 + p1) <= p0) {
                                        tauola_taukle.brk0 = 1.;
                                        tauola_taukle.brk0b = 0.;
                                } else {
                                        tauola_taukle.brk0 = 0.;
                                        tauola_taukle.brk0b = 1.;
                                }
                        } else {
                                tauola_taukle.brk0 = 0.;
                                tauola_taukle.brk0b = 0.;
                        }
                } else if ((mode == 16) || (mode == 19) || (mode == 22)) {
                        if (sub == 1) {
                                tauola_taukle.brk0 = 1.;
                                tauola_taukle.brk0b = 1.;
                        } else {
                                tauola_taukle.brk0 = 0.;
                                tauola_taukle.brk0b = 0.;
                        }
                }
        }
}

/* Select a decay channel in forward Monte Carlo. */
static enum alouette_return channel_select_forward(int * mode_ptr)
{
        enum alouette_return rc;
        int sub;
        if ((rc = channel_parse(mode_ptr, &sub)) == ALOUETTE_RETURN_SUCCESS) {
                channel_configure(*mode_ptr, sub);
        }
        return rc;
}

/* Select a decay channel for the BMC procedure. */
static enum alouette_return channel_select_backward(int daughter,
    int * mode_ptr, int * mother_ptr, int * multiplicity_ptr,
    double * weight_ptr)
{
        int converse = daughter_converse(daughter);
        int id = daughter_index(daughter);
        int ic = daughter_index(converse);

        if ((id < 0) && (ic < 0)) {
                return message_error(ALOUETTE_RETURN_VALUE_ERROR,
                    "bad daughter pid (%d)", daughter);
        }

        int mode = *mode_ptr, sub;
        enum alouette_return rc;
        if ((rc = channel_parse(&mode, &sub)) != ALOUETTE_RETURN_SUCCESS) {
                return rc;
        }

        int candidates[N_CHANNELS];
        int n_candidates = 0;
        if (mode == 0) {
                int i;
                for (i = 0; i < N_CHANNELS; i++) {
                        candidates[i] = i;
                }
                n_candidates = N_CHANNELS;
        } else if (mode == 5) {
                if ((!sub) || (sub == 1)) candidates[n_candidates++] = 4;
                if ((!sub) || (sub == 2)) candidates[n_candidates++] = 22;
        } else if (mode == 7) {
                if ((!sub) || (sub == 1)) candidates[n_candidates++] = 6;
                if ((!sub) || (sub == 2)) candidates[n_candidates++] = 23;
                if ((!sub) || (sub == 3)) candidates[n_candidates++] = 24;
        } else if (mode == 15) {
                if ((!sub) || (sub == 1)) candidates[n_candidates++] = 14;
                if ((!sub) || (sub == 2)) candidates[n_candidates++] = 25;
                if ((!sub) || (sub == 3)) candidates[n_candidates++] = 26;
        } else if (mode == 16) {
                if ((!sub) || (sub == 1)) candidates[n_candidates++] = 15;
                if ((!sub) || (sub == 2)) candidates[n_candidates++] = 27;
        } else if (mode == 19) {
                if ((!sub) || (sub == 1)) candidates[n_candidates++] = 18;
                if ((!sub) || (sub == 2)) candidates[n_candidates++] = 28;
        } else if (mode == 22) {
                if ((!sub) || (sub == 1)) candidates[n_candidates++] = 21;
                if ((!sub) || (sub == 2)) candidates[n_candidates++] = 29;
        } else {
                candidates[0] = mode - 1;
                n_candidates = 1;
        }

        /* Compute the total weight */
        int mother = *mother_ptr;
        double total_weight = 0.;
        if ((id >= 0) && ((mother == 0) || (mother == 15))) {
                int i;
                for (i = 0; i < n_candidates; i++) {
                        const int ii = candidates[i];
                        total_weight += _channels.weight[ii][id];
                }
        }
        if ((ic >= 0) && ((mother == 0) || (mother == -15))) {
                int i;
                for (i = 0; i < n_candidates; i++) {
                        const int ii = candidates[i];
                        total_weight += _channels.weight[ii][ic];
                }
        }

        if (total_weight == 0.) {
                return message_error(ALOUETTE_RETURN_VALUE_ERROR,
                    "inconsistent values for mother (%d), daughter (%d) "
                    "and mode (%d)", mother, daughter, mode);
        }

        /* Select the decay channel */
        double r;
        if (n_candidates == 1) {
                r = 0.;
        } else {
                r = alouette_random();
                if (r < 0) r = 0.;
                else if (r > 1) r = 1.;
                r *= total_weight;
        }

        double w = 0.;
        int channel = -1, multiplicity = 0;
        if ((id >= 0) && ((mother == 0) || (mother == 15))) {
                int i;
                for (i = 0; i < n_candidates; i++) {
                        const int ii = candidates[i];
                        w += _channels.weight[ii][id];
                        if (r <= w) {
                                *mother_ptr = 15;
                                channel = ii;
                                multiplicity =
                                     _channels.multiplicity[ii][id];
                                break;
                        }
                }
        }
        if ((channel == -1) && (ic >= 0) &&
            ((mother == 0) || (mother == -15))) {
                int i;
                for (i = 0; i < n_candidates; i++) {
                        const int ii = candidates[i];
                        w += _channels.weight[ii][ic];
                        if (r <= w) {
                                *mother_ptr = -15;
                                channel = ii;
                                multiplicity =
                                     _channels.multiplicity[ii][ic];
                                break;
                        }
                }
        }
        if (channel == -1) {
                /* This should not happen */
                return message_error(ALOUETTE_RETURN_VALUE_ERROR,
                    "unexpected channel (mother = %d, daughter = %d "
                    "and mode = %d)", mother, daughter, mode);
        }

        *mode_ptr = _channels.mode[channel];
        sub = _channels.subchannel[channel];
        channel_configure(*mode_ptr, sub);
        *multiplicity_ptr = multiplicity;

        double norm = 0.;
        int i;
        for (i = 0; i < n_candidates; i++) {
                const int ii = candidates[i];
                norm += _channels.branching_ratio[ii];
        }
        *weight_ptr = total_weight / (multiplicity * norm);

        return ALOUETTE_RETURN_SUCCESS;
}

/* Status flag for TAUOLA initialisation */
static int _tauola_initialised = 0;

/* Initialise TAUOLA and its wrapper. */
enum alouette_return alouette_initialise(double * xk0dec)
{
        message_reset();
        if (_tauola_initialised) {
                return message_error(ALOUETTE_RETURN_TAUOLA_ERROR,
                    "TAUOLA already initialised");
        }

#ifdef ALOUETTE_DEBUG
        /* Enable floating point exceptions. */
        feclearexcept(FE_ALL_EXCEPT);
        feenableexcept(
            FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW | FE_UNDERFLOW);
#endif

        /* Register a rally point in case of a TAUOLA error. */
        if (setjmp(alouette_context) != 0) {
                return message_error(ALOUETTE_RETURN_TAUOLA_ERROR, NULL);
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
                tauola_taurad.xk0dec = 1E-03;
        } else if (*xk0dec > 0.) {
                tauola_taurad.itdkrc = 1;
                tauola_taurad.xk0dec = *xk0dec;
        } else {
                tauola_taurad.itdkrc = 0;
                tauola_taurad.xk0dec = 1E-03; /* Setting to 0 results in FP
                                               * errors in TAUOLA fortran.
                                               */
        }

        /* PDG identifier for tau. */
        tauola_idfc.idff = 15;

        /* Configure the positions of taus in TAUOLA's stack. */
        extern struct {
                int npa, npb;
        } tauola_taupos;
        tauola_taupos.npa = 1;
        tauola_taupos.npb = 1;

        /* Initialise the built-in PRNG with an arbitrary seed for
         * TAUOLA warmup. First, backup the current random state.
         */
        random_cb * random = alouette_random;
        alouette_random = &random_uniform01;
        struct random_stream tmp_random;
        memcpy(&tmp_random, &_random_stream, sizeof tmp_random);

        unsigned long tmp_seed = 1357894;
        alouette_random_set(&tmp_seed);

        /* Initialise the decay routine. */
        double polarimeter[4];
        int state = -1;
        tauola_jaki.jak1 = 0;
        tauola_jaki.jak2 = 0;
        tauola_decay(&state, polarimeter);

        /* Restore the random state. */
        memcpy(&_random_stream, &tmp_random, sizeof tmp_random);
        alouette_random = random;

        /* Compute the branching ratios from TAUOLA's partial widths */
        extern struct {
                float gamprt[TAUOLA_MAX_CHANNELS];
                int jlist[TAUOLA_MAX_CHANNELS];
                int nchan;
        } tauola_taubra;

        _channels.n = tauola_taubra.nchan;
        int i;
        for (i = 0; i < _channels.n; i++) {
                const double d = tauola_taubra.gamprt[i];
                _channels.total_width += d;
                _channels.branching_ratio[i] = d;
        }
        for (i = 0; i < _channels.n; i++) {
                _channels.branching_ratio[i] /= _channels.total_width;
        }

        /* Backup the default BRs for subchannels. */
        memcpy(&_channels.default_taukle, &tauola_taukle,
            sizeof(struct tauola_taukle));

        /* Set the daughter multiplicities, needed by the BMC algorithm.
         *
         * This is hardcoded, since TAUOLA does not seem to explicitly export
         * this information.
         */
        if (tauola_taubra.nchan != 22) {
                /* We expect the `new-currents` version of `tauola-fortran` to
                 * be used. It includes 22 decay modes, with 4 channels
                 * exhibiting sub-channels, i.e. for modes 5, 7, 15 and 22.
                 */
                return message_error(ALOUETTE_RETURN_TAUOLA_ERROR,
                    "bad TAUOLA version");
        }

        /* 1st mode: tau- -> nu_tau e- nu_e_bar */
        _channels.mode[0] = 1;
        const int i_nutau = daughter_index(16);
        const int i_e = daughter_index(11);
        const int i_nueb = daughter_index(-12);
        _channels.multiplicity[0][i_nutau] = 1;
        _channels.multiplicity[0][i_e] = 1;
        _channels.multiplicity[0][i_nueb] = 1;

        /* 2nd mode: tau- -> nu_tau mu- nu_mu_bar */
        _channels.mode[1] = 2;
        const int i_mu = daughter_index(13);
        const int i_numub = daughter_index(-14);
        _channels.multiplicity[1][i_nutau] = 1;
        _channels.multiplicity[1][i_mu] = 1;
        _channels.multiplicity[1][i_numub] = 1;

        /* 3rd mode: tau- -> nu_tau pi- */
        _channels.mode[2] = 3;
        const int i_pim = daughter_index(-211);
        _channels.multiplicity[2][i_nutau] = 1;
        _channels.multiplicity[2][i_pim] = 1;

        /* 4th mode: tau- -> nu_tau pi- pi0 */
        _channels.mode[3] = 4;
        const int i_pi0 = daughter_index(111);
        _channels.multiplicity[3][i_nutau] = 1;
        _channels.multiplicity[3][i_pim] = 1;
        _channels.multiplicity[3][i_pi0] = 1;

        /* 5th mode: tau- -> nu_tau a1 (2 sub-channels, 501 and 502) */
        _channels.mode[4] = 5;
        _channels.subchannel[4] = 1;
        const int i_pip = daughter_index(211);
        _channels.multiplicity[4][i_nutau] = 1;
        _channels.multiplicity[4][i_pim] = 2;
        _channels.multiplicity[4][i_pip] = 1;

        _channels.mode[22] = 5;
        _channels.subchannel[22] = 2;
        _channels.multiplicity[22][i_nutau] = 1;
        _channels.multiplicity[22][i_pi0] = 2;
        _channels.multiplicity[22][i_pim] = 1;
        _channels.branching_ratio[22] = _channels.branching_ratio[4] *
            (1. - _channels.default_taukle.bra1);

        _channels.branching_ratio[4] *= _channels.default_taukle.bra1;

        /* 6th mode: tau- -> nu_tau K- */
        _channels.mode[5] = 6;
        const int i_km = daughter_index(-321);
        _channels.multiplicity[5][i_nutau] = 1;
        _channels.multiplicity[5][i_km] = 1;

        /* 7th mode: tau- -> nu_tau K*- (3 sub-channels, 701, 702 and 703) */
        _channels.mode[6] = 7;
        _channels.subchannel[6] = 1;
        const int i_ks = daughter_index(310);
        _channels.multiplicity[6][i_nutau] = 1;
        _channels.multiplicity[6][i_pim] = 1;
        _channels.multiplicity[6][i_ks] = 1;

        _channels.mode[23] = 7;
        _channels.subchannel[23] = 2;
        const int i_kl = daughter_index(130);
        _channels.multiplicity[23][i_nutau] = 1;
        _channels.multiplicity[23][i_pim] = 1;
        _channels.multiplicity[23][i_kl] = 1;
        _channels.branching_ratio[23] = _channels.branching_ratio[6] *
            _channels.default_taukle.brks *
            (1. - _channels.default_taukle.brk0);

        _channels.mode[24] = 7;
        _channels.subchannel[24] = 3;
        _channels.multiplicity[24][i_nutau] = 1;
        _channels.multiplicity[24][i_pi0] = 1;
        _channels.multiplicity[24][i_km] = 1;
        _channels.branching_ratio[24] = _channels.branching_ratio[6] *
            (1. - _channels.default_taukle.brks);

        _channels.branching_ratio[6] *=
            _channels.default_taukle.brks * _channels.default_taukle.brk0;

        /* 8th mode: tau- -> nu_tau 2 pi- pi+ pi0 */
        _channels.mode[7] = 8;
        _channels.multiplicity[7][i_nutau] = 1;
        _channels.multiplicity[7][i_pim] = 2;
        _channels.multiplicity[7][i_pip] = 1;
        _channels.multiplicity[7][i_pi0] = 1;

        /* 9th mode: tau- -> nu_tau 3 pi0 pi- */
        _channels.mode[8] = 9;
        _channels.multiplicity[8][i_nutau] = 1;
        _channels.multiplicity[8][i_pi0] = 3;
        _channels.multiplicity[8][i_pim] = 1;

        /* 10th mode: tau- -> nu_tau 2 pi- pi+ 2 pi0 */
        _channels.mode[9] = 10;
        _channels.multiplicity[9][i_nutau] = 1;
        _channels.multiplicity[9][i_pim] = 2;
        _channels.multiplicity[9][i_pip] = 1;
        _channels.multiplicity[9][i_pi0] = 2;

        /* 11th mode: tau- -> nu_tau 3 pi- 2 pi+ */
        _channels.mode[10] = 11;
        _channels.multiplicity[10][i_nutau] = 1;
        _channels.multiplicity[10][i_pim] = 3;
        _channels.multiplicity[10][i_pip] = 2;

        /* 12th mode: tau- -> nu_tau 3 pi- 2 pi+ pi0 */
        _channels.mode[11] = 12;
        _channels.multiplicity[11][i_nutau] = 1;
        _channels.multiplicity[11][i_pim] = 3;
        _channels.multiplicity[11][i_pip] = 2;
        _channels.multiplicity[11][i_pi0] = 1;

        /* 13th mode: tau- -> nu_tau 2 pi- 1 pi+ 3 pi0 */
        _channels.mode[12] = 13;
        _channels.multiplicity[12][i_nutau] = 1;
        _channels.multiplicity[12][i_pim] = 2;
        _channels.multiplicity[12][i_pip] = 1;
        _channels.multiplicity[12][i_pi0] = 3;

        /* 14th mode: tau- -> nu_tau K- pi- K+ */
        _channels.mode[13] = 14;
        const int i_kp = daughter_index(321);
        _channels.multiplicity[13][i_nutau] = 1;
        _channels.multiplicity[13][i_km] = 1;
        _channels.multiplicity[13][i_pim] = 1;
        _channels.multiplicity[13][i_kp] = 1;

        /* 15th mode: tau- -> nu_tau K0 pi- K0b (3 sub-channels, 1501, 1502 and
         * 1503)
         */
        _channels.mode[14] = 15;
        _channels.subchannel[14] = 1;
        _channels.multiplicity[14][i_nutau] = 1;
        _channels.multiplicity[14][i_ks] = 2;
        _channels.multiplicity[14][i_pim] = 1;

        _channels.mode[25] = 15;
        _channels.subchannel[25] = 2;
        _channels.multiplicity[25][i_nutau] = 1;
        _channels.multiplicity[25][i_ks] = 1;
        _channels.multiplicity[25][i_pim] = 1;
        _channels.multiplicity[25][i_kl] = 1;
        _channels.branching_ratio[25] = _channels.branching_ratio[14] * (
            _channels.default_taukle.brk0 *
            (1. - _channels.default_taukle.brk0b) +
            (1. - _channels.default_taukle.brk0) *
            _channels.default_taukle.brk0b);

        _channels.mode[26] = 15;
        _channels.subchannel[26] = 3;
        _channels.multiplicity[26][i_nutau] = 1;
        _channels.multiplicity[26][i_kl] = 2;
        _channels.multiplicity[26][i_pim] = 1;
        _channels.branching_ratio[26] = _channels.branching_ratio[14] *
            (1. - _channels.default_taukle.brk0) *
            (1. - _channels.default_taukle.brk0b);

        _channels.branching_ratio[14] *=
            _channels.default_taukle.brk0 * _channels.default_taukle.brk0b;

        /* 16th mode: tau- -> nu_tau K- pi0 K0 (2 sub-channels, 1601 and
         * 1602)
         */
        _channels.mode[15] = 16;
        _channels.subchannel[15] = 1;
        _channels.multiplicity[15][i_nutau] = 1;
        _channels.multiplicity[15][i_km] = 1;
        _channels.multiplicity[15][i_pi0] = 1;
        _channels.multiplicity[15][i_ks] = 1;

        _channels.mode[27] = 16;
        _channels.subchannel[27] = 2;
        _channels.multiplicity[27][i_nutau] = 1;
        _channels.multiplicity[27][i_km] = 1;
        _channels.multiplicity[27][i_pi0] = 1;
        _channels.multiplicity[27][i_kl] = 1;
        _channels.branching_ratio[27] = _channels.branching_ratio[15] *
            (1. - _channels.default_taukle.brk0);

        _channels.branching_ratio[15] *= _channels.default_taukle.brk0;

        /* 17th mode: tau- -> nu_tau 2 pi0 K- */
        _channels.mode[16] = 17;
        _channels.multiplicity[16][i_nutau] = 1;
        _channels.multiplicity[16][i_pi0] = 2;
        _channels.multiplicity[16][i_km] = 1;

        /* 18th mode: tau- -> nu_tau 2 K- pi- pi+ */
        _channels.mode[17] = 18;
        _channels.multiplicity[17][i_nutau] = 1;
        _channels.multiplicity[17][i_km] = 1;
        _channels.multiplicity[17][i_pim] = 1;
        _channels.multiplicity[17][i_pip] = 1;

        /* 19th mode: tau- -> nu_tau pi- K0 pi0 (2 sub-channels, 1901 and
         * 1902)
         */
        _channels.mode[18] = 19;
        _channels.subchannel[18] = 1;
        _channels.multiplicity[18][i_nutau] = 1;
        _channels.multiplicity[18][i_pim] = 1;
        _channels.multiplicity[18][i_ks] = 1;
        _channels.multiplicity[18][i_pi0] = 1;

        _channels.mode[28] = 19;
        _channels.subchannel[28] = 2;
        _channels.multiplicity[28][i_nutau] = 1;
        _channels.multiplicity[28][i_pim] = 1;
        _channels.multiplicity[28][i_kl] = 1;
        _channels.multiplicity[28][i_pi0] = 1;
        _channels.branching_ratio[28] = _channels.branching_ratio[18] *
            (1. - _channels.default_taukle.brk0);

        _channels.branching_ratio[18] *= _channels.default_taukle.brk0;

        /* 20th mode: tau- -> nu_tau eta pi- pi0 */
        _channels.mode[19] = 20;
        const int i_eta = daughter_index(221);
        _channels.multiplicity[19][i_nutau] = 1;
        _channels.multiplicity[19][i_eta] = 1;
        _channels.multiplicity[19][i_pim] = 1;
        _channels.multiplicity[19][i_pi0] = 1;

        /* 21st mode: tau- -> nu_tau pi- pi0 gamma */
        _channels.mode[20] = 21;
        _channels.multiplicity[20][i_nutau] = 1;
        _channels.multiplicity[20][i_pim] = 1;
        _channels.multiplicity[20][i_pi0] = 1;

        /* 22nd mode: tau- -> nu_tau K- K0 (2 sub-channels, 2201 and 2202 */
        _channels.mode[21] = 22;
        _channels.subchannel[21] = 1;
        _channels.multiplicity[21][i_nutau] = 1;
        _channels.multiplicity[21][i_km] = 1;
        _channels.multiplicity[21][i_ks] = 1;

        _channels.mode[29] = 22;
        _channels.subchannel[29] = 2;
        _channels.multiplicity[29][i_nutau] = 1;
        _channels.multiplicity[29][i_km] = 1;
        _channels.multiplicity[29][i_kl] = 1;
        _channels.branching_ratio[29] = _channels.branching_ratio[21] *
            (1. - _channels.default_taukle.brk0);

        _channels.branching_ratio[21] *= _channels.default_taukle.brk0;

        /* Compute the BMC sampling weights. */
        for (i = 0; i < N_CHANNELS; i++) {
                int j;
                for (j = 0; j < N_DAUGHTERS; j++) {
                        const double tmp = _channels.multiplicity[i][j] *
                            _channels.branching_ratio[i];
                        _channels.weight[i][j] = tmp;
                        _channels.total_weight[j] += tmp;
                }
        }

        /* Flag as initialised and return. */
        _tauola_initialised = 1;

        return ALOUETTE_RETURN_SUCCESS;
}
#undef TAUOLA_MAX_CHANNELS
#undef N_CHANNELS
#undef N_DAUGHTERS

/* Get the last (error) message(s). */
const char * alouette_message(void)
{
        static const char * msg[ALOUETTE_N_RETURNS] = {
            "Operation succeeded", "A value is out of range",
            "A Tauola error occured" };

        if (_message.data[0] == 0x0) {
                if (_message.code == ALOUETTE_RETURN_SUCCESS) {
                        return NULL;
                } else {
                        return msg[_message.code];
                }
        } else {
                return _message.data;
        }
}

/* Low level decay routine with TAUOLA. */
static enum alouette_return decay(
    int pid, int mode, struct alouette_products * products)
{
        /* Register a rally point in case of a TAUOLA error. */
        if (setjmp(alouette_context) != 0) {
                return message_error(ALOUETTE_RETURN_TAUOLA_ERROR, NULL);
        }

        tauola_jaki.jak1 = mode;
        tauola_jaki.jak2 = mode;
        _products = products;

        /* Call TAUOLA's decay routine. */
        int type = (pid > 0) ? 1 : 2;
        double h[4];
        tauola_decay(&type, h);
        type += 10;
        tauola_decay(&type, h);
        memcpy(_products->polarimeter, h, sizeof _products->polarimeter);

        return ALOUETTE_RETURN_SUCCESS;
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
static int rotate_direction(double cos_theta, double * direction)
{
        /* Check the numerical sine. */
        const double stsq = 1. - cos_theta * cos_theta;
        if (stsq <= DBL_EPSILON) {
                return EXIT_FAILURE;
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
        const double phi = M_PI * (1. - 2. * alouette_random());
        const double cp = cos(phi);
        const double sp = sin(phi);
        direction[0] = cos_theta * direction[0] + st * (cp * u0x + sp * u1x);
        direction[1] = cos_theta * direction[1] + st * (cp * u0y + sp * u1y);
        direction[2] = cos_theta * direction[2] + st * (cp * u0z + sp * u1z);

        return EXIT_SUCCESS;
}

/* Build a rotation matrix from vi to vf. */
static int build_rotation(
    const double * vi, const double * vf, double norm_f, double * R)
{
        /* Build the rotation axis. */
        double n[3] = { vi[1] * vf[2] - vi[2] * vf[1],
                vi[2] * vf[0] - vi[0] * vf[2], vi[0] * vf[1] - vi[1] * vf[0] };
        double nrm = n[0] * n[0] + n[1] * n[1] + n[2] * n[2];
        if (fabs(nrm) <= FLT_EPSILON) {
                if (vi[0] * vf[0] + vi[1] * vf[1] + vi[2] * vf[2] > 0.) {
                        /* The two vectors are aligned within numeric errors. */
                        return EXIT_FAILURE;
                } else {
                        /* The two vectors are back to back. */
                        R[0] = -1.;
                        R[1] = 0.;
                        R[2] = 0.;
                        R[3] = 0.;
                        R[4] = -1.;
                        R[5] = 0.;
                        R[6] = 0.;
                        R[7] = 0.;
                        R[8] = -1.;

                        return EXIT_SUCCESS;
                }
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

        return EXIT_SUCCESS;
}

/* Decay a tau with TAUOLA. */
enum alouette_return alouette_decay(int mode, int pid, const double momentum[3],
    const double * polarisation, struct alouette_products * products)
{
        /* Initialise the products container. */
        products_reset(products);

        /* Initialise the library, if not already done. */
        enum alouette_return rc;
        if (!_tauola_initialised) {
                if ((rc = alouette_initialise(NULL)) !=
                    ALOUETTE_RETURN_SUCCESS) {
                        return rc;
                }
        }
        message_reset();

        /* Parse and check the decay mode. */
        if ((rc = channel_select_forward(&mode)) != ALOUETTE_RETURN_SUCCESS) {
                return rc;
        }

        /* Check the mother pid */
        if (abs(pid) != 15) {
                return message_error(ALOUETTE_RETURN_VALUE_ERROR,
                    "bad mother pid (%d)", pid);
        }

        /* Decay a tau in its rest frame. */
        for (;;) {
                if ((rc = decay(pid, mode, products)) !=
                    ALOUETTE_RETURN_SUCCESS)
                        return rc;
                if ((!isnan(products->polarimeter[0])) &&
                    (!isinf(products->P[0][0]))) {
                        break;
                } else {
                        products_reset(products);
                }
        }

        if (polarisation != NULL) {
                const double * const hi = products->polarimeter;
                const double h2 =
                    hi[0] * hi[0] + hi[1] * hi[1] + hi[2] * hi[2];
                const double s2 = polarisation[0] * polarisation[0] +
                    polarisation[1] * polarisation[1] +
                    polarisation[2] * polarisation[2];
                if ((h2 > FLT_EPSILON) && (s2 > FLT_EPSILON)) {
                        /* Draw the direction of the polarimeter according to
                         * the spin polarisation.
                         */
                        double s = sqrt(s2);
                        if (s > 1.) s = 1.;
                        const double z = alouette_random();
                        double delta2 = 4 * s * z + (s - 1.) * (s - 1.);
                        if (delta2 <= FLT_EPSILON) delta2 = 0.;
                        double cos_theta = (sqrt(delta2) - 1.) / s;
                        if (cos_theta < -1.) cos_theta = -1.;
                        else if (cos_theta > 1.) cos_theta = 1.;

                        double u[3] = {polarisation[0] / s,
                            polarisation[1] / s, polarisation[2] / s};
                        if (rotate_direction(cos_theta, u) == EXIT_SUCCESS) {
                                /* Update the direction of decay products in
                                 * order to match the new polarimeter.
                                 */
                                const double h = sqrt(h2);
                                double R[9];
                                if (build_rotation(u, hi, h, R) ==
                                    EXIT_SUCCESS) {
                                        int i;
                                        for (i = 0; i < products->size; i++) {
                                                mv_multiply(
                                                    R, &products->P[i][0]);
                                        }

                                        for (i = 0; i < 3; i++) {
                                                products->polarimeter[i] =
                                                    h * u[i];
                                        }
                                }
                        }
                }
        }

        /* Boost the daughters to the lab frame. */
        const double amtau = tauola_parmas.amtau;
        const double tau[3] = { momentum[0] / amtau,
                momentum[1] / amtau, momentum[2] / amtau };
        const double gamma =
            sqrt(1. + tau[0] * tau[0] + tau[1] * tau[1] + tau[2] * tau[2]);
        products->weight = 1.;
        if (gamma <= 1. + FLT_EPSILON) {
                return ALOUETTE_RETURN_SUCCESS;
        }

        const double energy = gamma * amtau;
        double Psum[4] = {momentum[0], momentum[1], momentum[2], energy};
        int i;
        double * P;
        for (i = 0, P = &products->P[0][0]; i < products->size - 1;
            i++, P += 4) {
                /* Apply the boost. */
                const double ptau =
                    P[0] * tau[0] + P[1] * tau[1] + P[2] * tau[2];
                const double tmp = ptau / (gamma + 1.) + P[3];
                P[0] += tmp * tau[0];
                P[1] += tmp * tau[1];
                P[2] += tmp * tau[2];
                P[3] = gamma * P[3] + ptau;

                int j;
                for (j = 0; j < 4; j++) Psum[j] -= P[j];
        }

        /* Set last daughter using energy-momentum conservation. */
        int j;
        for (j = 0; j < 4; j++) P[j] = Psum[j];

        return ALOUETTE_RETURN_SUCCESS;
}

/* Mother particle(s) for backward decays. */
int alouette_undecay_mother = 0;


/** Tuning parameter for the spin bias in backward decays. */
double alouette_undecay_bias = 1.;

/* Backward decay from a tau neutrino to a tau. */
enum alouette_return alouette_undecay(int mode, int daughter,
    const double momentum[3], alouette_polarisation_cb * polarisation_cb,
    struct alouette_products * products)
{
        /* Initialise the products container. */
        products_reset(products);

        /* Initialise the library, if not already done. */
        enum alouette_return rc;
        if (!_tauola_initialised) {
                if ((rc = alouette_initialise(NULL)) !=
                    ALOUETTE_RETURN_SUCCESS) {
                        return rc;
                }
        }
        message_reset();

        /* Check zero momentum */
        const double momentum2 = momentum[0] * momentum[0] +
            momentum[1] * momentum[1] + momentum[2] * momentum[2];
        if (momentum2 < FLT_EPSILON) {
                return message_error(ALOUETTE_RETURN_VALUE_ERROR,
                    "bad daughter momentum (%g)",
                    momentum2);
        }

        /* Check the daughter pid. */
        const int valid_pids[] = {-16, 16, -11, 11, -12, 12, -13, 13, -14, 14,
            -211, 211, 111, 221, -321, 321, 310, 130};
        int is_valid = 0, i;
        for (i = 0; i < sizeof(valid_pids) / sizeof(valid_pids[0]); i++) {
                if (daughter == valid_pids[i]) {
                        is_valid = 1;
                        break;
                }
        }
        if (!is_valid) {
                return message_error(ALOUETTE_RETURN_VALUE_ERROR,
                        "bad daughter pid (%d)", daughter);
        }

        /* Check the mother PID. */
        int mother = alouette_undecay_mother;
        if (mother && (abs(mother) != 15)) {
                return message_error(ALOUETTE_RETURN_VALUE_ERROR,
                        "bad mother pid (%d)", mother);
        }

        /* Check the bias value */
        double bias = alouette_undecay_bias;
        if ((bias < -1.) || (bias > 1.)) {
                return message_error(ALOUETTE_RETURN_VALUE_ERROR,
                        "bad bias value (%g)", bias);
        }

        /* Select the decay mode. */
        double weight = 1.;
        int multiplicity = 1;
        if ((rc = channel_select_backward(
            daughter, &mode, &mother, &multiplicity, &weight)) !=
                ALOUETTE_RETURN_SUCCESS) {
                return rc;
        }

        /* Decay an unpolarised tau in its rest frame. */
        for (;;) {
                if ((rc = decay(mother, mode, products)) !=
                    ALOUETTE_RETURN_SUCCESS)
                        return rc;
                if ((!isnan(products->polarimeter[0])) &&
                    (!isinf(products->P[0][0]))) {
                        break;
                } else {
                        products_reset(products);
                }
        }

        /* Get the daughter's index. */
        int idx;
        if (multiplicity > 1) {
                idx = (int)(multiplicity * alouette_random());
                if (idx < 0) idx = 0;
                else if (idx >= multiplicity) idx = multiplicity - 1;
        } else {
                idx = 0;
        }
        double * pi;
        int idau;
        for (idau = 0, pi = &products->P[0][0]; idau < products->size;
            idau++, pi += 4) {
                if (products->pid[idau] == daughter) {
                        if (idx-- == 0) break;
                }
        }
        if (idau == products->size) {
                return message_error(ALOUETTE_RETURN_VALUE_ERROR,
                    "no such daughter (%d) in decay products", daughter);
        }

        if (polarisation_cb != NULL) {
                if (mother > 0) bias = -bias;

                const double * const hi = products->polarimeter;
                const double h2 =
                    hi[0] * hi[0] + hi[1] * hi[1] + hi[2] * hi[2];
                const double s = fabs(bias);
                if ((h2 > FLT_EPSILON) && (s > FLT_EPSILON)) {
                        /* Draw the direction of the polarimeter according to
                         * the daughter momentum and the bias polarisation.
                         */
                        const double z = alouette_random();
                        double delta2 = 4 * s * z + (s - 1.) * (s - 1.);
                        if (delta2 <= FLT_EPSILON) delta2 = 0.;
                        double cos_theta = (sqrt(delta2) - 1.) / s;
                        if (cos_theta < -1.) cos_theta = -1.;
                        else if (cos_theta > 1.) cos_theta = 1.;
                        weight /= 1. + s * cos_theta;

                        double norm = bias / sqrt(momentum2);
                        double u[3] = {momentum[0] * norm,
                            momentum[1] * norm, momentum[2] * norm};
                        if (rotate_direction(cos_theta, u) == EXIT_SUCCESS) {
                                /* Update the direction of decay products in
                                 * order to match the new polarimeter.
                                 */
                                const double h = sqrt(h2);
                                double R[9];
                                if (build_rotation(u, hi, h, R) ==
                                    EXIT_SUCCESS) {
                                        int i;
                                        for (i = 0; i < products->size; i++) {
                                                mv_multiply(R,
                                                    &products->P[i][0]);
                                        }

                                        for (i = 0; i < 3; i++) {
                                                products->polarimeter[i] =
                                                    h * u[i];
                                        }
                                }
                        }
                }
        }

        /* Compute the parameters of the frame transform. */
        const double mass2 = pi[3] * pi[3] - (pi[0] * pi[0] + pi[1] * pi[1] +
            pi[2] * pi[2]);
        const double mass = (mass2 > FLT_EPSILON) ? sqrt(mass2) : 0.;
        const double energy = sqrt(momentum[0] * momentum[0] +
            momentum[1] * momentum[1] + momentum[2] * momentum[2] +
            mass * mass);
        double beta[3] = {momentum[0] - pi[0], momentum[1] - pi[1],
            momentum[2] - pi[2]};
        const double gamma = 1. + (beta[0] * beta[0] + beta[1] * beta[1] +
            beta[2] * beta[2]) / (momentum[0] * pi[0] + momentum[1] * pi[1] +
            momentum[2] * pi[2] + energy * pi[3] + mass * mass);
        const double t0 = (gamma + 1.) / (gamma * (energy + pi[3]));
        int j;
        for (j = 0; j < 3; j++) {
                beta[j] *= -t0;
        }

        /* Boost the daughters to the lab frame. */
        double Et = 0., Pt[3] = { 0., 0., 0. };
        double * pj;
        for (j = 0, pj = &products->P[0][0]; j < products->size; j++, pj += 4) {
                if (j == idau) {
                        /* Update with the provided daughter data. */
                        Et += energy;
                        Pt[0] += momentum[0];
                        Pt[1] += momentum[1];
                        Pt[2] += momentum[2];
                        continue;
                }

                /* Apply the boost and update the sum. */
                const double pbeta =
                    pj[0] * beta[0] + pj[1] * beta[1] + pj[2] * beta[2];
                const double tmp = gamma * (pbeta * gamma / (gamma + 1.)
                    - pj[3]);
                pj[0] += tmp * beta[0];
                pj[1] += tmp * beta[1];
                pj[2] += tmp * beta[2];
                pj[3] = gamma * (pj[3] - pbeta);
                Pt[0] += pj[0];
                Pt[1] += pj[1];
                Pt[2] += pj[2];
                Et += pj[3];
        }

        /* Set the backward Monte-Carlo weight. */
        const double g1 = gamma + 1.;
        const double ee = energy + pi[3];
        const double amtau = tauola_parmas.amtau;
        weight *= amtau * amtau * amtau * gamma * g1 * g1 / (ee * ee * energy);

        /* Weight for the spin polarisation of the tau mother. */
        if (polarisation_cb != NULL) {
                double pol[3];
                polarisation_cb(mother, Pt, pol);
                weight *= 1. + products->polarimeter[0] * pol[0] +
                    products->polarimeter[1] * pol[1] +
                    products->polarimeter[2] * pol[2];
        }

        /* Insert the tau mother at top of the stack and remove the initial
         * daughter. In order to ensure the conservation of the energy-momentum,
         * the mother's 4 momentum has been re-computed from the boosted
         * products.
         */
        for (j = idau; j > 0; j--) {
                products->pid[j] = products->pid[j - 1];
                memcpy(&products->P[j][0], &products->P[j - 1][0],
                    4 * sizeof(products->P[0][0]));
        }

        products->pid[0] = mother;
        products->P[0][0] = Pt[0];
        products->P[0][1] = Pt[1];
        products->P[0][2] = Pt[2];
        products->P[0][3] = Et;

        products->weight = weight;

        return ALOUETTE_RETURN_SUCCESS;
}
