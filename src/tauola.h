/*
 * Copyright (C) 2022 Université Clermont Auvergne, CNRS/IN2P3, LPC
 * Author: Valentin NIESS (niess@in2p3.fr)
 *
 * This software is a C-library compliant reformating of TAUOLA.
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

#ifndef TAUOLA_H
#define TAUOLA_H

#ifdef __cplusplus
extern "C" {
#endif

/* Decay main function.
 *
 * This function behaves like the DEKAY function from TAUOLA Fortran. E.g., it
 * must first be called with kto=-1 for initialisation.
 */
void tauola_decay(int * kt0, double hx[4]);

/* Decay products callback function, to be implemented by the user. */
void tauola_filhep(int * n, int * status, int * pid, int * mother_first,
    int * mother_last, int * daughter_first, int * daughter_last, float p4[4],
    float * p_inv_mass, int * photos_flag);

/* Messaging callback function, to be implemented by the user. */
void tauola_print(const char * msg);

/* Exit callback function, to be be implemented by the user.
 *
 * This function is called by TAUOLA when an error occurs. Note that it must
 * return to TAUOLA's calling context, not to TAUOLA itself. This can be done,
 * e.g. by using setjmp before calling `tauola_decay`, and then longjmp inside
 * `tauola_stop`.
 */
void tauola_stop(void);

/* PRNG stream callback function, to be implemented by the user.
 *
 * Receives a size n array of floats. At return, the array must have been
 * filled with pseudo-random numbers uniformly distributed in (0,1).
 */
void tauola_random(float * r, int * n);

/* TAUOLA mass data.
 *
 * These data are set at TAUOLA initialisation. They should be considered as
 * read-only.
 */
struct tauola_parmas {
    float amtau, amnuta, amel, amnue, ammu, amnumu, ampiz, ampi, amro, gamro,
          ama1, gama1, amk, amkz, amkst, gamkst;
};

extern struct tauola_parmas tauola_parmas;

/* TAUOLA weight data.
 *
 * This is an addition to the initial TAUOLA library. The structures below
 * contain the maximum weights ($W_{max}$) computed during TAUOLA's
 * initialisation procedure.
 */
struct tauola_weight {
        float wtmax;
};

extern struct tauola_weight tauola_weight_dadmaa;
extern struct tauola_weight tauola_weight_dadmel;
extern struct tauola_weight tauola_weight_dadmks;
extern struct tauola_weight tauola_weight_dadmmu;
extern struct tauola_weight tauola_weight_dadmro;

struct tauola_weight_dadnew {
        float wtmax[15];
};
extern struct tauola_weight_dadnew tauola_weight_dadnew;

#ifdef __cplusplus
}
#endif
#endif
