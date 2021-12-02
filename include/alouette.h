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

#ifndef ALOUETTE_H
#define ALOUETTE_H

#ifdef __cplusplus
extern "C" {
#endif

/** Return codes for the API functions. */
enum alouette_return {
        /** Execution was successful. */
        ALOUETTE_RETURN_SUCCESS = 0,
        /** A parameter is out of its validity range. */
        ALOUETTE_RETURN_VALUE_ERROR,
        /** A floating point error occured. */
        ALOUETTE_RETURN_FLOATING_ERROR,
        /** A TAUOLA error occured. */
        ALOUETTE_RETURN_TAUOLA_ERROR,
        /** The number of ALOUETTE return codes.  */
        ALOUETTE_N_RETURNS
};

/**
 * Initialise TAUOLA and the ALOUETTE wrapper.
 *
 * @param xk0dec  Factor for radiative corrections or `NULL`.
 * @return On success `ALOUETTE_RETURN_SUCCESS` is returned otherwise an error
 * code is returned as detailed below.
 *
 * This function allows to initialise the library with custom settings. Note
 * that if not done explictly, the library is initialised on need with default
 * settings. Thus, this function must be called before other library functions
 * if custom settings are desired.
 *
 * __Note__ : if *xk0dec* is `NULL`, then a default value of 1E-03 is used,
 * following Jezabek et al., CPC 70 (1992) 69-76.
 *
 * __Error codes__
 *
 *     ALOUETTE_RETURN_TAUOLA_ERROR    A TAUOLA error occured.
 */
enum alouette_return alouette_initialise(
    double * xk0dec);

/**
 * Return the last (error) message(s).
 *
 * This function returns a string containing the last (error) message(s) issued
 * by TAUOLA, or by the ALOUETTE wrapper.
 */
const char * alouette_message(void);


/**
 * Container for TAUOLA decay products.
 */
struct alouette_products {
#define ALOUETTE_MAX_SIZE 16
        /** Number of decay products. */
        int size;
        /** PDG IDs of decay products. */
        int pid[ALOUETTE_MAX_SIZE];
        /** Four momenta (px, py, pz, E) of decay products, in GeV/c. */
        double P[ALOUETTE_MAX_SIZE][4];
        /** Polarimeter vector of the decay. */
        double polarimeter[3];
        /** Monte Carlo weight of the decay. */
        double weight;
};

/**
 * Perform a forward Monte-Carlo tau decay.
 *
 * @param mode            The tau decay mode, or `0`.
 * @param pid             The PDG ID of the decaying tau, i.e. `15` or `-15`.
 * @param momentum        The tau momentum at decay, in GeV/c.
 * @param polarisation    The tau polarisation vector, or `NULL`.
 * @param products        The decay products.
 * @return On success `ALOUETTE_RETURN_SUCCESS` is returned otherwise an error
 * code is returned as detailed below.
 *
 * Simulate a tau decay with TAUOLA. An optionnal polarisation 3-vector can
 * be provided. If `ǸULL` spin effects are not simulated.
 *
 * __Note__ : tau decay modes are indexed according to TAUOLA. If set to zero,
 * then the decay mode is randomised over al available modes, according to
 * branching ratios.
 *
 * __Error codes__
 *
 *     ALOUETTE_RETURN_VALUE_ERROR       The provided *pid* is not valid.
 *
 *     ALOUETTE_RETURN_FLOATING_ERROR    A floating point error occured.
 *
 *     ALOUETTE_RETURN_TAUOLA_ERROR      A TAUOLA error occured.
 */
enum alouette_return alouette_decay(
    int mode,
    int pid,
    const double momentum[3],
    const double * polarisation,
    struct alouette_products * products);

/**
 * Callback for the tau polarisation in backward decays.
 *
 * @param pid             The PDG ID of the tau mother, i.e. `15` or `-15`.
 * @param momentum        The tau's momentum at decay, in GeV/c.
 * @param polarisation    The tau polarisation.
 *
 * In a backward decay, the spin polarisation of the tau mother is not known
 * a priori. The user can supply an a posteriori value with this callback.
 *
 * __Note__ : tau decay modes are indexed according to TAUOLA. If set to zero,
 * then the decay mode is randomised over al available modes, according to
 * branching ratios.
 */
typedef void alouette_polarisation_cb(
    int pid,
    const double momentum[3],
    double * polarisation);

/**
 * Perform a backward Monte-Carlo tau decay.
 *
 * @param mode            The tau decay mode, or `0` for any.
 * @param daughter        The PDG ID of the daughter decaying product,
 *                        e.g. `16` for a tau neutrino.
 * @param mother          The PDG ID of the mother particle,
 *                        or `0` for any valid mother.
 * @param momentum        The daughter momentum after decay, in GeV/c.
 * @param polarisation    A callback for the spin polarisation of the mother or
 *                        `NULL`.
 * @param bias            Tuning parameter for the spin bias.
 * @param weight          The backward Monte-Carlo decay products.
 * @return On success `ALOUETTE_RETURN_SUCCESS` is returned otherwise an error
 * code is returned as detailed below.
 *
 * Simulate a backward tau decay from a given decay product. The spin
 * polarisation of the primary tau can be provided a posteriori. Set
 * *polarisation* to `NULL` in order to ignore spin effects.
 *
 * The spin *bias* parameters allows to control the biasing of the angular
 * distribution of decay products in the mother's rest frame. It must be in the
 * range [-1, 1]. It is expected to be a hint on the mother's longitudinal spin
 * polarisation. Thus, typically -1 for a &tau;^(-) or +1 for a &tau;^(+). Set
 * the bias to zero if the spin polarization is a priori unknown.
 *
 * __Error codes__
 *
 *     ALOUETTE_RETURN_VALUE_ERROR       The provided *pid* is not valid.
 *
 *     ALOUETTE_RETURN_FLOATING_ERROR    A floating point error occured.
 *
 *     ALOUETTE_RETURN_TAUOLA_ERROR      A TAUOLA error occured.
 */
enum alouette_return alouette_undecay(
    int mode,
    int daughter,
    int mother,
    const double momentum[3],
    alouette_polarisation_cb * polarisation,
    double bias,
    struct alouette_products * products);

/**
 * The library PRNG, uniform over (0,1).
 *
 * The ALOUETTE library uses a single random stream, replacing TAUOLA's `RANMAR`
 * random engine with a Merseene Twister algorithm. The PRNG can be seeded with
 * the `alouette_random_seed_set` function. The current seed value is returned
 * by the `alouette_random_seed_get` function.
 *
 * __Note__: custom PRNGs can be used, instead of the built-in one, by
 * overriding the `alouette_random` pointer. The replacement must return pseudo
 * random numbers uniformly distributed over (0, 1).
 */
extern float (*alouette_random)(void);

/**
 * Get the random seed of the built-in PRNG.
 *
 * @return The PRNG random seed.
 */
unsigned long alouette_random_seed(void);

/**
 * (Re)set the built-in PRNG.
 *
 * @param seed    The PRNG random seed or `NULL`.
 *
 * Reset the built-in PRNG with a new random *seed*. If *seed* is `NULL`, then
 * the PRNG is initialised using the OS entropy, e.g. /dev/urandom.
 */
void alouette_random_set(
    unsigned long * seed);

#ifdef __cplusplus
}
#endif
#endif
