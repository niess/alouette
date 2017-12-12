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
        ALOUETTE_RETURN_DOMAIN_ERROR,
        /** A floating point error occured. */
        ALOUETTE_RETURN_FLOATING_ERROR,
        /** A read /write error occured. */
        ALOUETTE_RETURN_IO_ERROR,
        /** A file couldn't be found. */
        ALOUETTE_RETURN_PATH_ERROR,
        /** A TAUOLA error occured. */
        ALOUETTE_RETURN_TAULOA_ERROR,
        /** The number of ALOUETTE return codes.  */
        ALOUETTE_N_RETURNS
};

/**
 * Initialise the wrapper and TAUOLA.
 *
 * @param mute    Flag to mute all low level messages from TAUOLA.
 * @param state   The pseudo random state vector for TAUOLA's builtin PRNG or
 *                  `NULL`.
 * @param offset  The offset to apply to the random sequence.
 * @return On success `ALOUETTE_RETURN_SUCCESS` is returned otherwise an error
 * code is returned as detailed below.
 *
 * Initialise the TAUOLA++ library. Call this function prior to any other
 * wrapper's routine. If *mute* is not zero all messages from TAUOLA are
 * redirected to /dev/null.
 *
 * __Note__ : if *state* is `ǸULL` the pseudo random engine is initialised from
 * /dev/urandom.
 *
 * __Error codes__
 *
 *     ALOUETTE_RETURN_IO_ERROR        Couldn't read from /dev/urandom.
 *
 *     ALOUETE_RETURN_PATH_ERROR       Couldn't open /dev/urandom.
 *
 *     ALOUETTE_RETURN_TAUOLA_ERROR    A TAUOLA error occured.
 */
enum alouette_return alouette_initialise(int mute, unsigned int state[3]);

/**
 * Finalise the wrapper.
 *
 * @return On success `ALOUETTE_RETURN_SUCCESS` is returned otherwise an error
 * code is returned as detailed below.
 *
 * __Error codes__
 *
 *     ALOUETTE_RETURN_TAUOLA_ERROR    A TAUOLA error occured.
 */
enum alouette_return alouette_finalise(void);

/**
 * Return a string describing an `alouette_return` code.
 *
 * This function is analog to the C89 `strerror` function but specific to
 * ALOUETTE return codes.
 */
const char * alouette_strerror(enum alouette_return rc);

/**
 * Perform a forward Monte-Carlo tau decay.
 *
 * @param pid             The PDG ID of the decaying tau, i.e. 15 or -15.
 * @param momentum        The tau momentum at decay, in GeV/c.
 * @param polarisation    The tau polarisation vector, or `NULL`.
 * @return On success `ALOUETTE_RETURN_SUCCESS` is returned otherwise an error
 * code is returned as detailed below.
 *
 * Simulate a tau decay with TAUOLA. An optionnal polarisation 3-vector can
 * be provided. If `ǸULL` spin effects are neglected.
 *
 * __Error codes__
 *
 *     ALOUETTE_RETURN_DOMAIN_ERROR      The provided *pid* is not valid.
 *
 *     ALOUETTE_RETURN_FLOATING_ERROR    A floating point error occured.
 *
 *     ALOUETTE_RETURN_TAUOLA_ERROR      A TAUOLA error occured.
 */
enum alouette_return alouette_decay(
    int pid, const double momentum[3], const double * polarisation);

/**
 * Callback for the tau polarisation in backward decays.
 *
 * @param pid             The PDG ID of the tau mother, i.e. 15 or -15.
 * @param momentum        The tau's momentum at decay, in GeV/c.
 * @param polarisation    The tau polarisation.
 *
 * In a backward decay, the spin polarisation of the tau mother is not known
 * a priori. The user can supply an a posteriori value with this callback.
 */
typedef void alouette_polarisation_cb(
    int pid, const double momentum[3], double * polarisation);

/**
 * Perform a backward Monte-Carlo tau decay.
 *
 * @param pid             The PDG ID of the tau neutrino decaying product,
 *                        i.e. 16 or -16.
 * @param momentum        The neutrino momentum after decay, in GeV/c.
 * @param polarisation    A callback for the primary tau spin polarisation or
 *                        `NULL`.
 * @param bias            Tuning parameter for the bias.
 * @param weight          The backward Monte-Carlo weight.
 * @return On success `ALOUETTE_RETURN_SUCCESS` is returned otherwise an error
 * code is returned as detailed below.
 *
 * Simulate a backward tau decay from a tau neutrino product. The polarisation
 * of the primary tau can be provided a posteriori. Set *polarisation* to
 * `NULL` in order to ignore spin effects.
 *
 * The *bias* parameters allows to control the biasing of the angular
 * distribution of decay products in the mother's rest frame. It must be in the
 * range ]-1, +inf[, though values close to `-1` or much higher than `10` are
 * likely to be numerically unstable. Set it to to `0` for an unbiased
 * distribution, i.e. isotropic.
 *
 * __Error codes__
 *
 *     ALOUETTE_RETURN_DOMAIN_ERROR      The provided *pid* is not valid.
 *
 *     ALOUETTE_RETURN_FLOATING_ERROR    A floating point error occured.
 *
 *     ALOUETTE_RETURN_TAUOLA_ERROR      A TAUOLA error occured.
 */
enum alouette_return alouette_undecay(int pid, const double momentum[3],
    alouette_polarisation_cb * polarisation, double bias, double * weight);

/**
 * Iterator over the decay products.
 *
 * @param pid         The PDG ID of the retrieved product.
 * @param momentum    The momentum of the retrieved product, in GeV/c.
 * @return On success `ALOUETTE_RETURN_SUCCESS` is returned otherwise an error
 * code is returned as detailed below.
 *
 * Loop over this routine after an `alouette_decay` in order to retrieve all
 * the decay products.
 *
 * __Warning__ : the decay products are consumed by the iterator.
 *
 * __Error codes__
 *
 *     ALOUETTE_RETURN_DOMAIN_ERROR    No product is available.
 */
enum alouette_return alouette_product(int * pid, double momentum[3]);

/**
 * Getter for the polarimetric vector of the last decay.
 *
 * @param polarimetric    The polarimetric vector to fill.
 * @return On success `ALOUETTE_RETURN_SUCCESS` is returned otherwise an error
 * code is returned as detailed below.
 *
 * Call this routine after an `alouette_decay` or `alouette_undecay` in order
 * to retrieve the polarimetric vector of the decay.
 *
 * __Error codes__
 *
 *     ALOUETTE_RETURN_DOMAIN_ERROR    No decay is available.
 */
enum alouette_return alouette_polarimetric(double polarimetric[3]);

/**
 * Get TAUOLA's builtin random state.
 *
 * @param state    The state vector to fill.
 */
void alouette_random_state(unsigned int state[3]);

#ifdef __cplusplus
}
#endif
#endif
