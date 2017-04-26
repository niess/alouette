/*
 * Copyright (C) 2017 CNRS/IN2P3
 * Author: Valentin NIESS (niess@in2p3.fr)
 *
 * This software is a C library whose purpose is to transport high energy
 * muons or taus in various media.
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

/**
 * Initialise the wrapper and TAUOLA++.
 *
 * @param mute    Flag to mute all low level messages from TAUOLA.
 * @param seed    The seed for TAUOLA's pseudo random engine or `NULL`.
 *
 * Initialise the TAUOLA++ library. Call this function prior to any other
 * wrapper's routine. If *mute* is not zero all messages from TAUOLA are
 * redirected to /dev/null.
 *
 * __Note__ : if *seed* is `ǸULL` the pseudo random engine is initialised from
 * /dev/urandom.
 */
void alouette_initialise(int mute, int * seed);

/**
 * Finalise the wrapper.
 */
void alouette_finalise(void);

/**
 * Perform a forward Monte-Carlo tau decay.
 *
 * @param pid             The PDG ID of the decaying tau, i.e. 15 or -15.
 * @param momentum        The tau momentum at decay, in GeV/c.
 * @param polarisation    The tau polarisation vector, or `NULL`.
 * @return                `0` if the decay failed. A non null integer
 *                        otherwise.
 *
 * Simulate a tau decay with TAUOLA. An optionnal polarisation 3-vector can
 * be provided. If `ǸULL` spin effects are neglected.
 */
int alouette_decay(
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
typedef void polarisation_cb(
    int pid, const double momentum[3], double * polarisation);

/**
 * Perform a backward Monte-Carlo tau decay.
 *
 * @param pid             The PDG ID of the tau neutrino decaying product,
 *                        i.e. 16 or -16.
 * @param momentum        The neutrino momentum after decay, in GeV/c.
 * @param polarisation    A callback for the primary tau spin polarisation or
 *                        `NULL`.
 * @param weight          The backward Monte-Carlo weight.
 * @return                `0` if the backward decay failed. A non null integer
 *                        otherwise.
 *
 * Simulate a backward tau decay from a tau neutrino product. The polarisation
 * of the primary tau can be provided a posteriori. Set *polarisation* to
 * `NULL` in order to ignore spin effects.
 */
int alouette_undecay(int pid, const double momentum[3],
   polarisation_cb * polarisation, double * weight);

/**
 * Iterator over the decay products.
 *
 * @param pid         The PDG ID of the retrieved product.
 * @param momentum    The momentum of the retrieved product, in GeV/c.
 * @return            `0` if no more product is available. A non null integer
 *                    otherwise.
 *
 * Loop over this routine after an `alouette_decay` in order to retrieve all
 * the decay products.
 *
 * __Warning__ : the decay products are consumed by the iterator.
 */
int alouette_product(int * pid, double momentum[3]);

#ifdef __cplusplus
}
#endif
#endif
