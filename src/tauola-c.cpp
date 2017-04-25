/*
 *  A minimalist C wrapper for tauola++
 *  Copyright (C) 2017  Valentin Niess
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRENTY; without even the implied warranty of
 *  MERCHENTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <float.h>

#include "Tauola/Tauola.h"
#include "Tauola/TauolaHEPEVTEvent.h"
#include "Tauola/TauolaHEPEVTParticle.h"

#include "tauola-c.h"
#include <fcntl.h>
#include <unistd.h>

using namespace std;
using namespace Tauolapp;

/* Static interface to HEPEVT */
static TauolaHEPEVTEvent event;

/* Particle cursor in the HEPEVT store. */
static int index = 0;

/* File identifiers for stdout, stderr and /dev/null. */
static int fd_out, fd_err, fd_null;
static fpos_t pos_out, pos_err;

/* Flag for muting the messages from TAUOLA++. */
static int mute = 0;

/* Restore the standard stdout and stderr. */
static void std_unmute(void)
{
        if (!mute) return;
        fflush(stdout);
        fflush(stderr);
        dup2(fd_out, fileno(stdout));
        dup2(fd_err, fileno(stderr));
        clearerr(stdout);
        clearerr(stderr);
        fsetpos(stdout, &pos_out);
        fsetpos(stderr, &pos_err);
}

/* Redirect stdout and stderr to /dev/null. */
static void std_mute(void)
{
        if (!mute) return;
        fflush(stdout);
        fflush(stderr);
        fgetpos(stdout, &pos_out);
        fgetpos(stderr, &pos_err);
        dup2(fd_null, fileno(stdout));
        dup2(fd_null, fileno(stderr));
}

/* Status flag for the wrapper's initialisation. */
static int initialised = 0;

/* Initialise TAUOLA and the wrapper. */
void tauola_initialise(int mute_, int * seed_p)
{
        if (initialised) return;

        /* Set the seed for the random engine. */
        int seed;
        if (seed_p == NULL) {
                seed = 0;
                FILE * stream = fopen("/dev/urandom", "rb");
                if (stream == NULL) goto seed_set;
                if (fread(&seed, sizeof(seed), 1, stream) <= 0) {
                        fclose(stream);
                        goto seed_set;
                }
                fclose(stream);
        } else {
                seed = *seed_p;
        }
seed_set:

        /* Set the redirection streams if muting. */
        if (mute_) {
                mute = 1;
                fd_null = open("/dev/null", O_RDWR | O_TRUNC);
                fd_out = dup(fileno(stdout));
                fd_err = dup(fileno(stderr));
        }

        /* Initialise TAUOLA library. */
        std_mute();
        Tauola::setSeed(seed, 0, 0);
        Tauola::initialize();
        std_unmute();
        initialised = 1;
}

/* Finalise the wrapper. */
void tauola_finalise(void)
{
        if (!initialised) return;
        if (mute) {
                close(fd_out);
                close(fd_err);
                mute = 0;
        }
        initialised = 0;
}

/* Uniform PRN over [0,1]. */
static double uniform01(void)
{
        return rand() / (double)(RAND_MAX);
}

/* Compute the rotation angle along the x or y *axis* in order to cancel the
 * y or x component while ensuring a positive z component.
 */
static double rotation_angle(int axis, const double p[3])
{
        const int coaxis = 1 - axis;
        if (p[2] == 0.) {
                if (p[coaxis] > 0.) return -0.5 * M_PI;
                else return 0.5 * M_PI;
        }
        else if (p[2] > 0.)
                return -atan(p[coaxis] / p[2]);
        else
                return M_PI - atan(p[coaxis] / p[2]);
}

/* Rotation along the X or Y axis. */
static void rotate(int axis, double phi, double v[3])
{
        const int coaxis = 1 - axis;
        const double vxy = v[coaxis];
        const double vz = v[2];
        const double c = cos(phi), s = sin(phi);
        v[coaxis] = c * vxy + s * vz;
        v[2] = -s * vxy + c * vz;
}

/* Compute the rotation angles w.r.t. TAUOLA's rest frame, i.e. in order to
 * align the tau's momentum with the z-axis.
 */
static void compute_angles(
    int pid, const double momentum[3], double * phi_x, double * phi_y)
{
        double p[3] = {momentum[0], momentum[1], momentum[2]};
        if (pid < 0) {
                *phi_x = rotation_angle(0, p);
                rotate(0, *phi_x, p);
                *phi_y = rotation_angle(1, p);
        }
        else {
                *phi_x = M_PI + rotation_angle(0, p);
                rotate(0, *phi_x, p);
                *phi_y = M_PI + rotation_angle(1, p);
        }
}

/* Decay a tau with TAUOLA. */
int tauola_decay(
    int pid, const double momentum[3], const double * polarisation)
{
        /* Build the event in the tau's rest frame. */
        event.clear();
        index = 0;
        if (abs(pid) != 15) return index;
        TauolaHEPEVTParticle * p = new TauolaHEPEVTParticle(
            pid, 1, 0., 0., 0., parmas_.amtau, parmas_.amtau, -1, -1, -1, -1);
        event.addParticle(p);

        /* Rotate the polarisation to TAUOLA's frame. */
        double pol[3];
        if (polarisation != NULL) {
                double phi_x, phi_y;
                pol[0] = polarisation[0];
                pol[1] = polarisation[1];
                pol[2] = polarisation[2];
                compute_angles(pid, momentum, &phi_x, &phi_y);
                rotate(1, -phi_y, pol);
                rotate(0, -phi_x, pol);
        }

        /* Decay the tau. */
        std_mute();
        for (;;) {
                p->decay();
                if (polarisation == NULL) break;
                const double w = 0.5 * (1. +
                    p->getPolarimetricX() * pol[0] +
                    p->getPolarimetricY() * pol[1] +
                    p->getPolarimetricZ() * pol[2]);
                if (uniform01() <= w) break;
        }
        p->addDecayToEventRecord();
        p->decayEndgame();
        std_unmute();

        /* Check the result. */
        p = event.getParticle(1);
        if (isinf(p->getPx())) {
                event.clear();
                index = 0;
                return index;
        }

        /* Boost the daughters to the lab frame. */
        const double tau[3] = {momentum[0] / parmas_.amtau,
            momentum[1] / parmas_.amtau, momentum[2] / parmas_.amtau};
        const double gamma =
            sqrt(1. + tau[0] * tau[0] + tau[1] * tau[1] + tau[2] * tau[2]);
        index = 1;
        int i;
        for (i = 0; i < event.getParticleCount(); i++) {
                p = event.getParticle(i);
                if (p->getStatus() != 1) {
                        continue;
                }
                const int aid = abs(p->getPdgID());
                if ((aid == 24) || (aid > 9999)) {
                        /* Correct the status of temporaries. */
                        p->setStatus(2);
                        continue;
                }
                if (gamma <= 1. + FLT_EPSILON) continue;

                /* Apply the boost. */
                const double P[3] = {p->getPx(), p->getPy(), p->getPz()};
                const double E = p->getE();
                const double ptau = P[0] * tau[0] + P[1] * tau[1]
                    + P[2] * tau[2];
                const double tmp = ptau / (gamma + 1.) + E;
                p->setE(gamma * E + ptau);
                p->setPx(P[0] + tmp * tau[0]);
                p->setPy(P[1] + tmp * tau[1]);
                p->setPz(P[2] + tmp * tau[2]);
        }

        return index;
}

/* Backward decay from a tau neutrino to a tau. */
int tauola_undecay(int pid, const double momentum[3],
    const double * polarisation, double * weight)
{
        *weight = 0.;
        event.clear();
        index = 0;
        if (abs(pid) != 16) return index;
        int tau_pid = pid > 0 ? 15 : -15;
        TauolaHEPEVTParticle * p = new TauolaHEPEVTParticle(tau_pid, 1, 0.,
            0., 0., parmas_.amtau, parmas_.amtau, -1, -1, -1, -1);
        event.addParticle(p);
        std_mute();
        if (polarisation != NULL) {
                Tauola::decayOne(p, true, polarisation[0], polarisation[1],
                    polarisation[2]);
        } else {
                event.decayTaus();
        }
        std_unmute();

        /* Check the result. */
        TauolaHEPEVTParticle * pi = event.getParticle(1);
        if (isinf(pi->getPx())) {
                event.clear();
                index = 0;
                return index;
        }

        /* Compute the parameters of the frame transform. */
        const double energy = sqrt(momentum[0] * momentum[0] +
            momentum[1] * momentum[1] + momentum[2] * momentum[2]);
        int i;
        for (i = 1; i < event.getParticleCount(); i++) {
                pi = event.getParticle(i);
                if (pi->getPdgID() == pid) break;
        }
        const double momentum0[3] = {pi->getPx(), pi->getPy(), pi->getPz()};
        const double energy0 = sqrt(momentum0[0] * momentum0[0] +
            momentum0[1] * momentum0[1] + momentum0[2] * momentum0[2]);

        const double d = energy * energy0 + momentum[0] * momentum0[0] +
            momentum[1] * momentum0[1] + momentum[2] * momentum0[2];
        const double ee = energy + energy0;
        const double gamma = ee * ee / d - 1.;
        const double t0 = ee / d;
        const double tau[3] = {t0 * (momentum[0] - momentum0[0]),
            t0 * (momentum[1] - momentum0[1]),
            t0 * (momentum[2] - momentum0[2])};

        /* Set the backward Monte-Carlo weight. */
        *weight = parmas_.amtau * parmas_.amtau * parmas_.amtau *
            fabs(t0 * t0 * gamma / energy);

        /* Apply the boost to the decay products */
        double Et = 0., Pt[3] = {0., 0., 0.};
        int j;
        for (j = 1; j < event.getParticleCount(); j++) {
                if (j == i) {
                        Et += energy;
                        Pt[0] += momentum[0];
                        Pt[1] += momentum[1];
                        Pt[2] += momentum[2];
                        continue;
                }

                TauolaHEPEVTParticle * pj = event.getParticle(j);
                if (pj->getStatus() != 1)
                        continue;

                const double P[3] = {pj->getPx(), pj->getPy(), pj->getPz()};
                const double E = pj->getE();
                const double ptau = P[0] * tau[0] + P[1] * tau[1]
                    + P[2] * tau[2];
                const double tmp = ptau / (gamma + 1.) + E;
                pj->setE(gamma * E + ptau);
                pj->setPx(P[0] + tmp * tau[0]);
                pj->setPy(P[1] + tmp * tau[1]);
                pj->setPz(P[2] + tmp * tau[2]);
                Et += pj->getE();
                Pt[0] += pj->getPx();
                Pt[1] += pj->getPy();
                Pt[2] += pj->getPz();
        }

        /* Replace the daughter particle with the tau mother. In order to ensure
         * the conservation of the energy-momentum, the mother's 4 momentum has
         * been computed from the boosted products.
         */
        pi->setPdgID(tau_pid);
        pi->setPx(Pt[0]);
        pi->setPy(Pt[1]);
        pi->setPz(Pt[2]);
        pi->setE(Et);
        pi->setMass(parmas_.amtau);

        index = 1;
        return index;
}

/* Iterator over the tau decay products. */
int tauola_product(int * pid, double momentum[3])
{
        if (!index) return 0;
        while (index < event.getParticleCount()) {
                TauolaHEPEVTParticle * p = event.getParticle(index);
                if (p->getStatus() != 1) {
                        index++;
                        continue;
                }
                *pid = p->getPdgID();
                momentum[0] = p->getPx();
                momentum[1] = p->getPy();
                momentum[2] = p->getPz();
                return index++;
        }
        event.clear();
        index = 0;
        return 0;
}
