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

/* Decay a tau with TAUOLA. */
void tauola_decay(int pid, double momentum[3], double * polarisation)
{
        event.clear();
        if (abs(pid) != 15) return;
        TauolaHEPEVTParticle * tau = new TauolaHEPEVTParticle(
            pid, 1, 0., 0., 0., parmas_.amtau, parmas_.amtau, -1, -1, -1, -1);
        event.addParticle(tau);
        tau->setPx(momentum[0]);
        tau->setPy(momentum[1]);
        tau->setPz(momentum[2]);
        const double p2 = momentum[0] * momentum[0] +
            momentum[1] * momentum[1] + momentum[2] * momentum[2];
        tau->setE(sqrt(p2 + parmas_.amtau * parmas_.amtau));
        std_mute();
        if (polarisation != NULL) {
                Tauola::decayOne(tau, false, polarisation[0], polarisation[1],
                    polarisation[2]);
        } else {
                event.decayTaus();
        }
        std_unmute();
        index = 1;
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
