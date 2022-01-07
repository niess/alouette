#! /usr/bin/env python3
import alouette
import numpy


def polarisation(pid, momentum):
    '''Callback for setting the polarisation in backward decays

       Returns a 100% longitudinal polarisation, i.e. a left handed tau- or a
       right handed tau+.
    '''
    polar = -1. if pid > 0 else 1.
    return polar * momentum / numpy.linalg.norm(momentum)


for i in range(3):
    product = alouette.undecay(momentum=(0, 0, 1), polarisation=polarisation)

    print(f'# Event {i + 1} ({product.weight:.5E})')
    print(product)
