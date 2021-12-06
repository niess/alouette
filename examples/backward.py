#! /usr/bin/env python3
import alouette
import numpy


polarisation = -1


def polarisation_cb(pid, momentum):
    '''Callback for setting the polarisation in backward decays
    '''
    return polarisation * momentum / numpy.linalg.norm(momentum)


for i in range(3):
    product = alouette.undecay(momentum=(0, 0, 1),
                               polarisation=polarisation_cb,
                               bias=polarisation)

    print(f'# Event {i + 1} ({product.weight:.5E})')
    print(product)
