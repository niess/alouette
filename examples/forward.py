#! /usr/bin/env python3
import alouette


for i in range(3):
    product = alouette.decay(momentum=(0, 0, 1), polarisation=(0, 0, -1))

    print(f'# Event {i + 1}')
    print(product)
