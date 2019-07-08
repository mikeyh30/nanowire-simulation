#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul  6 21:34:22 2019

@author: Domi
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul  6 15:21:30 2019

@author: Domi
"""

# Tutorial 2.8.2. 3D example: zincblende structure
# ================================================
#
# Physical background
# -------------------
#  - 3D Bravais lattices
#
# Kwant features highlighted
# --------------------------
#  - demonstrate different ways of plotting in 3D

import kwant

lat = kwant.lattice.general([(0, 0, 1), (0, 1, 0), (1, 0, 0)])
#a, b = lat.sublattices

def make_cuboid(a=1, b=7):
    def cuboid_shape(pos):
        x, y, z = pos
        return 0 <= x < a and 0 <= y**2 + z**2 < b

    syst = kwant.Builder()
    syst[lat.shape(cuboid_shape, (0, 0, 0))] = None
    syst[lat.neighbors()] = None

    return syst


def main():
    # the standard plotting style for 3D is mainly useful for
    # checking shapes:
    syst = make_cuboid()

    kwant.plot(syst, site_size=0.18, site_lw=0.01, hop_lw=0.05)

    # visualize the crystal structure better for a very small system
#    syst = make_cuboid(a=1.5, b=1.5)
#
#    def family_colors(site):
#        return 'r' if site.family == a else 'g'
#
#    kwant.plot(syst, site_size=0.18, site_lw=0.01, hop_lw=0.05,
#               site_color=family_colors)


# Call the main function if the script gets executed (as opposed to imported).
# See <http://docs.python.org/library/__main__.html>.
if __name__ == '__main__':
    main()