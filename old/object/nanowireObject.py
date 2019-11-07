#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 13 14:01:03 2019

@author: Domi
"""

import numpy as np
import tinyarray as ta

class SimpleNamespace(object):
    """A simple container for parameters."""
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)

s_0 = np.identity(2)
s_z = np.array([[1., 0.], [0., -1.]])
s_x = np.array([[0., 1.], [1., 0.]])
s_y = np.array([[0., -1j], [1j, 0.]])

tau_z = ta.array(np.kron(s_z, s_0))
tau_x = ta.array(np.kron(s_x, s_0))
tau_y = ta.array(np.kron(s_y, s_0))
sig_z = ta.array(np.kron(s_0, s_z))
sig_x = ta.array(np.kron(s_0, s_x))
sig_y = ta.array(np.kron(s_0, s_y))
tau_zsig_x = ta.array(np.kron(s_z, s_x))

def main():
    ...

if __name__ == '__main__':
    main()