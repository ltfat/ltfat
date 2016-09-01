#!/usr/bin/python3
# encoding: utf-8

from testhelper import lib, ffi
import numpy as np


def test_fftshift():
    inArr = np.array([1, 2, 3, 4, 5, 6, 7], dtype=np.float64)
    res = np.array([5, 6, 7, 1, 2, 3, 4],  dtype=np.float64)
    status = getattr(lib, 'ltfat_fftshift_d')(
        ffi.from_buffer(inArr), 7, ffi.from_buffer(inArr))
    assert all(inArr == res)
