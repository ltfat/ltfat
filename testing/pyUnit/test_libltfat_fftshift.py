#!/usr/bin/python3
# encoding: utf-8


from cffi import FFI
import os
import subprocess
import numpy as np
import matplotlib.pyplot as plt


def main():
    ffi = FFI()
    filedir = os.path.dirname(os.path.abspath(__file__))
    libpath = os.path.join(filedir, '..', '..', 'build')
    libltfat = os.path.join(libpath, 'libltfat.so')
    header = os.path.join(libpath, 'ltfat_flat.h')

    with open(header) as f_header:
        ffi.cdef(f_header.read())

    lib = ffi.dlopen(libltfat)

    inArr = np.array(range(1,11),dtype=np.float64)

    lib.ltfat_fftshift_d(ffi.from_buffer(inArr),10,ffi.from_buffer(inArr))
    plt.stem(inArr.real())
    plt.show()


if __name__ == '__main__':
    main()
