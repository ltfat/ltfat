# The Large Time Frequency Analysis Toolbox (LTFAT).

This is the official development repository of the Large Time Frequency Analysis Toolbox.
LTFAT is an open source toolbox for audio signal processing and harmonic analysis.
It is written in GNU Octave and C/C++. The development repository of the C/C++ backend library
can be found [here](https://github.com/ltfat/libltfat).

The phase retrieval toolbox, a collection of phase-reconstruction algorithms for time-frequency 
representations can be retrieved from [here](https://github.com/ltfat/phaseret).

# Usage
Download the release applicable to your operating system [here](https://github.com/ltfat/ltfat/releases).
To use the toolbox, in Matlab/Octave, 'cd' to the LTFAT directory and execute 'ltfatstart'.
In Octave you can put this command in your ~/.octaverc file.

# Directory structure.

The toolbox is organized in subdirectories including the following:

  fourier     - Fourier analysis, DCT/DST transforms and filtering.
  gabor       - Gabor, Wilson and WMDCT analysis/synthesis functions.
  wavelets    - Wavelet analysis/synthesis, wavelet filter banks
  filterbank  - General and auditory inspired filterbanks. 
  nonstatgab  - Non-stationary Gabor frames.
  quadratic   - Quadratic distributions
  frames      - Frame objects.
  operators   - Linear operators associated with frames.
  sigproc     - A collection of simple, signal processing tools.
  blockproc   - Methods for block processing and block-adapted transforms.
  auditory    - Auditory scales and common filter types.
  demos       - Demos showing applications or aspects of LTFAT functions.
  signals     - Test signals for use with the examples.
  comp        - Computational subroutines. These should not be called directly.
  src         - C implementation of the most computationally intensive routines.
  mex         - Mex files to speed up the toolbox (if compiled)
  oct         - Octave C++-files to speed up the toolbox (if compiled)

The file INSTALL contains instructions for compiling the C-library and 
the Octave and Matlab interfaces.
