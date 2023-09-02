% LTFAT - Gabor analysis
%
%  Peter L. Søndergaard, 2007 - 2023.
%
%  Basic Time/Frequency analysis
%    TCONV          -  Twisted convolution
%    DSFT           -  Discrete Symplectic Fourier Transform
%    ZAK            -  Zak transform
%    IZAK           -  Inverse Zak transform
%    COL2DIAG       -  Move columns of a matrix to diagonals
%    S0NORM         -  S0-norm
%
%  Gabor systems
%    DGT            -  Discrete Gabor transform
%    IDGT           -  Inverse discrete Gabor transform
%    ISGRAM         -  Iterative reconstruction from spectrogram
%    ISGRAMREAL     -  Iterative reconstruction from spectrogram for real-valued signals
%    DGT2           -  2D Discrete Gabor transform
%    IDGT2          -  2D Inverse discrete Gabor transform
%    DGTREAL        -  |dgt| for real-valued signals
%    IDGTREAL       -  |idgt| for real-valued signals
%    GABWIN         -  Evaluate Gabor window
%    PROJKERN	    -  Projection of Gabor coefficients onto kernel space
%    DGTLENGTH      -  Length of Gabor system to expand a signal
%
%  Multi-Gabor systems
%    MULTIDGTREALMP -  Matching pursuit decomposition in Multi-Gabor system 
%
%  Wilson bases and WMDCT
%    DWILT          -  Discrete Wilson transform
%    IDWILT         -  Inverse discrete Wilson transform
%    DWILT2         -  2-D Discrete Wilson transform
%    IDWILT2        -  2-D inverse discrete Wilson transform
%    WMDCT          -  Modified Discrete Cosine transform
%    IWMDCT         -  Inverse |wmdct|
%    WMDCT2         -  2-D |wmdct|
%    IWMDCT2        -  2-D inverse |wmdct|
%    WIL2RECT       -  Rectangular layout of Wilson coefficients
%    RECT2WIL       -  Inverse of |wil2rect|
%    WILWIN         -  Evaluate Wilson window
%    DWILTLENGTH    -  Length of Wilson/WMDCT system to expand a signal
%
%  Reconstructing windows
%    GABDUAL        -  Canonical dual window
%    GABTIGHT       -  Canonical tight window
%    GABFIRDUAL     -  FIR optimized dual window 
%    GABOPTDUAL     -  Optimized dual window
%    GABFIRTIGHT    -  FIR optimized tight window
%    GABOPTTIGHT    -  Optimized tight window
%    GABCONVEXOPT   -  Optimized window
%    GABPROJDUAL    -  Dual window by projection
%    GABMIXDUAL     -  Dual window by mixing windows
%    WILORTH        -  Window of Wilson/WMDCT orthonormal basis
%    WILDUAL        -  Riesz dual window of Wilson/WMDCT basis 
%
%  Conditions numbers
%    GABFRAMEBOUNDS -  Frame bounds of Gabor system
%    GABRIESZBOUNDS -  Riesz sequence/basis bounds of Gabor system
%    WILBOUNDS      -  Frame bounds of Wilson basis
%    GABDUALNORM    -  Test if two windows are dual
%    GABFRAMEDIAG   -  Diagonal of Gabor frame operator
%    WILFRAMEDIAG   -  Diagonal of Wilson/WMDCT frame operator
%
%  Phase gradient methods and reassignment
%    GABPHASEGRAD   -  Instantaneous time/frequency from signal
%    GABPHASEDERIV  -  Phase derivatives
%    GABPHASEDERIVREAL  -  Phase derivatives for |dgtreal|
%    GABREASSIGN    -  Reassign positive distribution
%    GABREASSIGNREAL    -  Reassign positive distribution for real signals
%    GABREASSIGNADJUST - Adjustable t-f reassignment
%
%  Phase reconstruction
%    CONSTRUCTPHASE     - Phase construction from abs. values of |dgt|
%    CONSTRUCTPHASEREAL - |constructphase| for |dgtreal|
%
%  Phase conversions
%    PHASELOCK       -  Phase Lock Gabor coefficients to time. inv.
%    PHASEUNLOCK     -  Undo phase locking
%    PHASELOCKREAL   -  Same as |phaselock| for |dgtreal|
%    PHASEUNLOCKREAL -  Same as |phaseunlock| for |idgtreal|
%    SYMPHASE        -  Convert to symmetric phase
%
%  Support for non-separable lattices
%    MATRIX2LATTICETYPE - Matrix form to standard lattice description
%    LATTICETYPE2MATRIX - Standard lattice description to matrix form
%    SHEARFIND      -  Shears to transform a general lattice to a separable
%    NOSHEARLENGTH  -  Next transform side not requiring a frequency side shear
%
%  Plots
%    TFPLOT         -  Plot coefficients on the time-frequency plane
%    PLOTDGT        -  Plot |dgt| coefficients
%    PLOTDGTREAL    -  Plot |dgtreal| coefficients
%    PLOTDWILT      -  Plot |dwilt| coefficients
%    PLOTWMDCT      -  Plot |wmdct| coefficients
%    SGRAM          -  Spectrogram based on |dgt|
%    GABIMAGEPARS   -  Choose parameters for nice Gabor image
%    RESGRAM        -  Reassigned spectrogram
%    INSTFREQPLOT   -  Plot of the instantaneous frequency
%    PHASEPLOT      -  Plot of STFT phase
%
%  For help, bug reports, suggestions etc. please visit 
%  http://github.com/ltfat/ltfat/issues
