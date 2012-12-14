% LTFAT - Gabor analysis
%
%  Peter L. Søndergaard, 2007 - 2012.
%
%  Basic Time/Frequency analysis
%    TCONV          -  Twisted convolution
%    DSFT           -  Discrete Symplectic Fourier Transform
%    ZAK            -  Zak transform
%    IZAK           -  Inverse Zak transform
%    COL2DIAG       -  Move columns of a matrix to diagonals
%    S0NORM         -  Compute the S0-norm
%    TFMAT          -  Matrix of transform or operator in LTFAT
%
%  Gabor systems
%    dgt            -  Discrete Gabor transform
%    IDGT           -  Inverse discrete Gabor transform
%    ISGRAM         -  Iterative reconstruction from spectrogram
%    ISGRAMREAL     -  Iterative reconstruction from spectrogram (real signal)
%    DGT2           -  2D Discrete Gabor transform
%    IDGT2          -  2D Inverse discrete Gabor transform
%    DGTREAL        -  DGT for real-valued signals
%    IDGTREAL       -  |idgt|_ for real-valued signals
%    GABWIN         -  Evaluate Gabor window
%    DGTLENGTH      -  Length of Gabor system to expand a signal
%
%  Wilson bases and WMDCT
%    DWILT          -  Discrete Wilson transform
%    IDWILT         -  Inverse discrete Wilson transform
%    DWILT2         -  2-D Discrete Wilson transform
%    IDWILT2        -  2-D inverse discrete Wilson transform
%    WMDCT          -  Modified Discrete Cosine transform
%    IWMDCT         -  Inverse MDCT
%    WMDCT2         -  2-D MDCT
%    IWMDCT2        -  2-D inverse MDCT
%    WIL2RECT       -  Rectangular layout of Wilson coefficients
%    RECT2WIL       -  Inverse of WIL2RECT
%    WILWIN         -  Evaluate Wilson window
%
%  Reconstructing windows
%    GABDUAL        -  Canonical dual window
%    GABTIGHT       -  Canonical tight window
%    GABPROJDUAL    -  Dual window by projection
%    GABMIXDUAL     -  Dual window by mixing windows
%    WILORTH        -  Window of Wilson/WMDCT orthonormal basis
%    WILDUAL        -  Riesz dual window of Wilson/WMDCT basis 
%
%  Time/Frequency operators
%    GABMUL         -  Gabor multiplier
%    GABMULEIGS     -  Eigenpairs of Gabor multiplier
%    GABMULAPPR     -  Best approximation by a Gabor mult.
%    SPREADOP       -  Spreading operator
%    SPREADINV      -  Apply inverse spreading operator
%    SPREADADJ      -  Symbol of adjoint spreading operator
%    SPREADFUN      -  Symbol of operator expressed as a matrix
%    SPREADEIGS     -  Eigenpairs of spreading operator
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
%    GABREASSIGN    -  Reassign positive distribution
%
%  Phase conversions
%    PHASELOCK      -  Phase Lock Gabor coefficients to time
%    PHASEUNLOCK    -  Undo phase locking
%    SYMPHASE       -  Convert to symmetric phase
%
%  Support for non-separable lattices
%    MATRIX2LATTICETYPE - Matrix form to standard lattice description
%    LATTICETYPE2MATRIX - Standard lattice description to matrix form
%
%  Plots
%    TFPLOT         -  Plot coefficients on the time-frequency plane
%    PLOTDGT        -  Plot DGT coefficients
%    PLOTDGTREAL    -  Plot DGTREAL coefficients
%    PLOTDWILT      -  Plot DWILT coefficients
%    PLOTWMDCT      -  Plot WMDCT coefficients
%    SGRAM          -  Spectrogram based on DGT
%    GABIMAGEPARS   -  Choose paramets for nice Gabor image
%    RESGRAM        -  Reassigned spectrogram
%    INSTFREQPLOT   -  Plot of the instantaneous frequency
%    PHASEPLOT      -  Plot of STFT phase
%
%  For help, bug reports, suggestions etc. please send email to
%  ltfat-help@lists.sourceforge.net
