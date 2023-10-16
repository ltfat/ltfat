% LTFAT - Signal processing tools
%
%  Peter L. Søndergaard, 2007 - 2023.
%
%  General
%    RMS            -  Root Mean Square norm of signal.
%    SETNORM        -  Normalize signal by specified norm.
%    GAINDB         -  Scale input signal.
%    CRESTFACTOR    -  Compute the crest factor of a signal.
%    LOWDISCREPANCY -  Compute a low discrepancy sequence.
%    UQUANT         -  Simulate uniform quantization.
%    POSTPAD        -  Pad or truncate a vector.
%
%  Window functions
%    FIRWIN         -  FIR windows (Hanning,Hamming,Blackman,...).
%    FIRKAISER      -  FIR Kaiser-Bessel window.
%    FIR2LONG       -  Extend FIR window to LONG window.
%    LONG2FIR       -  Cut LONG window to FIR window.
%    FREQWIN        -  Freq responses (Gauss,Gammatone,Butterworth,...)
%
%  Wavelet functions
%    FREQWAVELET    -  Frequency responses of wavelets (Cauchy, Morse)
%
%  Filtering
%    FIRFILTER      -  Construct an FIR filter.
%    BLFILTER       -  Construct a band-limited filter.
%    WARPEDBLFILTER -  Warped, band-limited filter.
%    FREQFILTER     -  Construct a full length frequency side filter.
%    PFILT          -  Apply filter with periodic boundary conditions.
%    MAGRESP        -  Magnitude response plot.
%    TRANSFERFUNCTION - Compute the transfer function of a filter.
%    PGRPDELAY      -  Periodic Group Delay
%
%  Ramping
%    RAMPUP         -  Rising ramp.
%    RAMPDOWN       -  Falling ramp.
%    RAMPSIGNAL     -  Ramp a signal.
%
%  Thresholding methods
%    THRESH         -  Coefficient thresholding.
%    LARGESTR       -  Keep largest ratio of coefficients.
%    LARGESTN       -  Keep N largest coefficients.
%    DYNLIMIT       -  Limit the dynamical range.
%    GROUPTHRESH    -  Group thresholding.
%
%  Image processing
%    RGB2JPEG       -  Convert RGB values to the JPEG colour model
%    JPEG2RGB       -  Convert values from the JPEG colour model to RGB
%
%  Tools for OFDM
%    QAM4           -  Quadrature amplitude modulation, order 4
%    IQAM4          -  Inverse QAM of order 4
%
%  Tonal-transient separation
%    tfjigsawsep      - Tonal-transient-residual separation using the T-F jigsaw puzzle algorithm.
%    plottfjigsawsep  - Plot the separated layers.
%
%  For help, bug reports, suggestions etc. please visit 
%  http://github.com/ltfat/ltfat/issues
%
