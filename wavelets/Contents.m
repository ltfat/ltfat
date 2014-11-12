% LTFAT - Wavelets
%
%   Zdenek Prusa, 2013 - 2014.
%
%   Basic analysis/synthesis
%      FWT               - Fast Wavelet Transform 
%      IFWT              - Inverse Fast Wavelet Transform
%      FWT2              - 2D Fast Wavelet Transform 
%      IFWT2             - 2D Inverse Fast Wavelet Transform
%      UFWT              - Undecimated Fast Wavelet Transform
%      IUFWT             - Inverse Undecimated Fast Wavelet Transform 
%      FWTLENGTH         - Length of Wavelet system to expand a signal
%      FWTCLENGTH        - Lengths of the wavelet coefficient subbands
%
%   Advanced analysis/synthesis
%      WFBT              - Transform using general Wavelet Filterbank Tree 
%      IWFBT             - Inverse transform using general Wavelet Filterbank Tree
%      UWFBT             - Undecimated transform using general Wavelet Filterbank Tree 
%      IUWFBT            - Inverse Undecimated transform using general Wavelet Filterbank Tree
%      WPFBT             - Wavelet Packet Transform using general Wavelet Filterbank Tree 
%      IWPFBT            - Inverse Wavelet Packet Transform using general Wavelet Filterbank Tree
%      UWPFBT            - Undecimated Wavelet Packet Transform using general Wavelet Filterbank Tree 
%      IUWPFBT           - Inverse Undecimated Wavelet Packet Transform using general Wavelet Filterbank Tree
%      WPBEST            - Best Tree selection
%      WFBTLENGTH        - Length of Wavelet filterbank system to expand a signal
%      WFBTCLENGTH       - Lengths of Wavelet filterbank coefficent subbands
%      WPFBTCLENGTH      - Lengths of Wavelet Packet transform coefficent subbands
%
%   Dual-tree complex wavelet transform
%      DTWFB             - Dual-Tree Wavelet FilterBank
%      IDTWFB            - Inverse Dual-Tree Wavelet FilterBank
%      DTWFBREAL         - Dual-Tree Wavelet FilterBank for real signals
%      IDTWFBREAL        - Inverse Dual-Tree Wavelet FilterBank for real signals
%
%   Wavelet Filterbank trees manipulation
%      WFBTINIT          - Wavelet Filterbank tree structure initialization
%      DTWFBINIT         - Dual-Tree wavelet filterbank structure initialization
%      WFBTPUT           - Puts node (basic filterbank) to the specific  tree coordinates
%      WFBTREMOVE        - Removes node (basic filterbank) from the specific tree coordinates
%      WFBT2FILTERBANK   - WFBT or FWT non-iterated filterbank using the multirate identity
%      WPFBT2FILTERBANK  - WPFBT non-iterated filterbank using the multirate identity
%      DTWFB2FILTERBANK  - DTWFB or DTWFBREAL non-iterated filterbank
%      FWTINIT           - Basic Wavelet Filters structure initialization
%
%   Frame properties of wavelet filterbanks:
%      WFBTBOUNDS        - Frame bounds of WFBT or FWT
%      WPFBTBOUNDS       - Frame bounds of WPFBT
%      UWFBTBOUNDS       - Frame bounds of UWFBT or UFWT
%      UWPFBTBOUNDS      - Frame bounds of UWPFBT
%      DTWFBBOUNDS       - Frame bounds of DTWFB
%  
%   Plots
%      PLOTWAVELETS      - Plot wavelet coefficients
%      WFILTINFO         - Plot wavelet filters impulse and frequency responses and approximation of scaling and wavelet functions
%      WFILTDTINFO       - Plot the same as WFILTINFO but for dual-tree wavelet transform
%
%   Auxilary
%      WAVFUN            - Aproximate of the continuous scaling and wavelet functions
%      WAVCELL2PACK      - Changes wavelet coefficient storing format
%      WAVPACK2CELL      - Changes wavelet coefficient storing format back
%
%   Wavelet Filters defined in the time-domain
%      WFILT_ALGMBAND    - An ALGebraic construction of orthonormal M-BAND wavelets with perfect reconstruction
%      WFILT_CMBAND      - M-Band cosine modulated wavelet filters
%      WFILT_COIF        - Coiflets
%      WFILT_DB          - DauBechies orthogonal filters (ortonormal base)
%      WFILT_DDEN        - Double-DENsity dwt filters (tight frame)
%      WFILT_DGRID       - Dense GRID framelets (tight frame, symmetric)
%      WFILT_HDEN        - Higher DENsity dwt filters (tight frame, frame)  
%      WFILT_LEMARIE       - Battle and Lemarie quadrature filters
%      WFILT_MATLABWRAPPER - Wrapper of the wfilters function from the Matlab Wavelet Toolbox 
%      WFILT_MBAND           - M-band filters
%      WFILT_REMEZ           - Wavelet orthonogal filters based on the Remez Exchange algorithm
%      WFILT_SYMDS           - SYMmetric wavelet Dyadic Siblings (frames)
%      WFILT_SPLINE          - Biorthogonal spline wavelet filters
%      WFILT_SYM             - Least asymmetric Daubechies wavelet filters
%      WFILT_SYMDDEN         - Symmetric Double-DENsity dwt filters (tight frame)
%      WFILT_SYMORTH         - Symmetric nearly-orthogonal and orthogonal nearly-symmetric wav. filters
%      WFILT_SYMTIGHT        - Symmetric nearly shift-invariant tight frame wavelets
%      WFILT_QSHIFTA         - First tree filters from WFILTDT_QSHIFT 
%      WFILT_QSHIFTB         - Second tree filters from WFILTDT_QSHIFT 
%      WFILT_ODDEVENA        - First tree filters from WFILTDT_ODDEVEN 
%      WFILT_ODDEVENB        - Second tree filters from WFILTDT_ODDEVEN 
%      WFILT_OPTSYMA         - First tree filters from WFILTDT_OPTSYM 
%      WFILT_OPTSYMB         - Second tree filters from WFILTDT_OPTSYM 
%      WFILT_DDENA           - First tree filters from WFILTDT_DDEN 
%      WFILT_DDENB           - Second tree filters from WFILTDT_DDEN 
%
%   Dual-Tree Filters
%      WFILTDT_QSHIFT        - Kingsbury's quarter-shift filters
%      WFILTDT_OPTSYM        - Optimizatized Symmetric Self-Hilbertian Filters
%      WFILTDT_ODDEVEN       - Kingsbury's symmetric odd and even biorthogonal filters
%      WFILTDT_DDEN          - Double-density dual-tree filters
%
%  For help, bug reports, suggestions etc. please send an email to
%  ltfat-help@lists.sourceforge.net

