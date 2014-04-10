% LTFAT - Filterbanks
%
%  Peter L. Søndergaard, 2011 - 2014
%
%  Transforms and basic routines
%    FILTERBANK             - Filter bank
%    UFILTERBANK            - Uniform Filter bank
%    IFILTERBANK            - Inverse normal/uniform filter bank
%    FILTERBANKWIN          - Evaluate filterbank window
%    FILTERBANKLENGTH       - Length of filterbank to expand signal
%    FILTERBANKLENGTHCOEF   - Length of filterbank to expand coefficients
%
%  Auditory inspired filterbanks
%    CQT                    - Constant Q transform
%    ICQT                   - Inverse constant Q transform
%    ERBLETT                - Erb-let transform
%    IERBLETT               - Inverse Erb-let transform
%    INSDGFB                - Say what?
%
%  Filter generators
%    CQTFILTERS             - Logaritmically spaced filters
%    ERBFILTERS             - ERB-spaced filters
%  
%  Window construction and bounds
%    FILTERBANKDUAL         - Canonical dual filters
%    FILTERBANKTIGHT        - Canonical tight filters
%    FILTERBANKREALDUAL     - Canonical dual filters for real signals
%    FILTERBANKREALTIGHT    - Canonical tight filters for real signals
%    FILTERBANKBOUNDS       - Frame bounds of filter bank
%    FILTERBANKREALBOUNDS   - Frame bounds of filter bank for real signals
%    FILTERBANKRESPONSE     - Total frequency response
%
%  Auxilary
%    NONU2UFILTERBANK       - Non-uni. to uniform filterbank transformation
%    U2NONUCFMT             - Change format of coefficients
%    NONU2UCFMT             - Change format of coefficients back
%
%  Plots
%    PLOTFILTERBANK         - Plot normal/uniform filter bank coefficients
%
%  For help, bug reports, suggestions etc. please send email to
%  ltfat-help@lists.sourceforge.net

