% LTFAT - Frames
%
%  Peter L. Søndergaard, 2012.
%
%  Basic methods
%    FRAME             - Construct a new frame
%    FRAMEPAIR         - Construct a pair of frames
%    FRAMEDUAL         - The canonical dual frame
%    FRAMETIGHT        - The canonical tight frame
%    FRANA             - Frame analysis
%    FRSYN             - Frame synthesis
%    PLOTFRAME         - Plot frame coefficients
%    FRAMEGRAM         - Plot energy of signal in frame space
%    FRAMEACCEL        - Precompute arrays for faster application
%
%  Information about a frame
%    FRAMEBOUNDS       - Frame bounds
%    FRAMERED          - Redundancy of frame
%    FRAMEMATRIX       - Frame analysis operator matrix
%    FRAMELENGTH       - Length of frame to expand signal
%    FRAMELENGTHCOEF   - Length of frame given a set of coefficients
%
%  Coefficients conversions
%    FRAMECOEF2NATIVE  - Convert to native transform format
%    FRAMENATIVE2COEF  - Convert native to column format
%    FRAMECOEF2TF      - Convert to time-frequency plane layout
%    FRAMETF2COEF      - Convert TF-plane layout to native
%
%  Advanced methods on frames
%    FRSYNABS          - Frame synthesis from magnitude of coefficients
%    FRANAITER         - Iterative inversion of the synthesis frame
%    FRSYNITER         - Iterative inversion of the analysis frame
%    FRAMEMULEIGS      - Eigenpairs of a frame multiplier
%    FRAMELASSO        - LASSO threshholding using Landweber iterations.
%    FRAMEGROUPLASSO   - Group LASSO threshholding.
%
%  For help, bug reports, suggestions etc. please send email to
%  ltfat-help@lists.sourceforge.net
