% LTFAT - Frames
%
%  Peter L. Søndergaard, 2012 - 2023.
%
%  Creation of a frame object
%    FRAME             - Construct a new frame
%    FRAMEPAIR         - Construct a pair of frames
%    FUSIONFRAME       - Construct a new fusion frame
%    FRAMEDUAL         - The canonical dual frame
%    FRAMETIGHT        - The canonical tight frame
%    FRAMEACCEL        - Precompute arrays for faster application
%
%  Linear operators
%    FRANA             - Frame analysis
%    FRSYN             - Frame synthesis
%    FRSYNMATRIX       - Frame synthesis operator matrix
%    FRGRAMIAN	       - Frame Gramian operator
%    FRAMEOPERATOR     - Frame operator
%    FRAMEDIAG         - Diagonal of frame operator
%    FRANAITER         - Iterative perfect reconstruction analysis
%    FRSYNITER         - Iterative perfect reconstruction synthesis
%
%  Visualization
%    PLOTFRAME         - Plot frame coefficients
%    FRAMEGRAM         - Plot energy of signal in frame space
%
%  Information about a frame
%    FRAMEBOUNDS       - Frame bounds
%    FRAMERED          - Redundancy of frame
%    FRAMELENGTH       - Length of frame to expand signal
%    FRAMELENGTHCOEF   - Length of frame given a set of coefficients
%    FRAMECLENGTH      - Number of coefficients given input signal length
%    FRAMEVECTORNORMS  - Norms of the frame vectors
%
%  Coefficients conversions
%    FRAMECOEF2NATIVE  - Convert to native transform format
%    FRAMENATIVE2COEF  - Convert native to column format
%    FRAMECOEF2TF      - Convert to time-frequency plane layout
%    FRAMETF2COEF      - Convert TF-plane layout to native
%    FRAMECOEF2TFPLOT  - Convert to time-frequency plane layout for plotting
%
%  Non-linear analysis and synthesis
%    FRANABP           - Basis pursuit using the SALSA algorithm.
%    FRANAMP           - Orthogonal matching pursuit
%    FRANALASSO        - LASSO thresholding using Landweber iterations.
%    FRANAGROUPLASSO   - Group LASSO thresholding.
%    FRSYNABS          - Frame synthesis from magnitude of coefficients
%
%  For help, bug reports, suggestions etc. please visit 
%  http://github.com/ltfat/ltfat/issues
