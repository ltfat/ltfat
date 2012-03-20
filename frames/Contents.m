% LTFAT - Frames
%
%  Peter L. Soendergaard, 2012.
%
%  Basic methods
%    NEWFRAME          - Construct a new frame
%    FRANA             - Frame analysis
%    FRSYN             - Frame synthesis
%    FRANAADJ          - Frame analysis adjoint operator
%    FRSYNADJ          - Frame synthesis adjoint operator
%    PLOTFRAME         - Plot frame coefficients
%    FRAMEACCEL        - Precompute arrays for faster application
%
%  Information about a frame
%    FRAMETYPE         - Type of frame
%    FRAMEBOUNDS       - Frame bounds
%    FRAMERED          - Redundancy of frame
%    FRAMEMATRIX       - Frame transform matrix
%    FRAMELENGTHSIGNAL - Length of frame to expand signal
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
%    FRSYNITER         - Iterative frame inversion
%    FRAMEMULEIGS      - Eigenpairs of a frame multiplier
%    FRAMELASSO        - LASSO threshholding using Landweber iterations.
%    FRAMEGROUPLASSO   - Group LASSO threshholding.
%
%  For help, bug reports, suggestions etc. please send email to
%  ltfat-help@lists.sourceforge.net
