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
%  Advanced methods on frames
%    FRSYNABS          - Frame synthesis from magnitude of coefficients
%    FRSYNITER         - Iterative frame inversion
%    FRAMEMULEIGS      - Eigenpairs of a frame multiplier
%
%  For help, bug reports, suggestions etc. please send email to
%  ltfat-help@lists.sourceforge.net
