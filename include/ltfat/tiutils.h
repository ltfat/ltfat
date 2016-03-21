LTFAT_EXTERN void
LTFAT_NAME(fftcircshift)(const LTFAT_COMPLEX* in, int L, double shift,
                         LTFAT_COMPLEX* out);

LTFAT_EXTERN void
LTFAT_NAME(fftfftshift)(const LTFAT_COMPLEX* in, int L, LTFAT_COMPLEX* out);

LTFAT_EXTERN void
LTFAT_NAME(fftifftshift)(const LTFAT_COMPLEX* in, int L, LTFAT_COMPLEX* out);

LTFAT_EXTERN void
LTFAT_NAME(fftrealcircshift)( const LTFAT_COMPLEX* in, int L, double shift,
                              LTFAT_COMPLEX* out);

LTFAT_EXTERN void
LTFAT_NAME(fftrealfftshift)(const LTFAT_COMPLEX* in, int L, LTFAT_COMPLEX* out);

LTFAT_EXTERN void
LTFAT_NAME(fftrealifftshift)(const LTFAT_COMPLEX* in, int L, LTFAT_COMPLEX* out);


