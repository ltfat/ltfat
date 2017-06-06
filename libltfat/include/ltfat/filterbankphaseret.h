LTFAT_EXTERN void
LTFAT_NAME(fbmagphasegrad)(const LTFAT_REAL* logs,
                              const LTFAT_REAL* sqtfr,
                              ltfatInt N[], double a[], double fc[], ltfatInt M,
                              ltfatInt neigh[], double posInfo[], LTFAT_REAL gderivweight,
                              LTFAT_REAL tgrad[], LTFAT_REAL fgrad[]);

LTFAT_EXTERN void
LTFAT_NAME(logarray)(const LTFAT_REAL* in, ltfatInt L, LTFAT_REAL* out);
