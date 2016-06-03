#ifndef _LTFAT_ERRNO_H
#define _LTFAT_ERRNO_H

/*
 * The error handling system is adopted from the GSL
 * http://www.gnu.org/software/gsl/manual/html_node/Error-Handling.html#Error-Handling
 * and from
 *
 */

#ifdef __cplusplus
extern "C"
{
#endif

enum
{
    LTFATERR_SUCCESS     =  0,
    LTFATERR_FAILED      =  -1,
    LTFATERR_NULLPOINTER =  -2,
    LTFATERR_NOMEM       =  1,
    LTFATERR_NOTAFRAME   =  3,
    LTFATERR_NOTPAINLESS =  4,
    LTFATERR_NOTPOSARG   =  5,
    LTFATERR_NOTINRANGE  =  6,
    LTFATERR_BADARG      =  7,
    LTFATERR_OVERFLOW    = 10,
    LTFATERR_UNDERFLOW   = 11,
    LTFATERR_CANNOTHAPPEN = 1000,

    LTFATERR_INITFAILED     =  -3
};

void
ltfat_error (const char * reason, const char * file, int line,
             int ltfat_errno);


#ifdef NDEBUG
#define DEBUG(M, ...)
#else
#define DEBUG(M, ...) fprintf(stderr, "DEBUG %s:%d: " M "\n", __FILE__, __LINE__, ##__VA_ARGS__)
#endif

// #define LTFAT_ERROR_PRE(M, ...) fprintf(stderr,"[ERROR] (%s:%d:) " M "\n", __FILE__, __LINE__, ##__VA_ARGS__)
#define LTFAT_ERROR_PRE(M) fprintf(stderr,"[ERROR] (%s:%d:) " M "\n", __FILE__, __LINE__)


#ifdef NCHECKS
#define CHECK(errcode, condition, message, ...)
#define CHECKSTATUS(errcode, message, ...)
#define CHECKMEM(A)
#define CHECKINIT(A, message, ...)
#else

// #define CHECK(errcode, A, M, ...) do{if(!(A)) { LTFAT_ERROR_PRE(M, ##__VA_ARGS__);status=errcode;goto error;}}while(0)
#define CHECK(errcode, A, M ) do{if(!(A)){status=errcode; LTFAT_ERROR_PRE(M); goto error;}}while(0)

// #define CHECKSTATUS(errcode, message, ...) CHECK(errcode,!(errcode),message, ##__VA_ARGS__)
#define CHECKSTATUS(errcode, M) CHECK(errcode,!(errcode), M)

#define CHECKMEM(A) CHECK(LTFATERR_NOMEM,(A), "Out of memory.")

// #define CHECKINIT(A, message, ...) CHECK(LTFATERR_INITFAILED,(A), message, ##__VA_ARGS__)
#define CHECKINIT(A, M) CHECK(LTFATERR_INITFAILED,(A), M)

// #define CHECKCANTHAPPEN(message, ...) do{LTFAT_ERROR_PRE(message, ##__VA_ARGS__); status=LTFATERR_CANNOTHAPPEN; goto error;}while(0)
#define CHECKCANTHAPPEN(M) do{status=LTFATERR_CANNOTHAPPEN; LTFAT_ERROR_PRE(M);  goto error;}while(0)

#endif




#ifdef __cplusplus
}
#endif



#endif
