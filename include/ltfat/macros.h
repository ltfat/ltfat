#ifndef _LTFAT_MACROS_H
#define _LTFAT_MACROS_H

#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif /* defined(M_PI) */

// To help muting the unused variable compiler warning
// Only works for GCC and Clang
#ifdef __GNUC__
#  define UNUSED(x) UNUSED_ ## x __attribute__((__unused__))
#else
#  define UNUSED(x) UNUSED_ ## x
#endif

// "Vectorizes" a function call
#define LTFAT_APPLYFN(type,fn,...) do{ \
   const type list[] = {(const type)0,__VA_ARGS__}; \
   size_t len = sizeof(list)/sizeof(*list) - 1; \
   for(size_t ii=0;ii<len;ii++) \
      fn((const type)list[ii+1]); \
}while(0)

// Vectorized free
#define LTFAT_SAFEFREEALL(...) LTFAT_APPLYFN(void*,ltfat_safefree,__VA_ARGS__)


#ifdef NDEBUG
#define DEBUG( M, ... )
#define DEBUGNOTE( M )
#else
#define DEBUG(M, ...) fprintf(stderr, "[DEBUG]: (%s:%d:) " M "\n", __FILE__, __LINE__, __VA_ARGS__)
#define DEBUGNOTE(M) fprintf(stderr, "[DEBUG]: (%s:%d:) " M "\n", __FILE__, __LINE__)

#endif


#define CHECK(errcode, A, ...) do{if(!(A)){status=(errcode); ltfat_error(status, __FILE__, __LINE__,__func__ , __VA_ARGS__); goto error;}}while(0)

// The following cannot be just
// #define CHECKSTATUS(errcode, M)  CHECK(errcode,!(errcode), M)
// it evaluates errcode twice!
#define CHECKSTATUS(errcode, M) do{ status = (errcode); CHECK(status,!(status), M);}while(0)
//#define CHECKSTATUSNOMESG(errcode) do{if((errcode)){status=errcode; goto error;}}while(0)

#define CHECKMEM(A) CHECK(LTFATERR_NOMEM,(A), "Out of memory.")
#define CHECKNULL(A) CHECK(LTFATERR_NULLPOINTER,(A), "%s is a null-pointer.",#A)
#define CHECKINIT(A, M) CHECK(LTFATERR_INITFAILED,(A), M)
#define CHECKCANTHAPPEN(M) CHECK(LTFATERR_CANNOTHAPPEN, 0, M)


#endif /* _LTFAT_MACROS_H */
