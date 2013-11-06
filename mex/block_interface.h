#include "mex.h"

#define COMMAND_LENGTH 20

struct biEntry;

typedef struct 
{
char name[COMMAND_LENGTH];
mxArray* var;
void (*setter)(struct biEntry* obj,mxArray* in);
mxArray* (*getter)(struct biEntry* obj);
void (*reseter)(struct biEntry* obj, mxArray* in);
}
biEntry;


void defaultSetter(biEntry* obj,mxArray* in);
mxArray* defaultGetter(biEntry* obj);


biEntry* lookupEntry(const char* name, biEntry* dict, size_t dictLen);

/*
 * Custom functions
 * */
mxArray* getSource(biEntry* obj);
mxArray* getEnqBufCount(biEntry* obj);
mxArray* getToPlay(biEntry* obj);
void incPageNo();
