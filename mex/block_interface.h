#include "mex.h"

#define COMMAND_LENGTH 20

struct biEntry;

typedef struct 
{
char name[COMMAND_LENGTH];
mxArray* var;
void (*setter)(biEntry* obj,mxArray* in);
mxArray* (*getter)(biEntry* obj);
//const mxArray* defaultVal;
}
biEntry;

void defaultSetter(biEntry* obj,mxArray* in);
mxArray* defaultGetter(biEntry* obj);


