#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

// Nothing to register:
static const R_CMethodDef CEntries[] = {
    {NULL, NULL, 0}
};

void R_init_ergm_count(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, TRUE);
}
