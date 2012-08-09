#include "wtMHproposal.h"
#include "R_ext/Rdynload.h"
void WtMH_init(WtMHproposal *MH,char *MHproposaltype, char *MHproposalpackage,double *inputs,int fVerbose,WtNetwork *nwp){
static void (*fun)(WtMHproposal *,char *,char *,double *,int,WtNetwork *) = NULL;
if(fun==NULL) fun = (void (*)(WtMHproposal *,char *,char *,double *,int,WtNetwork *)) R_FindSymbol("WtMH_init", "ergm", NULL);
fun(MH,MHproposaltype,MHproposalpackage,inputs,fVerbose,nwp);
}
void WtMH_free(WtMHproposal *MH){
static void (*fun)(WtMHproposal *) = NULL;
if(fun==NULL) fun = (void (*)(WtMHproposal *)) R_FindSymbol("WtMH_free", "ergm", NULL);
fun(MH);
}
