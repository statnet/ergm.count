#include "wtmodel.h"
#include "R_ext/Rdynload.h"
WtModel* WtModelInitialize(char *fnames, char *sonames, double **inputs,int n_terms){
static WtModel* (*fun)(char *,char *,double **,int) = NULL;
if(fun==NULL) fun = (WtModel* (*)(char *,char *,double **,int)) R_FindSymbol("WtModelInitialize", "ergm", NULL);
return fun(fnames,sonames,inputs,n_terms);
}
void WtModelDestroy(WtModel *m){
static void (*fun)(WtModel *) = NULL;
if(fun==NULL) fun = (void (*)(WtModel *)) R_FindSymbol("WtModelDestroy", "ergm", NULL);
fun(m);
}
void WtChangeStats(unsigned int ntoggles, Vertex *toggletail, Vertex *togglehead, double *toggleweight, WtNetwork *nwp, WtModel *m){
static void (*fun)(unsigned int,Vertex *,Vertex *,double *,WtNetwork *,WtModel *) = NULL;
if(fun==NULL) fun = (void (*)(unsigned int,Vertex *,Vertex *,double *,WtNetwork *,WtModel *)) R_FindSymbol("WtChangeStats", "ergm", NULL);
fun(ntoggles,toggletail,togglehead,toggleweight,nwp,m);
}
