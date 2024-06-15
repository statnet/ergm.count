/*  File src/wtMHproposals.c in package ergm.count, part of the
 *  Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2008-2024 Statnet Commons
 */

#include "ergm_wtMHproposal.h"
#include "ergm_dyadgen.h"
#include "ergm_MHstorage.h"

/*********************
 WtMH_P_FN(MH_ZIPoisson

 MH algorithm for Poisson-reference ERGM with zero-inflating terms.
 Ocassionally proposes jumps to 0.
*********************/
WtMH_P_FN(MH_ZIPoisson){  
  double edgestate, p0=MHp->inputs[0];
  
  if(MHp->ntoggles == 0) { // Initialize Poisson 
    MHp->ntoggles=1;
    return;
  }
  
  GetRandDyad(Mtail, Mhead, nwp);

  edgestate = WtGetEdge(Mtail[0],Mhead[0],nwp);

  const double fudge = 0.5; // Mostly comes in when proposing from 0.

  if(edgestate!=0 && unif_rand()<p0) Mweight[0] = 0;
  else do{
      Mweight[0] = rpois(edgestate + fudge);
    }while(Mweight[0]==edgestate);
 
  // This could probably be done in a numerically-better way:
  // jumping from 0;
  // using (-fudge + fudge*Mweight[0]) in place of dpois(Mweight[0],fudge,1) incorporates the reference measure
  if(edgestate==0)
    MHp->logratio += log(p0+(1-p0)*dpois(0,Mweight[0]+fudge,0)/(1-dpois(Mweight[0],Mweight[0]+fudge,0))) - (-fudge + log(fudge)*Mweight[0]) + log(1-dpois(0,fudge,0)) ;
  else if(Mweight[0]==0)
    MHp->logratio -= log(p0+(1-p0)*dpois(0,edgestate+fudge,0)/(1-dpois(edgestate,edgestate+fudge,0))) - (-fudge + log(fudge)*edgestate) + log(1-dpois(0,fudge,0));
  else // otherwise
    MHp->logratio += (1 + log(Mweight[0]+fudge))*edgestate - (1 + log(edgestate+fudge))*Mweight[0] + log(1-dpois(edgestate,edgestate+fudge,0)) - log(1-dpois(Mweight[0],Mweight[0]+fudge,0)); // Note that (1-p0)s cancel
}


/*********************
Helper functions for discrete references implemented.
*********************/

#define ENSURE_CHANGED(gen) {double to; do{ to = gen; }while(to == from); return to;}

static double rpoisj(double *param, double from, double fudge){
  double mu = from + fudge;
  ENSURE_CHANGED(rpois(mu));
}
static double ldpoisj(double *param, double from, double to, double fudge){
  double mu = from + fudge;
  return(dpois(to, mu, 1) - log1p(-dpois(from, mu, 0)));
}
static double lhpois(double *param, double to){
  return(-lgamma1p(to));
}

static double rgeomj(double *param, double from, double fudge){
  double p = 1/(from + 1 + fudge);
  ENSURE_CHANGED(rgeom(p));
}
static double ldgeomj(double *param, double from, double to, double fudge){
  double p = 1/(from + 1 + fudge);
  return(dgeom(to, p, 1) - log1p(-dgeom(from, p, 0)));
}
static double lhgeom(double *param, double to){
  return(0);
}

static double rbinomj(double *param, double from, double fudge){
  double p = (from + fudge) / (param[0] + fudge*2);
  ENSURE_CHANGED(rbinom(param[0], p));
}
static double ldbinomj(double *param, double from, double to, double fudge){
  double p = (from + fudge) / (param[0] + fudge*2);
  return(dbinom(to, param[0], p, 1) - log1p(-dbinom(from, param[0], p, 0)));
}
static double lhbinom(double *param, double to){
  double result = lchoose(param[0], to);
  return(result);
}

static double rdunifj(double *param, double from, double fudge){
  ENSURE_CHANGED(floor(runif(param[0], param[1]+1)));
}
static double lddunifj(double *param, double from, double to, double fudge){
  return(-log(param[1]-param[0]));
}
static double lhdunif(double *param, double to){
  return(0);
}

/*********************
Disc

MH algorithm for discrete-reference ERGMs.
*********************/

WtMH_I_FN(Mi_Disc){
  MH_STORAGE = DyadGenInitializeR(MHp->R, nwp, FALSE);
  MHp->ntoggles = ((DyadGen *) MH_STORAGE)->ndyads!=0 ? 1 : MH_FAILED;
}


WtMH_P_FN(Mp_Disc){
  const double fudge = 0.5; // Mostly comes in when proposing from 0.
  DyadGen *gen = (DyadGen *) MH_STORAGE;

  double (*rj)(double *param, double from, double fudge),
    (*ldj)(double *param, double from, double to, double fudge),
    (*lh)(double *param, double to);
  switch(MH_IINPUTS[0]){
  case 0: rj = rpoisj; ldj = ldpoisj; lh = lhpois; break;
  case 1: rj = rgeomj; ldj = ldgeomj; lh = lhgeom; break;
  case 2: rj = rbinomj; ldj = ldbinomj; lh = lhbinom; break;
  case 3: rj = rdunifj; ldj = lddunifj; lh = lhdunif; break;
  default: error("Unknown discrete distribution requested.");
  }

  DyadGenRandDyad(Mtail, Mhead, gen);
  double edgestate = WtGetEdge(Mtail[0], Mhead[0], nwp);
  Mweight[0] = rj(MH_DINPUTS, edgestate, fudge);

  double ldjft = ldj(MH_DINPUTS, edgestate, Mweight[0], fudge),
    ldjtf = ldj(MH_DINPUTS, Mweight[0], edgestate, fudge);
  MHp->logratio += ldjtf - ldjft;
  // h(y)
  MHp->logratio += lh(MH_DINPUTS, Mweight[0]) - lh(MH_DINPUTS, edgestate);
}


WtMH_F_FN(Mf_Disc){
  DyadGenDestroy(MH_STORAGE);
  MH_STORAGE = NULL;
}


/*********************
DiscTNT

MH algorithm for discrete-reference ERGMs in sparse networks. Uses TNT weighting.
*********************/

WtMH_I_FN(Mi_DiscTNT){
  MH_STORAGE = DyadGenInitializeR(MHp->R, nwp, TRUE);
  MHp->ntoggles = ((DyadGen *) MH_STORAGE)->ndyads!=0 ? 1 : MH_FAILED;
}


WtMH_P_FN(Mp_DiscTNT){
  const double fudge = 0.5; // Mostly comes in when proposing from 0.
  DyadGen *gen = (DyadGen *) MH_STORAGE;

  double P = MH_DINPUTS[0], Q = 1-P;
  double DP = P*gen->ndyads, DO = DP/Q;
  Edge E = DyadGenEdgecount(gen);

  double (*rj)(double *param, double from, double fudge),
    (*ldj)(double *param, double from, double to, double fudge),
    (*lh)(double *param, double to);
  switch(MH_IINPUTS[0]){
  case 0: rj = rpoisj; ldj = ldpoisj; lh = lhpois; break;
  case 1: rj = rgeomj; ldj = ldgeomj; lh = lhgeom; break;
  case 2: rj = rbinomj; ldj = ldbinomj; lh = lhbinom; break;
  case 3: rj = rdunifj; ldj = lddunifj; lh = lhdunif; break;
  default: error("Unknown discrete distribution requested.");
  }

  double edgestate;
  if (unif_rand() < P && E > 0) { /* Select a tie at random from the network of eligibles */
    DyadGenRandWtEdge(Mtail, Mhead, &edgestate, gen);
    edgestate = WtGetEdge(Mtail[0], Mhead[0], nwp);
    Mweight[0] = 0;
  }else{ /* Select a dyad at random */
    DyadGenRandDyad(Mtail, Mhead, gen);
    edgestate = WtGetEdge(Mtail[0], Mhead[0], nwp);
    Mweight[0] = rj(MH_DINPUTS+1, edgestate, fudge);
  }

  double ldjft = ldj(MH_DINPUTS+1, edgestate, Mweight[0], fudge),
    ldjtf = ldj(MH_DINPUTS+1, Mweight[0], edgestate, fudge);
  if(edgestate==0){
    MHp->logratio += log(exp(ldjtf)+DO/(E+1)) - ldjft + (E==0 ? log(Q) : 0);
  }else if(Mweight[0]==0){
    MHp->logratio += ldjtf - log(exp(ldjft)+DO/E) - (E==1 ? log(Q) : 0);
  }else{
    MHp->logratio += ldjtf - ldjft;
  }
  // h(y)
  MHp->logratio += lh(MH_DINPUTS+1, Mweight[0]) - lh(MH_DINPUTS+1, edgestate);
}


WtMH_F_FN(Mf_DiscTNT){
  DyadGenDestroy(MH_STORAGE);
  MH_STORAGE = NULL;
}
