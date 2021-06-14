/*  File src/wtMHproposals.c in package ergm.count, part of the
 *  Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2008-2021 Statnet Commons
 */
#include "wtMHproposals.h"

/* Shorthand. */
#define Mtail (MHp->toggletail)
#define Mhead (MHp->togglehead)
#define Mweight (MHp->toggleweight)

/*********************
 WtMH_P_FN(MH_Poisson

 Default MH algorithm for Poisson-reference ERGM
*********************/
WtMH_P_FN(MH_Poisson){  
  double edgestate;
  
  if(MHp->ntoggles == 0) { // Initialize Poisson 
    MHp->ntoggles=1;
    return;
  }

  GetRandDyad(Mtail, Mhead, nwp);
  
  edgestate = WtGetEdge(Mtail[0],Mhead[0],nwp);

  const double fudge = 0.5; // Mostly comes in when proposing from 0.

  do{
    Mweight[0] = rpois(edgestate + fudge);
  }while(Mweight[0]==edgestate);
    
  MHp->logratio += (1 + log(Mweight[0]+fudge))*edgestate - (1 + log(edgestate+fudge))*Mweight[0] + log(1-dpois(edgestate,edgestate+fudge,0)) - log(1-dpois(Mweight[0],Mweight[0]+fudge,0));
}

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
 WtMH_P_FN(MH_ZIPoisson

 MH algorithm for Poisson-reference ERGM with zero-inflating terms.
 Ocassionally proposes jumps to 0.
*********************/
WtMH_P_FN(MH_PoissonTNT){  
  Edge nedges=nwp->nedges;
  double edgestate;
  static double comp, odds;
  static Dyad ndyads;
  const double fudge = 0.5; // Mostly comes in when proposing from 0.
  
  if(MHp->ntoggles == 0) { // Initialize Poisson 
    MHp->ntoggles=1;
    comp = MHp->inputs[0];
    odds = comp/(1-comp);
    ndyads = DYADCOUNT(nwp->nnodes, nwp->bipartite, nwp->directed_flag);
    return;
  }

  if (unif_rand() < comp && nedges > 0) { /* Select a tie at random */
    WtGetRandEdge(Mtail, Mhead, &edgestate, nwp);
    Mweight[0] = 0;
  }else{ /* Select a dyad at random */
    GetRandDyad(Mtail, Mhead, nwp);
    edgestate = WtGetEdge(Mtail[0],Mhead[0],nwp);
    do{
      Mweight[0] = rpois(edgestate + fudge);
    }while(Mweight[0]==edgestate);
  }

  // Log-probability of a Poisson jump from from to to, given that we didn't stay put:
#define ldpoisj(from,to) (dpois(to, from+fudge, 1) - log1p(-dpois(from, from+fudge, 0)))
#define dpoisj(from,to) exp(dpois(to, from+fudge, 1) - log1p(-dpois(from, from+fudge, 0)))
  
  if(edgestate==0){
    MHp->logratio += log(dpoisj(Mweight[0],edgestate)+odds*ndyads/(nedges+1)) - ldpoisj(edgestate,Mweight[0]) + (nedges==0 ? log(1-comp) : 0);
  }else if(Mweight[0]==0){
    MHp->logratio += ldpoisj(Mweight[0],edgestate) - log(dpoisj(edgestate,Mweight[0])+odds*ndyads/nedges) - (nedges==1 ? log(1-comp) : 0);
  }else{
    MHp->logratio += ldpoisj(Mweight[0],edgestate) - ldpoisj(edgestate,Mweight[0]);
  }
  // h(y)
  MHp->logratio += -lgamma1p(Mweight[0]) - -lgamma1p(edgestate);
  
#undef ldpoisj
#undef dpoisj
}

  
/*********************
 WtMH_P_FN(MH_PoissonNonObserved

 Missing data MH algorithm for Poisson-reference ERGM on bipartite networks.
*********************/
WtMH_P_FN(MH_PoissonNonObserved){  
  Edge nmissing = MHp->inputs[0];

  if(MHp->ntoggles == 0) { /* Initialize */
    MHp->ntoggles=1;
    return;
  }

  if(nmissing==0){
    *Mtail = MH_FAILED;
    *Mhead = MH_IMPOSSIBLE;
    return;
  }


  // Note that missing edgelist is indexed from 0 but the first
  // element of MHp->inputs is the number of missing edges.
  Edge rane = 1 + unif_rand() * nmissing;

  Mtail[0]=MHp->inputs[rane];
  Mhead[1]=MHp->inputs[nmissing+rane];

  double edgestate = WtGetEdge(Mtail[0],Mhead[0],nwp);

  const double fudge = 0.5; // Mostly comes in when proposing from 0.

  do{
    Mweight[0] = rpois(edgestate + fudge);
  }while(Mweight[0]==edgestate);
    
  MHp->logratio += (1 + log(Mweight[0]+fudge))*edgestate - (1 + log(edgestate+fudge))*Mweight[0] + log(1-dpois(edgestate,edgestate+fudge,0)) - log(1-dpois(Mweight[0],Mweight[0]+fudge,0));
}



/*********************
 WtMH_P_FN(MH_Binomial

 Default MH algorithm for binomial-reference ERGM
*********************/
WtMH_P_FN(MH_Binomial){  
  double edgestate;
  static unsigned int trials;
  
  if(MHp->ntoggles == 0) { // Initialize Binomial 
    MHp->ntoggles=1;
    trials = MHp->inputs[0];
    return;
  }
  
  GetRandDyad(Mtail, Mhead, nwp);
  
  edgestate = WtGetEdge(Mtail[0],Mhead[0],nwp);

  const double fudge = 0.5;

  // Use proposal success probability between fudge/trials and (trials-fudge)/trials.
  double p = (fudge + edgestate * (trials-fudge*2)/trials)/trials;

  do{
    Mweight[0] = rbinom(trials, p);
  }while(Mweight[0]==edgestate);

  // What the proposal success probability would have been starting from Mweight[0].
  double revp = (fudge + Mweight[0] * (trials-fudge*2)/trials)/trials;

  MHp->logratio += (lchoose(trials, Mweight[0]) - lchoose(trials, edgestate)) + // h(y*)/h(y)
    (dbinom(edgestate, trials, revp, 1) - log1p(-dbinom(Mweight[0], trials, revp, 0))) - // q(y|y*)
    (dbinom(Mweight[0], trials, p, 1) - log1p(-dbinom(edgestate, trials, p, 0))); // q(y*|y)
}


/*********************
 WtMH_P_FN(MH_BinomialNonObserved

 Missing data MH algorithm for binomial-reference ERGM on bipartite networks.
*********************/
WtMH_P_FN(MH_BinomialNonObserved){  
  Edge nmissing = MHp->inputs[1];

  static unsigned int trials;
  
  if(MHp->ntoggles == 0) { // Initialize Binomial 
    MHp->ntoggles=1;
    trials = MHp->inputs[0];
    return;
  }

  if(nmissing==0){
    *Mtail = MH_FAILED;
    *Mhead = MH_IMPOSSIBLE;
    return;
  }


  // Note that missing edgelist is indexed from 0 but the first
  // element of MHp->inputs is the number of missing edges.
  Edge rane = 1 + unif_rand() * nmissing;

  Mtail[0]=MHp->inputs[rane+1];
  Mhead[1]=MHp->inputs[nmissing+rane+1];

  double edgestate = WtGetEdge(Mtail[0],Mhead[0],nwp);


  const double fudge = 0.5;

  // Use proposal success probability between fudge/trials and (trials-fudge)/trials.
  double p = (fudge + edgestate * (trials-fudge*2)/trials)/trials;

  do{
    Mweight[0] = rbinom(trials, p);
  }while(Mweight[0]==edgestate);

  // What the proposal success probability would have been starting from Mweight[0].
  double revp = (fudge + Mweight[0] * (trials-fudge*2)/trials)/trials;

  MHp->logratio += (lchoose(trials, Mweight[0]) - lchoose(trials, edgestate)) + // h(y*)/h(y)
    (dbinom(edgestate, trials, revp, 1) - log1p(-dbinom(Mweight[0], trials, revp, 0))) - // q(y|y*)
    (dbinom(Mweight[0], trials, p, 1) - log1p(-dbinom(edgestate, trials, p, 0))); // q(y*|y)
}


/*********************
 WtMH_P_FN(MH_Geometric

 Default MH algorithm for geometric-reference ERGM
*********************/
WtMH_P_FN(MH_Geometric){  
  double edgestate;
  
  if(MHp->ntoggles == 0) { // Initialize Geometric 
    MHp->ntoggles=1;
    return;
  }
  
  GetRandDyad(Mtail, Mhead, nwp);
  
  edgestate = WtGetEdge(Mtail[0],Mhead[0],nwp);

  const double fudge = 0.5; // Mostly comes in when proposing from 0.

  double p = 1/(edgestate + 1 + fudge);

  do{
    Mweight[0] = rgeom(p);    
  }while(Mweight[0]==edgestate);

  double revp = 1/(Mweight[0] + 1 + fudge);
    
  MHp->logratio += // h(y) is uniform
    (dgeom(edgestate, revp, 1) - log1p(-dgeom(Mweight[0], revp, 0))) - // q(y|y*)
    (dgeom(Mweight[0], p, 1) - log1p(-dgeom(edgestate, p, 0))); // q(y*|y)
}


/*********************
 WtMH_P_FN(MH_GeometricNonObserved

 Missing data MH algorithm for geometric-reference ERGM.
*********************/
WtMH_P_FN(MH_GeometricNonObserved){  
  Edge nmissing = MHp->inputs[0];

  if(MHp->ntoggles == 0) { /* Initialize */
    MHp->ntoggles=1;
    return;
  }

  if(nmissing==0){
    *Mtail = MH_FAILED;
    *Mhead = MH_IMPOSSIBLE;
    return;
  }


  // Note that missing edgelist is indexed from 0 but the first
  // element of MHp->inputs is the number of missing edges.
  Edge rane = 1 + unif_rand() * nmissing;

  Mtail[0]=MHp->inputs[rane];
  Mhead[1]=MHp->inputs[nmissing+rane];

  double edgestate = WtGetEdge(Mtail[0],Mhead[0],nwp);

  const double fudge = 0.5; // Mostly comes in when proposing from 0.

  double p = 1/(edgestate + fudge);

  do{
    Mweight[0] = rgeom(p);    
  }while(Mweight[0]==edgestate);

  double revp = 1/(Mweight[0] + fudge);
    
  MHp->logratio += // h(y) is uniform
    (dgeom(edgestate, revp, 1) - log1p(-dgeom(Mweight[0], revp, 0))) - // q(y|y*)
    (dgeom(Mweight[0], p, 1) - log1p(-dgeom(edgestate, p, 0))); // q(y*|y)
}
