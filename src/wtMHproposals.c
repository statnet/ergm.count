/*  File src/wtMHproposals.c in package ergm.count, part of the Statnet suite
 *  of packages for network analysis, http://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  http://statnet.org/attribution
 *
 *  Copyright 2008-2016 Statnet Commons
 */
/*  File src/wtMHproposals.c in package ergm.count, part of the Statnet suite
 *  of packages for network analysis, http://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  http://statnet.org/attribution
 *
 *  Copyright 2008-2016 Statnet Commons
 */
#include "wtMHproposals.h"

/* Shorthand. */
#define Mtail (MHp->toggletail)
#define Mhead (MHp->togglehead)
#define Mweight (MHp->toggleweight)

/*********************
 void MH_Poisson

 Default MH algorithm for Poisson-reference ERGM
*********************/
void MH_Poisson(WtMHproposal *MHp, WtNetwork *nwp)  {  
  double oldwt;
  
  if(MHp->ntoggles == 0) { // Initialize Poisson 
    MHp->ntoggles=1;
    return;
  }

  GetRandDyad(Mtail, Mhead, nwp);
  
  oldwt = WtGetEdge(Mtail[0],Mhead[0],nwp);

  const double fudge = 0.5; // Mostly comes in when proposing from 0.

  do{
    Mweight[0] = rpois(oldwt + fudge);    
  }while(Mweight[0]==oldwt);
    
  MHp->logratio += (1 + log(Mweight[0]+fudge))*oldwt - (1 + log(oldwt+fudge))*Mweight[0] + log(1-dpois(oldwt,oldwt+fudge,0)) - log(1-dpois(Mweight[0],Mweight[0]+fudge,0));
}

/*********************
 void MH_ZIPoisson

 MH algorithm for Poisson-reference ERGM with zero-inflating terms.
 Ocassionally proposes jumps to 0.
*********************/
void MH_ZIPoisson(WtMHproposal *MHp, WtNetwork *nwp)  {  
  double oldwt, p0=MHp->inputs[0];
  
  if(MHp->ntoggles == 0) { // Initialize Poisson 
    MHp->ntoggles=1;
    return;
  }
  
  GetRandDyad(Mtail, Mhead, nwp);

  oldwt = WtGetEdge(Mtail[0],Mhead[0],nwp);

  const double fudge = 0.5; // Mostly comes in when proposing from 0.

  if(oldwt!=0 && unif_rand()<p0) Mweight[0] = 0;
  else do{
      Mweight[0] = rpois(oldwt + fudge);    
    }while(Mweight[0]==oldwt);
 
  // This could probably be done in a numerically-better way:
  // jumping from 0;
  // using (-fudge + fudge*Mweight[0]) in place of dpois(Mweight[0],fudge,1) incorporates the reference measure
  if(oldwt==0)
    MHp->logratio += log(p0+(1-p0)*dpois(0,Mweight[0]+fudge,0)/(1-dpois(Mweight[0],Mweight[0]+fudge,0))) - (-fudge + log(fudge)*Mweight[0]) + log(1-dpois(0,fudge,0)) ;
  else if(Mweight[0]==0)
    MHp->logratio -= log(p0+(1-p0)*dpois(0,oldwt+fudge,0)/(1-dpois(oldwt,oldwt+fudge,0))) - (-fudge + log(fudge)*oldwt) + log(1-dpois(0,fudge,0));
  else // otherwise
    MHp->logratio += (1 + log(Mweight[0]+fudge))*oldwt - (1 + log(oldwt+fudge))*Mweight[0] + log(1-dpois(oldwt,oldwt+fudge,0)) - log(1-dpois(Mweight[0],Mweight[0]+fudge,0)); // Note that (1-p0)s cancel
}

/*********************
 void MH_ZIPoisson

 MH algorithm for Poisson-reference ERGM with zero-inflating terms.
 Ocassionally proposes jumps to 0.
*********************/
void MH_PoissonTNT(WtMHproposal *MHp, WtNetwork *nwp)  {  
  Edge nedges=nwp->nedges;
  double oldwt;
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
    WtGetRandEdge(Mtail, Mhead, &oldwt, nwp);
    Mweight[0] = 0;
  }else{ /* Select a dyad at random */
    GetRandDyad(Mtail, Mhead, nwp);
    oldwt = WtGetEdge(Mtail[0],Mhead[0],nwp);
    do{
      Mweight[0] = rpois(oldwt + fudge);    
    }while(Mweight[0]==oldwt);
  }

  // Log-probability of a Poisson jump from from to to, given that we didn't stay put:
#define ldpoisj(from,to) (dpois(to, from+fudge, 1) - log1p(-dpois(from, from+fudge, 0)))
#define dpoisj(from,to) exp(dpois(to, from+fudge, 1) - log1p(-dpois(from, from+fudge, 0)))
  
  if(oldwt==0){
    MHp->logratio += log(dpoisj(Mweight[0],oldwt)+odds*ndyads/(nedges+1)) - ldpoisj(oldwt,Mweight[0]) + (nedges==0 ? log(1-comp) : 0);
  }else if(Mweight[0]==0){
    MHp->logratio += ldpoisj(Mweight[0],oldwt) - log(dpoisj(oldwt,Mweight[0])+odds*ndyads/nedges) - (nedges==1 ? log(1-comp) : 0);
  }else{
    MHp->logratio += ldpoisj(Mweight[0],oldwt) - ldpoisj(oldwt,Mweight[0]);
  }
  // h(y)
  MHp->logratio += -lgamma1p(Mweight[0]) - -lgamma1p(oldwt);
  
#undef ldpoisj
#undef dpoisj
}

  
/*********************
 void MH_PoissonNonObserved

 Missing data MH algorithm for Poisson-reference ERGM on bipartite networks.
*********************/
void MH_PoissonNonObserved(WtMHproposal *MHp, WtNetwork *nwp)  {  
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

  double oldwt = WtGetEdge(Mtail[0],Mhead[0],nwp);

  const double fudge = 0.5; // Mostly comes in when proposing from 0.

  do{
    Mweight[0] = rpois(oldwt + fudge);    
  }while(Mweight[0]==oldwt);
    
  MHp->logratio += (1 + log(Mweight[0]+fudge))*oldwt - (1 + log(oldwt+fudge))*Mweight[0] + log(1-dpois(oldwt,oldwt+fudge,0)) - log(1-dpois(Mweight[0],Mweight[0]+fudge,0));
}



/*********************
 void MH_Binomial

 Default MH algorithm for binomial-reference ERGM
*********************/
void MH_Binomial(WtMHproposal *MHp, WtNetwork *nwp)  {  
  double oldwt;
  static unsigned int trials;
  
  if(MHp->ntoggles == 0) { // Initialize Binomial 
    MHp->ntoggles=1;
    trials = MHp->inputs[0];
    return;
  }
  
  GetRandDyad(Mtail, Mhead, nwp);
  
  oldwt = WtGetEdge(Mtail[0],Mhead[0],nwp);

  const double fudge = 0.5;

  // Use proposal success probability between fudge/trials and (trials-fudge)/trials.
  double p = (fudge + oldwt * (trials-fudge*2)/trials)/trials;

  do{
    Mweight[0] = rbinom(trials, p);
  }while(Mweight[0]==oldwt);

  // What the proposal success probability would have been starting from Mweight[0].
  double revp = (fudge + Mweight[0] * (trials-fudge*2)/trials)/trials;

  MHp->logratio += (lchoose(trials, Mweight[0]) - lchoose(trials, oldwt)) + // h(y*)/h(y)
    (dbinom(oldwt, trials, revp, 1) - log1p(-dbinom(Mweight[0], trials, revp, 0))) - // q(y|y*)
    (dbinom(Mweight[0], trials, p, 1) - log1p(-dbinom(oldwt, trials, p, 0))); // q(y*|y)
}


/*********************
 void MH_BinomialNonObserved

 Missing data MH algorithm for binomial-reference ERGM on bipartite networks.
*********************/
void MH_BinomialNonObserved(WtMHproposal *MHp, WtNetwork *nwp)  {  
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

  double oldwt = WtGetEdge(Mtail[0],Mhead[0],nwp);


  const double fudge = 0.5;

  // Use proposal success probability between fudge/trials and (trials-fudge)/trials.
  double p = (fudge + oldwt * (trials-fudge*2)/trials)/trials;

  do{
    Mweight[0] = rbinom(trials, p);
  }while(Mweight[0]==oldwt);

  // What the proposal success probability would have been starting from Mweight[0].
  double revp = (fudge + Mweight[0] * (trials-fudge*2)/trials)/trials;

  MHp->logratio += (lchoose(trials, Mweight[0]) - lchoose(trials, oldwt)) + // h(y*)/h(y)
    (dbinom(oldwt, trials, revp, 1) - log1p(-dbinom(Mweight[0], trials, revp, 0))) - // q(y|y*)
    (dbinom(Mweight[0], trials, p, 1) - log1p(-dbinom(oldwt, trials, p, 0))); // q(y*|y)
}


/*********************
 void MH_Geometric

 Default MH algorithm for geometric-reference ERGM
*********************/
void MH_Geometric(WtMHproposal *MHp, WtNetwork *nwp)  {  
  double oldwt;
  
  if(MHp->ntoggles == 0) { // Initialize Geometric 
    MHp->ntoggles=1;
    return;
  }
  
  GetRandDyad(Mtail, Mhead, nwp);
  
  oldwt = WtGetEdge(Mtail[0],Mhead[0],nwp);

  const double fudge = 0.5; // Mostly comes in when proposing from 0.

  double p = 1/(oldwt + 1 + fudge);

  do{
    Mweight[0] = rgeom(p);    
  }while(Mweight[0]==oldwt);

  double revp = 1/(Mweight[0] + 1 + fudge);
    
  MHp->logratio += // h(y) is uniform
    (dgeom(oldwt, revp, 1) - log1p(-dgeom(Mweight[0], revp, 0))) - // q(y|y*)
    (dgeom(Mweight[0], p, 1) - log1p(-dgeom(oldwt, p, 0))); // q(y*|y)
}


/*********************
 void MH_GeometricNonObserved

 Missing data MH algorithm for geometric-reference ERGM.
*********************/
void MH_GeometricNonObserved(WtMHproposal *MHp, WtNetwork *nwp)  {  
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

  double oldwt = WtGetEdge(Mtail[0],Mhead[0],nwp);

  const double fudge = 0.5; // Mostly comes in when proposing from 0.

  double p = 1/(oldwt + fudge);

  do{
    Mweight[0] = rgeom(p);    
  }while(Mweight[0]==oldwt);

  double revp = 1/(Mweight[0] + fudge);
    
  MHp->logratio += // h(y) is uniform
    (dgeom(oldwt, revp, 1) - log1p(-dgeom(Mweight[0], revp, 0))) - // q(y|y*)
    (dgeom(Mweight[0], p, 1) - log1p(-dgeom(oldwt, p, 0))); // q(y*|y)
}
