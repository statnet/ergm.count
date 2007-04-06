#include "MHproposals_bipartite.h" 

/*********************
 void MH_bipartite

 Default MH algorithm for bipartite networks
*********************/
void MH_Bipartiterandomtoggle (MHproposal *MHp, DegreeBound *bd, Network *nwp)  {  
  
  if(MHp->ntoggles == 0) { /* Initialize */
    MHp->ntoggles=1;
    return;
  }
  MHp->ratio = 1.0;
  MHp->togglehead[0] = 1 + unif_rand() * nwp->bipartite;
  MHp->toggletail[0] = 1 + nwp->bipartite + 
                       unif_rand() * (nwp->nnodes - nwp->bipartite);
}

void MH_BipartiteConstantEdges (MHproposal *MHp, DegreeBound *bd, Network *nwp)  {  
  Vertex head, tail;
  int valid;
  
  if(MHp->ntoggles == 0) { /* Initialize */
    MHp->ntoggles=2;    
    return;
  } /* Note:  This proposal cannot be used for full or empty observed graphs.
       If desired, we could check for this at initialization phase. 
       (For now, however, no way to easily return an error message and stop.)*/
  MHp->ratio=1.0;   
  /* First, select edge at random */
  FindithEdge(MHp->togglehead, MHp->toggletail, 1+(nwp->nedges)*unif_rand(), nwp);
  /* Second, select dyad at random until it has no edge */
  valid=0;
   while (valid==0) {
//    head = 1 + unif_rand() * nwp->nnodes;
//    tail = 1 + unif_rand() * nwp->nnodes;
    head = 1 + unif_rand() * nwp->bipartite;
    tail = 1 + nwp->bipartite + 
           unif_rand() * (nwp->nnodes - nwp->bipartite);
    if (EdgetreeSearch(head, tail, nwp->outedges) == 0) {
      valid=1;
    }
   }
  MHp->togglehead[1]=head;
  MHp->toggletail[1]=tail;
}

/********************
   void MH_BipartiteTNT
   Tie/no tie:  Gives at least 50% chance of
   proposing a toggle of an existing edge, as opposed
   to simple random toggles that rarely do so in sparse 
   networks
***********************/
void MH_BipartiteTNT (MHproposal *MHp, DegreeBound *bd, Network *nwp) 
{  
  Vertex head, tail;
  Edge rane, nedges=nwp->nedges;
  int fvalid, trytoggle;
  static double comp=0.5;
  static double odds;
  static Edge ndyads;
  static Edge nnodes;
  static Edge nactors;
  
  if(MHp->ntoggles == 0) { /* Initialize */
    MHp->ntoggles=1;
    odds = comp/(1.0-comp);
    nnodes = nwp[0].nnodes;
    nactors = nwp[0].bipartite;
    ndyads = (nnodes-nactors)*nactors;  
    return;
  }
  
  fvalid = 0;
  trytoggle = 0;
  while(fvalid==0 && trytoggle < 5000){

  if (unif_rand() < comp && nedges > 0) { /* Select a tie at random */
    rane = 1 + unif_rand() * nedges;
    FindithEdge(MHp->togglehead, MHp->toggletail, rane, nwp);
    MHp->ratio = nedges  / (odds*ndyads + nedges);
  }else{ /* Select a dyad at random */
    head = 1 + unif_rand() * nactors;
    tail = 1 + nactors + unif_rand() * (nnodes - nactors);
    MHp->togglehead[0] = head;
    MHp->toggletail[0] = tail;
    if(EdgetreeSearch(MHp->togglehead[0],MHp->toggletail[0],nwp->outedges)!=0){
      MHp->ratio = nedges / (odds*ndyads + nedges);
    }else{
      MHp->ratio = 1.0 + (odds*ndyads)/(nedges + 1);
    }
  }
  fvalid=CheckTogglesValid(MHp, bd, nwp);
  }
}

/********************
   void MH_BipartiteHammingConstantEdges
   Gives at least 50% chance of
   proposing a toggle of an existing edge, as opposed
   to simple random toggles that rarely do so in sparse 
   networks
***********************/
void MH_BipartiteHammingConstantEdges (MHproposal *MHp, DegreeBound *bd, Network *nwp) 
{  
  Vertex head, tail;
  Edge rane, nedges=nwp[0].nedges, nddyads=nwp[1].nedges;
  int valid;
  int nde, ndn, nce, ncn;
  static double comp=0.5;
  static double odds;
  static Edge ndyads;
  static Edge nnodes;
  static Edge nactors;
  
  if(MHp->ntoggles == 0) { /* Initialize */
    MHp->ntoggles=2;
    odds = comp/(1.0-comp);
    nnodes = nwp[0].nnodes;
    nactors = nwp[0].bipartite;
    ndyads = (nnodes-nactors)*nactors;  
    return;
  }
  
  if (unif_rand() < comp && nddyads > 0) { /* Select a discordant pair of tie/nontie at random */
    /* First, select discord edge at random */
    valid=0;
    while (valid==0) {
      rane = 1 + unif_rand() * nddyads;
      FindithEdge(MHp->togglehead, MHp->toggletail, rane, &nwp[1]);
      if (EdgetreeSearch(MHp->togglehead[0], MHp->toggletail[0], nwp[0].outedges) != 0) {
        valid=1;
      }
    }
    head=MHp->togglehead[0];
    tail=MHp->toggletail[0];
    /* Next, select discord non-edge at random */
    valid=0;
    while (valid==0) {
      rane = 1 + unif_rand() * nddyads;
      FindithEdge(MHp->togglehead, MHp->toggletail, rane, &nwp[1]);
      if (EdgetreeSearch(MHp->togglehead[0], MHp->toggletail[0], nwp[0].outedges) == 0) {
        valid=1;
      }
    }
    MHp->togglehead[1]=MHp->togglehead[0];
    MHp->toggletail[1]=MHp->toggletail[0];
    MHp->togglehead[0]=head;
    MHp->toggletail[0]=tail;
    nde = nddyads / 2;
    ndn = nddyads / 2;
    nce = nedges-nde;
    ncn = ndyads-nedges-ndn;
//    MHp->ratio = (nddyads*nddyads) / (odds*(nnodes-nddyads-2)*(nnodes-nddyads-2));
    MHp->ratio = (nde*ndn*1.0) / (odds*(nce+1)*(ncn+1));
//    MHp->ratio = (1.0*(nce-1)*(ncn-1)) / (nde*ndn*odds);
//    MHp->ratio = 1.0;
//   Rprintf("disconcord nde %d nce %d ndn %d ncn %d nddyads %d MHp->ratio %f\n",
//	    nde, nce, ndn,ncn,nddyads, MHp->ratio);
  }else{
    /* First, select concordant edge at random */
    valid=0;
    while (valid==0) {
      rane = 1 + unif_rand() * nedges;
      FindithEdge(MHp->togglehead, MHp->toggletail, rane, &nwp[0]);
      if (EdgetreeSearch(MHp->togglehead[0], MHp->toggletail[0], nwp[1].outedges) == 0) {
        valid=1;
      }
    }
    /* Next, select concord non-edge at random */
    head = 1 + unif_rand() * nactors;
    tail = 1 + nactors + unif_rand() * (nnodes - nactors);
    while((EdgetreeSearch(head,tail,nwp[0].outedges)!=0) ||
          (EdgetreeSearch(head,tail,nwp[1].outedges)!=0)
	 ){
     head = 1 + unif_rand() * nactors;
     tail = 1 + nactors + unif_rand() * (nnodes - nactors);
    }
    MHp->togglehead[1]=head;
    MHp->toggletail[1]=tail;
    nde = nddyads / 2;
    ndn = nddyads / 2;
    nce = nedges-nde;
    ncn = ndyads-nedges-ndn;
//    MHp->ratio = ((nnodes-nddyads)*(nnodes-nddyads)) / (odds*(nddyads+2)*(nddyads+2));
    if(nddyads > 4){
      MHp->ratio = (odds*nce*ncn) / ((nde+1)*(ndn+1)*1.0);
//    MHp->ratio = ((nde+1)*(ndn+1)*odds) / (1.0*nce*ncn);
    }else{
      MHp->ratio = 100000000.0;
    }
//   Rprintf("concord nde %d nce %d ndn %d ncn %d nddyads %d MHp->ratio %f\n",
//	    nde, nce, ndn,ncn,nddyads, MHp->ratio);
  }
//   Rprintf("h0 %d t0 %d h1 %d t1 %d\n", MHp->togglehead[0],  MHp->toggletail[0], 
//                                        MHp->togglehead[1],  MHp->toggletail[1]); 
}

/********************
   void MH_BipartiteHamming
   Gives at least 50% chance of
   proposing a toggle of an existing edge, as opposed
   to simple random toggles that rarely do so in sparse 
   networks
***********************/
void MH_BipartiteHamming (MHproposal *MHp, DegreeBound *bd, Network *nwp) 
{  
  Vertex head, tail;
  Edge rane, nddyads=nwp[1].nedges;
  int nd, nc;
  static double comp=0.5;
  static double odds;
  static Edge ndyads;
  static Edge nnodes;
  static Edge nactors;
  
  if(MHp->ntoggles == 0) { /* Initialize */
    MHp->ntoggles=1;
    odds = comp/(1.0-comp);
    nnodes = nwp[0].nnodes;
    nactors = nwp[0].bipartite;
    ndyads = (nnodes-nactors)*nactors;  
    return;
  }
  
  if (unif_rand() < comp && nddyads > 0) { /* Select a discordant dyad at random */
    rane = 1 + unif_rand() * nddyads;
    FindithEdge(MHp->togglehead, MHp->toggletail, rane, &nwp[1]);
    nd = nddyads;
    nc = ndyads-nd;
    MHp->ratio = (nd*1.0) / (odds*(nc+1));
//    MHp->ratio = (1.0*(nce-1)*(ncn-1)) / (nde*ndn*odds);
//    MHp->ratio = 1.0;
//   Rprintf("disconcord nd %d nc %d nddyads %d MHp->ratio %f\n",
//	    nd, nc, nddyads, MHp->ratio);
  }else{
    /* select a concordant dyad at random */
    head = 1 + unif_rand() * nactors;
    tail = 1 + nactors + unif_rand() * (nnodes - nactors);
    while(EdgetreeSearch(head,tail,nwp[1].outedges)!=0){
     head = 1 + unif_rand() * nactors;
     tail = 1 + nactors + unif_rand() * (nnodes - nactors);
    }
    MHp->togglehead[0]=head;
    MHp->toggletail[0]=tail;
    nd = nddyads;
    nc = ndyads-nd;
    MHp->ratio = (odds*nc) / ((nd+1)*1.0);
//   Rprintf("concord nd %d nc %d nddyads %d MHp->ratio %f\n",
//	    nd, nc, nddyads, MHp->ratio);
  }
//   Rprintf("h0 %d t0 %d h1 %d t1 %d\n", MHp->togglehead[0],  MHp->toggletail[0], 
//                                        MHp->togglehead[1],  MHp->toggletail[1]); 
}


/*********************
 void MH_BipartiteCondDegreeDist
 Pick three nodes -- head, tail, alter -- such that
   * H shares an edge with T
   * H does not share an edge with A
   * deg(T) = deg(A) + 1
 Then, propose swapping the (H,T) edge for a (H,A) edge so that
 deg(H) stays the same while deg(T) and deg(A) swap
 with one another.
*********************/
void MH_BipartiteCondDegreeDist (MHproposal *MHp, DegreeBound *bd, Network *nwp) {  
  int valid, count;
  Edge e;
  Vertex H, T, A, Hin, Hout, Tdeg, Adeg, minA, maxA, i, k;
  double u;
  TreeNode *H_edges;
  Vertex *TA_degrees; 

  if(MHp->ntoggles == 0) { /* Initialize */
    MHp->ntoggles=2;    
    return;
  }
  
  MHp->ratio = 1.0; /* By symmetry:  P(choosing H,T,A) must equal  */
  /*                   P(choosing H,A,T after T and A swap roles), */
  /*                   which makes the ratio equal to 1.           */
  
  for(valid = count = 0; valid == 0 && count<500; count++) {
    Hin = Hout = 0;  
    /* choose a node at random; ensure it has some edges */
    while (Hin + Hout == 0) {
      u = unif_rand();
      if (u < .5) { /* Pick "male" and "female" nodes with equal prob */
        H = 1 + unif_rand() * nwp->bipartite;
      } else {
        H = 1 + nwp->bipartite + unif_rand() * (nwp->nnodes - nwp->bipartite);
      }
      Hin = nwp->indegree[H];
      Hout = nwp->outdegree[H];
    }
    
    /* select an edge to/from H at random */
    k = (int)(unif_rand() * (Hout + Hin));  
    if (k < Hout) { /* we chose an outedge */
      H_edges = nwp->outedges;
      i = k;
      TA_degrees = nwp->indegree;
    } else { /* we chose an inedge */
      H_edges = nwp->inedges;
      i = k - Hout;
      TA_degrees = nwp->outdegree;
    }
    /* Find ith edge in correct edgetree for H; this will be (H,T) or (T,H) */
    e=EdgetreeMinimum(H_edges, H);
    while (i-- > 0) {
      e=EdgetreeSuccessor(H_edges, e);
    } 
    T = H_edges[e].value; 
    Tdeg = nwp->directed_flag ? TA_degrees[T] : nwp->indegree[T] + nwp->outdegree[T];
    
    /* Now choose alter at random */
    /* Now search for eligible alters */
    minA = (u<.5) ? 1 + nwp->bipartite : 1 ;
    maxA = (u<.5) ? nwp->nnodes : nwp->bipartite;    
    i=0;
    for (A=minA; A<=maxA; A++) {
      Adeg = nwp->directed_flag ? TA_degrees[A] : nwp->indegree[A] + nwp->outdegree[A];
      /* To be a valid choice, alter (A) must not be tied to H and it must have */
      /* a degree one less than the tail (T) */
      if (Tdeg == Adeg + 1 && EdgetreeSearch(H, A, H_edges) == 0) {
        if (nwp->directed_flag || EdgetreeSearch(A, H, H_edges) ==0) {
          i++;
        }
      }
    } /* After for-loop, i is # of eligible alters.  */
    if (i>0) {
      valid = 1;
      i = 1 + unif_rand()*i; /* Pick an eligible alter at random */
      for (A=minA; i>0; A++) {
        Adeg = nwp->directed_flag ? TA_degrees[A] : nwp->indegree[A] + nwp->outdegree[A];
        /* To be a valid choice, alter (A) must not be tied to H and it must have */
        /* a degree one less than the tail (T) */
        if (Tdeg == Adeg + 1 && EdgetreeSearch(H, A, H_edges) == 0) {
          if (nwp->directed_flag || EdgetreeSearch(A, H, H_edges) ==0) {
            i--; /* By counting down, when i==0 we have the selected A. */
          }
        }
      } 
    }    
  }
      
  if ( (!nwp->directed_flag && H > T) ||
    (nwp->directed_flag && k < Hout) )
  {
    MHp->togglehead[0] = T;
    MHp->toggletail[0] = H;
  }else{
    MHp->togglehead[0] = H;
    MHp->toggletail[0] = T;
  }

  if(!valid) {
    MHp->togglehead[1] = MHp->togglehead[0];
    MHp->toggletail[1] = MHp->togglehead[0];
  } else {
    if ( (!nwp->directed_flag && H > A) ||
      (nwp->directed_flag && k < Hout) )
    {
      MHp->togglehead[0] = A;
      MHp->toggletail[0] = H;
    }else{
      MHp->togglehead[0] = H;
      MHp->toggletail[0] = A;
    }
  }
}

/********************
   void MH_BipartiteFormation
   Propose ONLY edges not in the reference graph
***********************/
void MH_BipartiteFormation (MHproposal *MHp, DegreeBound *bd, Network *nwp) 
{  
  Vertex head, tail;
  static Edge nnodes;
  static Edge nactors;
  
  if(MHp->ntoggles == 0) { /* Initialize */
    MHp->ntoggles=1;
    nnodes = nwp[0].nnodes;
    nactors = nwp[0].bipartite;
    return;
  }
  
  /* select a dyad not in the reference network at random */
  head = 1 + unif_rand() * nactors;
  tail = 1 + nactors + unif_rand() * (nnodes - nactors);
  while(EdgetreeSearch(head,tail,nwp[0].outedges)!=0 &&
        EdgetreeSearch(head,tail,nwp[1].outedges)==0){
   head = 1 + unif_rand() * nactors;
   tail = 1 + nactors + unif_rand() * (nnodes - nactors);
  }
  MHp->togglehead[0]=head;
  MHp->toggletail[0]=tail;
  MHp->ratio = 1.0;
//   Rprintf("reference nddyads %d MHp->ratio %f\n",
//	    nwp[1].nedges, MHp->ratio);
}

/********************
   void MH_BipartiteFormationTNT
   Tie/no tie:  Gives at least 50% chance of
   proposing a toggle of an existing (non-reference) edge, as opposed
   to simple random toggles that rarely do so in sparse 
   networks
   Propose ONLY edges not in the reference graph
***********************/
void MH_BipartiteFormationTNT (MHproposal *MHp, DegreeBound *bd, Network *nwp) 
{  
  Vertex head, tail;
  Edge rane, nedges, ndedges, nddyads;
  double comp, odds;
  int fvalid, trytoggle;
//  static double comp=0.99;
//  static double odds;
//  static Edge ndedges;
  static Edge ndyads;
  static Edge nnodes;
  static Edge nactors;
  
  if(MHp->ntoggles == 0) { /* Initialize */
    MHp->ntoggles=1;
    nnodes = nwp[0].nnodes;
    nactors = nwp[0].bipartite;
    ndyads = (nnodes-nactors)*nactors;  
    return;
  }
//  Rprintf("nactors %d\n",  nactors);
  
  nedges  = nwp[0].nedges;
  ndedges = nwp[1].nedges;
  if (ndedges > 0) {
    comp = (ndedges*5.0)/(1.0*nedges);
    if(comp > 0.5){comp = 0.5;}
    odds = comp/(1.0-comp);
  }else{
    odds = 0.0;
  }
//  nddyads = ndyads-nedges+ndedges;
//  Rprintf("comp %f nwp[0].nedges %d nwp[1].nedges %d %d %d\n",  comp,
//		       nwp[0].nedges,
//		       nwp[1].nedges, ndedges, nedges);

  fvalid = 0;
  trytoggle = 0;
  while(fvalid==0 && trytoggle < 5000){

  if (ndedges > 0 && unif_rand() < comp) { /* Select a new tie at random */
    rane = 1 + unif_rand() * ndedges;
    FindithEdge(MHp->togglehead, MHp->toggletail, rane, &nwp[1]);
    /* select a dyad not in the reference network at random */
//   Rprintf("nontie h %d t %d nwp[0].nedges %d nwp[1].nedges %d\n",  MHp->togglehead[0],  
//		       MHp->toggletail[0], 
//		       nwp[0].nedges,
//		       nwp[1].nedges);
    MHp->ratio = nedges  / (odds*ndyads + nedges);
  }else{ /* select a dyad not an edge in the reference network at random */
    head = 1 + unif_rand() * nactors;
    tail = 1 + nactors + unif_rand() * (nnodes - nactors);
    while(EdgetreeSearch(head,tail,nwp[0].outedges)!=0 &&
          EdgetreeSearch(head,tail,nwp[1].outedges)==0){
     head = 1 + unif_rand() * nactors;
     tail = 1 + nactors + unif_rand() * (nnodes - nactors);
    }
    MHp->togglehead[0] = head;
    MHp->toggletail[0] = tail;
    if(EdgetreeSearch(MHp->togglehead[0],MHp->toggletail[0],nwp[1].outedges)!=0){
      MHp->ratio = nedges / (odds*ndyads + nedges);
    }else{
      MHp->ratio = 1.0 + (odds*ndyads)/(nedges + 1);
    }
//   Rprintf("nontie h %d t %d nwp[0].nedges %d nwp[1].nedges %d\n",  MHp->togglehead[0],  
//		       MHp->toggletail[0], 
//		       nwp[0].nedges,
//		       nwp[1].nedges);
  }
  fvalid=CheckTogglesValid(MHp, bd, nwp);
  }
}
