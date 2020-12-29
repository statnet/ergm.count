/*  File src/wtchangestats.c in package ergm.count, part of the Statnet suite
 *  of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution
 *
 *  Copyright 2008-2019 Statnet Commons
 */
#include "wtchangestats.h"

/********************  changestats:   N    ***********/


/*****************
 stat: sumlogfactorial
 Turns a Poisson-reference or geometric-reference ERGM into a Conway-Maxwell-Poisson distribution
*****************/
WtC_CHANGESTAT_FN(c_sumlogfactorial){
  CHANGE_STAT[0] = lgamma1p(weight)-lgamma1p(GETWT(tail,head));
}

WtC_CHANGESTAT_FN(c_CMB){
  double oldweight = GETWT(tail,head), ntrials = INPUT_PARAM[0];
  CHANGE_STAT[0] = lgamma1p(weight)-lgamma1p(oldweight) + lgamma1p(ntrials-weight)-lgamma1p(ntrials-oldweight);
}

WtC_CHANGESTAT_FN(c_CMB2){
  double oldweight = GETWT(tail,head), ntrials = INPUT_PARAM[0];
  CHANGE_STAT[0] = lgamma1p(weight)-lgamma1p(oldweight);
  CHANGE_STAT[1] = lgamma1p(ntrials-weight)-lgamma1p(ntrials-oldweight);
}
