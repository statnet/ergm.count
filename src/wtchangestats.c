/*  File src/wtchangestats.c in package ergm.count, part of the
 *  Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2008-2022 Statnet Commons
 */
#include "wtchangestats.h"

/********************  changestats:   N    ***********/


/*****************
 stat: sumlogfactorial
 Turns a Poisson-reference or geometric-reference ERGM into a Conway-Maxwell-Poisson distribution
*****************/
WtC_CHANGESTAT_FN(c_sumlogfactorial){
  CHANGE_STAT[0] = lgamma1p(weight)-lgamma1p(edgestate);
}

WtC_CHANGESTAT_FN(c_CMB){
  double ntrials = INPUT_PARAM[0];
  CHANGE_STAT[0] = lgamma1p(weight)-lgamma1p(edgestate) + lgamma1p(ntrials-weight)-lgamma1p(ntrials-edgestate);
}

WtC_CHANGESTAT_FN(c_CMB2){
  double ntrials = INPUT_PARAM[0];
  CHANGE_STAT[0] = lgamma1p(weight)-lgamma1p(edgestate);
  CHANGE_STAT[1] = lgamma1p(ntrials-weight)-lgamma1p(ntrials-edgestate);
}
