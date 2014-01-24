#include "wtchangestats.h"

/********************  changestats:   N    ***********/


/*****************
 stat: sumlogfactorial
 Turns a Poisson-reference or geometric-reference ERGM into a Conway-Maxwell-Poisson distribution
*****************/
WtD_CHANGESTAT_FN(d_sumlogfactorial){
  EXEC_THROUGH_TOGGLES({
      CHANGE_STAT[0] += lgamma1p(NEWWT)-lgamma1p(OLDWT);
  });
}

