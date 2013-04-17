/*  File src/wtMHproposals.h in package ergm.count, part of the Statnet suite
 *  of packages for network analysis, http://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  http://statnet.org/attribution
 *
 *  Copyright 2003-2013 Statnet Commons
 */
#ifndef WTMHPROPOSALS_H
#define WTMHPROPOSALS_H

#include "wtMHproposal.h"

void MH_Poisson(WtMHproposal *MHp, WtNetwork *nwp);
void MH_ZIPoisson(WtMHproposal *MHp, WtNetwork *nwp);
void MH_PoissonNonObserved(WtMHproposal *MHp, WtNetwork *nwp);
void MH_Binomial(WtMHproposal *MHp, WtNetwork *nwp);
void MH_BinomialNonObserved(WtMHproposal *MHp, WtNetwork *nwp);
void MH_Geometric(WtMHproposal *MHp, WtNetwork *nwp);
void MH_GeometricNonObserved(WtMHproposal *MHp, WtNetwork *nwp);

#endif 



