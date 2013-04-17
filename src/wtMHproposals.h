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



