
//warning, lumiunc is kinda a bad name. it's really lumi+jetID+trigger
//really, lumiunc == sqrt( 2.3% ^2 + 1% ^2 + 1% ^2)
const float lumiunc=sqrt(0.023*0.023+2*0.01*0.01);



#include "makeSMPInclJetXsec_strs.h"
#include "makeSMPInclJetXsec_ptbins.h"
#include "makeSMPInclJetXsec_methods.h"



////TO DO: write code in such a way that i can get rid of this stupid const definition that causes me plenty o headaches
//const std::string thyname="NLO_NNPDF_NLO_R04_jtpt";//11.26.19 results
const std::string thyname="NLO_CT14_NLO_R04_jtpt";  //03.18.20 results
const std::string nlothyfile=CT14nnlo;


#include "makeSMPInclJetXsec_data.h"

#include "makeSMPInclJetXsec_PY8.h"
#include "makeSMPInclJetXsec_PY8unf.h"
#include "makeSMPInclJetXsec_NLOandNPCs.h"
#include "makeSMPInclJetXsec_NLOunf.h"











