#include "PhotonVariableHelper.h"

#include "DataFormats/PatCandidates/interface/Photon.h"

#include "FWCore/Framework/interface/MakerMacros.h"

typedef PhotonVariableHelper<pat::Photon>     PatPhotonVariableHelper;
DEFINE_FWK_MODULE(PatPhotonVariableHelper);
