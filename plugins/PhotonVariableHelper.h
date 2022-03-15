#ifndef _ELECTRONVARIABLEHELPER_H
#define _ELECTRONVARIABLEHELPER_H

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/Common/interface/ValueMap.h"

//#include "DataFormats/VertexReco/interface/Vertex.h"
//#include "DataFormats/VertexReco/interface/VertexFwd.h"

//#include "DataFormats/L1Trigger/interface/L1EmParticle.h"
//#include "DataFormats/L1Trigger/interface/L1EmParticleFwd.h"
//#include "DataFormats/L1Trigger/interface/EGamma.h"
//#include "DataFormats/HLTReco/interface/TriggerFilterObjectWithRefs.h"

//#include "DataFormats/Math/interface/deltaR.h"

//#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
//#include "DataFormats/Candidate/interface/CandidateFwd.h"
//#include "DataFormats/Candidate/interface/Candidate.h"

#include <DataFormats/PatCandidates/interface/Photon.h>

//#include "DataFormats/EgammaCandidates/interface/Conversion.h"
//#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"

#include "EgammaAnalysis/TnPTreeProducer/plugins/WriteValueMap.h"
#include "EgammaAnalysis/TnPTreeProducer/plugins/isolations.h"

#include "RecoEcal/EgammaCoreTools/interface/EcalTools.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"

#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"

#include "TMath.h"


template <class T>
class PhotonVariableHelper : public edm::EDProducer {
 public:
  explicit PhotonVariableHelper(const edm::ParameterSet & iConfig);
  virtual ~PhotonVariableHelper() ;

  virtual void produce(edm::Event & iEvent, const edm::EventSetup & iSetup) override;

private:
  edm::EDGetTokenT<reco::PhotonCoreCollection> probesToken_;
  edm::EDGetTokenT<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> > > ebRecHitsToken_;
  edm::EDGetTokenT<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> > > eeRecHitsToken_;

  bool isMiniAODformat;
};

template<class T>
PhotonVariableHelper<T>::PhotonVariableHelper(const edm::ParameterSet & iConfig) :
  probesToken_(consumes<reco::PhotonCoreCollection>(iConfig.getParameter<edm::InputTag>("probes"))),
  ebRecHitsToken_(consumes<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> > >(iConfig.getParameter<edm::InputTag>("ebRecHits"))),
  eeRecHitsToken_(consumes<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> > >(iConfig.getParameter<edm::InputTag>("eeRecHits")))
  {

  produces<edm::ValueMap<float>>("phoSMajor");
  produces<edm::ValueMap<float>>("phoSMinor");

}

template<class T>
PhotonVariableHelper<T>::~PhotonVariableHelper()
{}


template<class T>
void PhotonVariableHelper<T>::produce(edm::Event & iEvent, const edm::EventSetup & iSetup) {

  // read input
  //edm::Handle<pat::Photon> probes;
  //iEvent.getByToken(probesToken_, probes);

  //edm::Handle<reco::PhotonCoreCollection>& photonCoreHandle;
  edm::Handle<reco::PhotonCoreCollection> photonCoreHandle;
  iEvent.getByToken(probesToken_, photonCoreHandle);

  //const reco::SuperCluster& phosuperClus = *probes.superCluster();
  //const reco::CaloCluster &phoseedCluster = *phosuperClus.seed();

  //const bool phoiseb = phoseedCluster.hitsAndFractions()[0].first.subdetId() == EcalBarrel;

  edm::Handle<EcalRecHitCollection> phoEBRecHits_;
  iEvent.getByToken(ebRecHitsToken_, phoEBRecHits_);
  edm::Handle<EcalRecHitCollection> phoEERecHits_;
  iEvent.getByToken(eeRecHitsToken_, phoEERecHits_);



  // prepare vector for output
  std::vector<float> phoSMajorVals;
  std::vector<float> phoSMinorVals;

  /*
  typename std::vector<T>::const_iterator probe, endprobes = probes->end();

  for (probe = probes->begin(); probe != endprobes; ++probe) {
  */

  const EcalRecHitCollection* hits = nullptr;
  for (unsigned int PhotonSize = 0; PhotonSize < photonCoreHandle->size(); PhotonSize++){
    reco::PhotonCoreRef phoCoreRef(reco::PhotonCoreRef(photonCoreHandle, PhotonSize));
    reco::SuperClusterRef scRef = phoCoreRef->superCluster();

    int subdet = scRef->seed()->hitsAndFractions()[0].first.subdetId();
    if (subdet == EcalBarrel) {
      hits = phoEBRecHits_.product();
    }
    else if (subdet == EcalEndcap) {
      hits = phoEERecHits_.product();
    }
    else {
      edm::LogWarning("") << "GEDPhotonProducer: do not know if it is a barrel or endcap SuperCluster: " << subdet;
    }

    if (hits) {
      Cluster2ndMoments moments = EcalClusterTools::cluster2ndMoments(*(scRef->seed()), *hits);

      phoSMajorVals.push_back(moments.sMaj);
      phoSMinorVals.push_back(moments.sMin);
    }

  }


/*
  for (edm::View<pat::Photon>::const_iterator probe = probes->begin(); probe != probes->end(); ++probe) {

    Cluster2ndMoments moments = EcalClusterTools::cluster2ndMoments(*(scRef->seed()), *hits);

    phoSMajorVals.push_back(moments.sMaj);
    phoSMinorVals.push_back(moments.sMin);
  }
  */

  // convert into ValueMap and store
  //writeValueMap(iEvent, probes, phoSMajorVals, "phoSMajor");
  //writeValueMap(iEvent, probes, phoSMinorVals, "phoSMinor");
  writeValueMap(iEvent, photonCoreHandle, phoSMajorVals, "phoSMajor");
  writeValueMap(iEvent, photonCoreHandle, phoSMinorVals, "phoSMinor");

}

#endif
