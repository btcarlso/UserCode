#ifndef PileUpDep_H
#define PileUpDep_H

/*
 * \file PileUpDep.h
 *
 * \author Ben Carlson - CMU
 *
*/

#include <FWCore/Framework/interface/EDAnalyzer.h>
#include <FWCore/Framework/interface/Event.h>
#include <FWCore/Framework/interface/Run.h>
#include <FWCore/Framework/interface/MakerMacros.h>
#include <FWCore/ParameterSet/interface/ParameterSet.h>
#include <FWCore/ParameterSet/interface/Registry.h>
#include <FWCore/ServiceRegistry/interface/Service.h>
#include <FWCore/MessageLogger/interface/MessageLogger.h>

#include <DQMServices/Core/interface/DQMStore.h>
#include <DQMServices/Core/interface/MonitorElement.h>
#include <DQMServices/Core/interface/DQMEDAnalyzer.h>

#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/ESDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/EgammaReco/interface/PreshowerCluster.h"
#include "DataFormats/EgammaReco/interface/PreshowerClusterFwd.h"
#include "DataFormats/EgammaReco/interface/BasicClusterFwd.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalTools.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalSeverityLevelAlgo.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalCleaningAlgo.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalRecHitLess.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"

#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

//DataFormats
#include "DataFormats/Scalers/interface/DcsStatus.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"

class PileUpDep: public DQMEDAnalyzer{

public:

  /// Constructor
  PileUpDep(const edm::ParameterSet& ps);

  /// Destructor
  virtual ~PileUpDep();

protected:

  /// Analyze
  void analyze(const edm::Event& e, const edm::EventSetup& c);
  void bookHistograms(DQMStore::IBooker &, edm::Run const &, edm::EventSetup const &) override;
  void endLuminosityBlock(const edm::LuminosityBlock& l, const edm::EventSetup& c);
	edm::EDGetTokenT<reco::BasicClusterCollection> basicClusterCollection_EB_;
	edm::EDGetTokenT<reco::BasicClusterCollection> basicClusterCollection_EE_;

	edm::EDGetTokenT<reco::SuperClusterCollection> superClusterCollection_EB_;
	edm::EDGetTokenT<reco::SuperClusterCollection> superClusterCollection_EE_;

private:
	
   // histograms

};

#endif
