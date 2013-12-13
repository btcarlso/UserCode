/*
 * \file PileUpDep.cc
 * \author Ben Carlson - CMU
 * Last Update:
 *
 */

#include "DQMOffline/Ecal/interface/PileUpDep.h"
#include "DataFormats/Scalers/interface/DcsStatus.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GtFdlWord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"

// Framework

const static int XBINS=2000;

PileUpDep::PileUpDep(const edm::ParameterSet& ps)
{
	//initialize parameters
	basicClusterCollection_EE_ = consumes<reco::BasicClusterCollection>(ps.getParameter<edm::InputTag>("basicClusterCollection_EE"));

}

PileUpDep::~PileUpDep(){
}

void PileUpDep::bookHistograms(DQMStore::IBooker & ibooker,
                                edm::Run const & /* iRun */,
                                edm::EventSetup const & /* iSetup */) {

  // Fetch GlobalTag information and fill the string/ME.
   // initialize
 
}

void PileUpDep::analyze(const edm::Event& e, const edm::EventSetup& c){

	//Analaysis Code
 
	edm::Handle<reco::BasicClusterCollection> basicClusters_EB_h;
	e.getByToken( basicClusterCollection_EB_, basicClusters_EB_h );
//	const reco::BasicClusterCollection* theBarrelBasicClusters = basicClusters_EB_h.product () ;
	if ( ! basicClusters_EB_h.isValid() ) {
		edm::LogWarning("EBRecoSummary") << "basicClusters_EB_h not found"; 
	}
	edm::LogInfo("PileUpOut") << "BasicClusters: " << basicClusters_EB_h->size(); 
	

  return;
}

void
PileUpDep::endLuminosityBlock(const edm::LuminosityBlock& l, const edm::EventSetup& c)
{
  // int nlumi = l.id().luminosityBlock();

  // fill dcs vs lumi
  /* set those bins 0 for which bits are ON
     needed for merge off lumi histograms across files */

  return;
}

DEFINE_FWK_MODULE(PileUpDep);

