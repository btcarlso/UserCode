
// -*- C++ -*-
//
// Package:    DemoAnalyzer
// Class:      DemoAnalyzer
// 
/**\class DemoAnalyzer DemoAnalyzer.cc Demo/DemoAnalyzer/src/DemoAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  benjamin carlson
//         Created:  Wed Dec  4 12:51:20 CST 2013
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Common/interface/EventBase.h"
#include "DataFormats/Provenance/interface/EventAuxiliary.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"
#include "Geometry/CaloTopology/interface/CaloSubdetectorTopology.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "Geometry/EcalAlgo/interface/EcalBarrelGeometry.h"
#include "Geometry/EcalAlgo/interface/EcalEndcapGeometry.h"

#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "DataFormats/EcalDetId/interface/ESDetId.h"
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

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"

#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"


#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TH1.h"
#include "TProfile.h"

//
// class declaration
//

class DemoAnalyzer : public edm::EDAnalyzer {
   public:
      explicit DemoAnalyzer(const edm::ParameterSet&);
      ~DemoAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
  // ----member data -----
  unsigned int minTracks_;
  TH1F *demohisto;
  TProfile *bcEB_pv;
  TProfile *bcEE_pv;
  TProfile *scEE_pv;
  TProfile *scEB_pv;
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

      // ----------member data ---------------------------
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
DemoAnalyzer::DemoAnalyzer(const edm::ParameterSet& iConfig) :
  minTracks_(iConfig.getUntrackedParameter<unsigned int>("minTracks",0))
{
   //now do what ever initialization is needed
  edm::Service<TFileService> fs;
  demohisto = fs->make<TH1F>("tracks", "Tracks", 100,0,5000);
  bcEE_pv = fs->make<TProfile>("basicEE_clusters_pv","Basic clusters EE vs. PV; NPV; multi5x5EndcapBasicClusters", 50, 0,50); 
  bcEB_pv = fs->make<TProfile>("basicEB_clusters_pv","Basic Clusters EB vs. PV; NPV; hybridBarrelBasicClusters",50, 0, 50); 
  scEE_pv = fs->make<TProfile>("superEE_clusters_pv","Super Clusters EE vs. PV; NPV; correctedHybridSuperClusters",50, 0, 50);
  scEB_pv = fs->make<TProfile>("superEB_clusters_pv","Super Clusters EB vs. PV; NPV; multi5x5EndcapSuperClusters",50, 0, 50);
}


DemoAnalyzer::~DemoAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
DemoAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   Handle<reco::TrackCollection> tracks;
   iEvent.getByLabel("generalTracks", tracks); 

   edm::Handle<reco::VertexCollection> PVCollection;
   iEvent.getByLabel("offlinePrimaryVertices",PVCollection);
   int nPVcount=0; 
   for (reco::VertexCollection::const_iterator pv = PVCollection->begin(); pv != PVCollection->end(); ++pv) {
              if (pv->isFake() || pv->tracksSize()==0)  continue;
           nPVcount++; // count non fake pv:
   }
   //   LogInfo("Demo") << "Number of PV: " << nPVcount;

   Handle<reco::BasicClusterCollection> EB_BC;
   try {
     iEvent.getByLabel("hybridSuperClusters","hybridBarrelBasicClusters",EB_BC);
   }catch (cms::Exception& ex) {
     edm::LogError("EnergyScaleAnalyzer") << "Can't get collection with producer hybridBarrelBasicClusters";
   }
   //   LogInfo("Demo") << "Number of Basic Clusters EB: " << EB_BC->size();
   bcEB_pv->Fill(nPVcount,EB_BC->size());

   Handle<reco::BasicClusterCollection> EE_BC;
   try {
     iEvent.getByLabel("multi5x5SuperClusters","multi5x5EndcapBasicClusters",EE_BC);
   }catch (cms::Exception& ex) {
     edm::LogError("EnergyScaleAnalyzer") << "Can't get collection with producer multi5x5EndcapBasicClusters";
   }
   //   LogInfo("Demo") << "Number of Basic Clusters EB: " << EB_BC->size();
   bcEE_pv->Fill(nPVcount,EE_BC->size());


    Handle<reco::SuperClusterCollection> correctedHybridSC;
      try {
            iEvent.getByLabel("correctedHybridSuperClusters","",correctedHybridSC);
       }catch (cms::Exception& ex) {
         edm::LogError("EnergyScaleAnalyzer") << "Can't get collection with producer correctedHybridSuperClusters.";
    }
      //      LogInfo("Demo") << "Number of Corrected hybrid SC: " << correctedHybridSC->size();
      
      scEB_pv->Fill(nPVcount,correctedHybridSC->size());
      
      Handle<reco::SuperClusterCollection> SC_EE;
      try {
	iEvent.getByLabel("correctedMulti5x5SuperClustersWithPreshower",SC_EE);
      }catch (cms::Exception& ex) {
	edm::LogError("EnergyScaleAnalyzer") << "Can't get collection with producer correctedHybridSuperClusters.";
      }
      scEE_pv->Fill(nPVcount,SC_EE->size());
      
   demohisto->Fill(tracks->size());
   if( minTracks_ <= tracks->size() ) {
     //LogInfo("Demo") << "number of tracks "<<tracks->size();
   }



#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
}


// ------------ method called once each job just before starting event loop  ------------
void 
DemoAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
DemoAnalyzer::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
DemoAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
DemoAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
DemoAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
DemoAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
DemoAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(DemoAnalyzer);
