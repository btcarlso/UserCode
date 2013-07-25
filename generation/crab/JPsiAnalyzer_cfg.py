import FWCore.ParameterSet.Config as cms

process = cms.Process("ANAPAT")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.GeometryExtended_cff')
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#process.GlobalTag.globaltag = "GR_R_42_V14::All"
#process.GlobalTag.globaltag = "FT_R_42_V13A::All"
#process.GlobalTag.globaltag = "GR_P_V20::All"
process.GlobalTag.globaltag = "START42_V14A::All"

#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

process.source = cms.Source("PoolSource",
    duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
    fileNames = cms.untracked.vstring("file:dcap:///pnfs/cms/WAX/11/store/user/btcarlso/UpsilonPATs/bcarlson/UpsilonGun_1S/UpsilonGun_1S_Onia2MuMuPAT/8f65a5ae4dfc7f75f445421b9cc08be9/onia2MuMuPAT_34_1_ZQ5.root"),
)


# filter on good vertex
process.primaryVertexFilter = cms.EDFilter("GoodVertexFilter",
                                           vertexCollection = cms.InputTag('offlinePrimaryVertices'),
                                           minimumNDOF = cms.uint32(4),
                                           maxAbsZ = cms.double(24),	
                                           maxd0 = cms.double(2)	
                                           )

process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 10000

#read the muon scale correction parameters from the DB
process.poolDBESSource = cms.ESSource("PoolDBESSource",
     BlobStreamerName = cms.untracked.string('TBufferBlobStreamingService'),
     DBParameters = cms.PSet(
          messageLevel = cms.untracked.int32(2),
          authenticationPath = cms.untracked.string('/afs/cern.ch/cms/DB/conddb')
      ),
      timetype = cms.untracked.string('runnumber'),
#      connect = cms.string('oracle://cms_orcoff_prod/CMS_COND_31X_PHYSICSTOOLS'), #needs to point to the correct file
      connect = cms.string('frontier://FrontierProd/CMS_COND_31X_PHYSICSTOOLS'),
      toGet = cms.VPSet(cms.PSet(
          record = cms.string('MuScleFitDBobjectRcd'),
          tag = cms.string('MuScleFit_Scale_JPsi_19_invPb_innerTrack') #adjust the tag once a new one is available
      ))
  )


process.demo = cms.EDAnalyzer('JPsiAnalyzerPAT',

    src = cms.InputTag("onia2MuMuPatTrkTrk"),
    srcWithCaloMuons = cms.InputTag("onia2MuMuPatGlbCal"),

    writeTree = cms.bool(True),
    treeFileName = cms.string("onia2MuMu_tree.root"),

    writeDataSet = cms.bool(True),         
    dataSetName = cms.string("dataSet.root"),
    triggersForDataset = cms.vstring(
                                     "HLT_Dimuon5_Upsilon_Barrel_v3",
                                     "HLT_Dimuon5_Upsilon_Barrel_v5",
                                     "HLT_Dimuon7_Upsilon_Barrel_v1",
                                     "HLT_Dimuon9_Upsilon_Barrel_v1",
                                      ),
     
    massMin = cms.double(8.5),
    massMax = cms.double(11.5),
    pTBinRanges = cms.vdouble(0,10,20,30,40,50.0,60, 70.0,80,90, 100.0),
    etaBinRanges = cms.vdouble(0,0.25, 0.5, 0.75, 1.0),
    onlyTheBest = cms.bool(True),		
    applyCuts = cms.bool(True),
    applyExpHitCuts = cms.untracked.bool(False),#used to be false
    applyDiMuonCuts = cms.untracked.bool(True), #Any reason why this should be false?                         
    useBeamSpot = cms.bool(False),
    useCaloMuons = cms.untracked.bool(False),
    removeSignalEvents = cms.untracked.bool(False),
    removeTrueMuons = cms.untracked.bool(False),
    writeOutCandidates = cms.untracked.bool(False),
    massCorrectionMode = cms.int32(0),    # mode 0 no correction,
                                          # mode 1 constant corr,
                                          # mode 2 pt dependent corr,
                                          # mode 3 pt and eta dependent corr
    oniaPDG = cms.int32(553),
    genParticles = cms.InputTag("genMuons"),
    isMC = cms.untracked.bool(True),
    storeAllMCEvents = cms.untracked.bool(True),
    isPromptMC = cms.untracked.bool(True),
                              
    # Configuration for the extrapolation at the muon system 
    propagatorStation1 = cms.PSet(
        useStation2 = cms.bool(False), 
        useTrack = cms.string("tracker"),
        useState = cms.string("atVertex"),  # for AOD
        useSimpleGeometry = cms.bool(True), 
    ),
    propagatorStation2 = cms.PSet(
        useStation2 = cms.bool(True), 
        useTrack = cms.string("tracker"),
        useState = cms.string("atVertex"),  # for AOD
        useSimpleGeometry = cms.bool(True), 
    ),

    # Configuration of trigger matching; used to be RECO                           
    triggerResultsLabel = cms.InputTag("TriggerResults","","RECO"),
                              
    HLTBitNames_SingleMu = cms.vstring(),
    # ONE FILTER NAME PER PATH
    HLTLastFilterNames_SingleMu = cms.vstring(),
                              
    HLTBitNames_DoubleMu = cms.vstring(
                       "HLT_Dimuon5_Upsilon_Barrel_v3",
                       "HLT_Dimuon5_Upsilon_Barrel_v5",
                       "HLT_Dimuon7_Upsilon_Barrel_v1",
                       "HLT_Dimuon9_Upsilon_Barrel_v1",
                              ),

   # ONE FILTER NAME PER PATH    
   HLTLastFilterNames_DoubleMu = cms.vstring(
                  #"hltVertexmumuFilterUpsilonBarrel",#hltVertexmumuFilterDimuon5UpsiLonbarrelv1
                  #"hltVertexmumuFilterUpsilonBarrel",#No version or dimuon5 in HLTBrowser?
                  "hltVertexmumuFilterDimuon5UpsilonBarrelV3",
                  "hltVertexmumuFilterDimuon5UpsilonBarrelV5",
                  "hltVertexmumuFilterDimuon7UpsilonBarrel",
                  "hltVertexmumuFilterDimuon9UpsilonBarrel",
                  #"hltVertexmumuFilterDimuon7UpsilonBarrelV4",
                  #"hltVertexmumuFilterDimuon9UpsilonBarrelV4"
                  ),



    HLTBitNames_MuL2Mu = cms.vstring(),
    # TWO FILTER NAMES PER PATH (FIRST is L3, SECOND is L2)                           
    HLTLastFilterNames_MuL2Mu = cms.vstring(),

    HLTBitNames_MuTrack = cms.vstring(), 
    # ONE FILTER NAME PER PATH
    HLTLastFilterNames_MuTrack = cms.vstring(),
                              
    HLTBitNames_MuTkMu = cms.vstring(),
    # ONE FILTER NAME PER PATH
    HLTLastFilterNames_MuTkMu = cms.vstring(),

)

## no filter
# process.p = cms.Path(process.demo)

## filter on vertex
#process.p = cms.Path(process.primaryVertexFilter*process.demo)
process.p = cms.Path(process.demo)

#import PhysicsTools.PythonAnalysis.LumiList as LumiList
#import FWCore.ParameterSet.Types as CfgTypes
#myLumis = LumiList.LumiList(filename = 'Cert_160404-163869_7TeV_May10ReReco_Collisions11_JSON_MuonPhys.txt').getCMSSWString().split(',')
#myLumis = LumiList.LumiList(filename = 'Cert_160404-178677_7TeV_PromptReco_Collisions11_JSON_MuonPhys.txt').getCMSSWString().split(',')
#process.source.lumisToProcess = CfgTypes.untracked(CfgTypes.VLuminosityBlockRange())
#process.source.lumisToProcess.extend(myLumis)
## filter on vertex and HLT
# process.p = cms.Path(process.primaryVertexFilter*process.hltMuF*process.demo)
