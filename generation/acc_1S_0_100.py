# ./Ups1S_MuMu-r00-0000018-0000.py with 200000 events, skipEvents = 0
# ######################################################################
import os
import FWCore.ParameterSet.Config as cms

process = cms.Process("HFA")

# ----------------------------------------------------------------------
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.MessageLogger.categories.append('HLTrigReport')
process.MessageLogger.categories.append('L1GtTrigReport')
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )


# ----------------------------------------------------------------------
# -- Database configuration
process.load("CondCore.DBCommon.CondDBCommon_cfi")
process.load("CondCore.DBCommon.CondDBSetup_cfi")

# -- Conditions
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("RecoVertex.BeamSpotProducer.BeamSpot_cfi")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = "START42_V14A::All"
#process.GlobalTag.globaltag = "GR_P_V20::All"

# ----------------------------------------------------------------------
process.source = cms.Source(
 "PoolSource",
  fileNames = cms.untracked.vstring(
    'rfio:/castor/cern.ch/user/b/bora/UpsGun/Ups1S/upsilon_pgun_1S_RECO_NewVtx_11.root',
    )
 
 )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

# ----------------------------------------------------------------------
process.load("HeavyFlavorAnalysis.Bs2MuMu.HFTree_cff")
process.load("HeavyFlavorAnalysis.Bs2MuMu.HFRecoStuff_cff")
process.load("HeavyFlavorAnalysis.Bs2MuMu.HFTruthCandidates_cff")
process.load("HeavyFlavorAnalysis.Bs2MuMu.HFDimuonsCandidates_cff")
process.load("HeavyFlavorAnalysis.Bs2MuMu.HFPhysicsDeclared_cff")
process.load("HeavyFlavorAnalysis.Bs2MuMu.HFSkipEvents_cff")
process.load("HeavyFlavorAnalysis.Bs2MuMu.HFMCTruth_cff")
process.load("HeavyFlavorAnalysis.Bs2MuMu.HFBmm_cff")
process.load("RecoParticleFlow.PFProducer.particleFlow_cfi")
#--------------------------------------------------------
process.hltrep = cms.EDAnalyzer(
    "HLTrigReport",
    HLTriggerResults = cms.InputTag("TriggerResults","","HLT") #REDIGI36X->HLT for collision data
    )


# ----------------------------------------------------------------------
process.HepPDTESSource = cms.ESSource(
    "HepPDTESSource",
    pdtFileName = cms.FileInPath('SimGeneral/HepPDTESSource/data/particle.tbl')
    )
process.genParticles = cms.EDProducer(
    "GenParticleProducer",
    saveBarCodes          = cms.untracked.bool(True),
    src                   = cms.InputTag("generator"),
    abortOnUnknownPDGCode = cms.untracked.bool(False)
    )


# ----------------------------------------------------------------------
process.genDump = cms.EDAnalyzer(
    "HFDumpGenerator",
    generatorCandidates = cms.untracked.string('genParticles'),
    generatorEvent = cms.untracked.string('generator')
    )


# ----------------------------------------------------------------------
try:
    rootFileName = os.environ["JOB"] + ".root"
except KeyError:
    rootFileName = "Acc_1S_0_100_NewVtx_11.root"

process.tree = cms.EDAnalyzer(
    "HFTree",
    verbose      = cms.untracked.int32(0),
    requireCand  =  cms.untracked.bool(False),
    fileName     = cms.untracked.string(rootFileName)
    )


# ----------------------------------------------------------------------
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load("SimTracker.TrackAssociation.TrackAssociatorByHits_cfi")

# ----------------------------------------------------------------------
process.triggerDump = cms.EDAnalyzer(
    "HFDumpTrigger",
    verbose                 = cms.untracked.int32(0),
    L1GTReadoutRecordLabel  = cms.untracked.InputTag("gtDigis"), 
    hltL1GtObjectMap        = cms.untracked.InputTag("hltL1GtObjectMap"), 
    L1MuonsLabel            = cms.untracked.InputTag("hltL1extraParticles"), 
    HLTResultsLabel         = cms.untracked.InputTag("TriggerResults::HLT"), 
    HLTProcessName          = cms.untracked.string('HLT'),
    TriggerEventLabel       = cms.untracked.InputTag("hltTriggerSummaryAOD::HLT"), 
    hltLabel                = cms.untracked.InputTag("TriggerResults::HLT"), 
    )

# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
process.stuffDump = cms.EDAnalyzer(
    "HFDumpStuff",
    verbose                  = cms.untracked.int32(0),
    PrimaryVertexLabel       = cms.untracked.InputTag("offlinePrimaryVertices"),
    PrimaryVertexTracksLabel = cms.untracked.InputTag("generalTracks")
    )


# ----------------------------------------------------------------------
muonDump = cms.EDAnalyzer(
    "HFDumpMuons",
    verbose         = cms.untracked.int32(0),
    muonsLabel      = cms.untracked.InputTag("muons"),
    calomuonsLabel  = cms.untracked.InputTag("calomuons"),
    doTruthMatching = cms.untracked.int32(0),
    runOnAOD        = cms.untracked.bool(True),
    # Configuration for the extrapolation at the muon system 
    propM1 = cms.PSet(
                      useStation2 = cms.bool(False), 
                      useTrack = cms.string("tracker"),
                      useState = cms.string("atVertex"),  # in AOD
                      useSimpleGeometry = cms.bool(True), # use just one cylinder and two planes, not all the fancy chambers
                      ),
    propM2 = cms.PSet(
                      useStation2 = cms.bool(True), 
                      useTrack = cms.string("tracker"),
                      useState = cms.string("atVertex"),  # in AOD
                      useSimpleGeometry = cms.bool(True), # use just one cylinder and two planes, not all the fancy chambers
                      )
    )



# ----------------------------------------------------------------------
process.truthDump = cms.EDAnalyzer(
    "HFTruthCandidate",
    verbose      = cms.untracked.int32(0), 
    tracksLabel  = cms.untracked.InputTag('generalTracks'),
    motherID     = cms.untracked.int32(531),
    type         = cms.untracked.int32(67),
    GenType      = cms.untracked.int32(-67),
    daughtersID  = cms.untracked.vint32(443, 13, -13, 333, 321, -321)
    )

# ----------------------------------------------------------------------
process.p = cms.Path(
#    process.MCTruthSequence*
    process.genDump*
    process.stuffDump*
    process.trkDump*
#    process.muonDump*
#    process.photonDump*
#    process.triggerDump*
#    process.bdpsikstarDump*
    process.bmmDump*
#    process.bmtDump*
#    process.truthDump*
    process.tree
#    process.hltrep
)




