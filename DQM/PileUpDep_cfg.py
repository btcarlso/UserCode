import FWCore.ParameterSet.Config as cms

process = cms.Process("PileUpOut")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.categories.append('PileUpOut')
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
        limit = cms.untracked.int32(-1)
        )
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
         'file:dcap:///pnfs/cms/WAX/11/store/relval/CMSSW_7_0_0_pre9/RelValZEE/GEN-SIM-RECO/PU_START70_V2-v4/00000/42A074AC-605C-E311-BC4C-003048D37666.root'

    )
)

process.load("DQMOffline.Ecal.PileUpDep_cfi")

process.p = cms.Path(process.pileUp)

