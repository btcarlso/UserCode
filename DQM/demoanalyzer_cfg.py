import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.categories.append('Demo')
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
        limit = cms.untracked.int32(-1)
        )
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'file:dcap://cmsdca1.fnal.gov:24140/pnfs/fnal.gov/usr/cms/WAX/11/store/data/Run2012A/SingleMu/AOD/22Jan2013-v1/20000/009E3EE0-716F-E211-81F2-0017A4770C28.root'
    )
)

process.load("Demo.DemoAnalyzer.demoanalyzer_cfi")
process.demo.minTracks=0

process.TFileService = cms.Service("TFileService",
                                                                          fileName = cms.string('histodemo.root')
                                                                      )


#process.dump=cms.EDAnalyzer('EventContentAnalyzer')
#process.Tracer = cms.Service("Tracer")
process.p = cms.Path(process.demo)
#process.p = cms.Path(process.demo*process.dump)
