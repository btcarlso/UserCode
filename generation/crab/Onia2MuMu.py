# for the list of used tags please see:
# https://twiki.cern.ch/twiki/bin/view/CMS/Onia2MuMuSamples

import FWCore.ParameterSet.Config as cms

process = cms.Process("Onia2MuMuPAT")

from HeavyFlavorAnalysis.Onia2MuMu.onia2MuMuPAT_cff import *

onia2MuMuPAT(process, GlobalTag="START42_V14A::All", MC=True, HLT="RECO", Filter=False)

process.source.fileNames = cms.untracked.vstring(
          'file:upsilon_pgun_RECO_test.root'
          )

#file:upsilongun_10.root

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1) )
