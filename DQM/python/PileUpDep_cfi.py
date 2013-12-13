import FWCore.ParameterSet.Config as cms
pileUp = cms.EDAnalyzer('PileUpDep',
                           basicClusterCollection_EE = cms.InputTag("multi5x5SuperClusters","multi5x5EndcapBasicClusters"),
			   basicClusterCollection_EB = cms.InputTag("hybridSuperClusters","hybridBarrelBasicClusters"),
                           superClusterCollection_EE = cms.InputTag("correctedMulti5x5SuperClustersWithPreshower"),
			   superClusterCollection_EB = cms.InputTag("correctedHybridSuperClusters"),
                      )

