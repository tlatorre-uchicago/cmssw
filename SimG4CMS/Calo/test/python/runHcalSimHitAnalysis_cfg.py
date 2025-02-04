import FWCore.ParameterSet.Config as cms

#from Configuration.Eras.Era_Run3_DDD_cff import Run3_DDD
#process = cms.Process('Dump',Run3_DDD)
from Configuration.Eras.Era_Run3_dd4hep_cff import Run3_dd4hep
process = cms.Process('Dump',Run3_dd4hep)

# import of standard configurations
process.load('FWCore.MessageService.MessageLogger_cfi')
#process.load('Configuration.Geometry.GeometryExtended2021_cff') 
#process.load('Configuration.Geometry.GeometryExtended2021Reco_cff') 
process.load('Configuration.Geometry.GeometryDD4hepExtended2021_cff') 
process.load('Configuration.Geometry.GeometryDD4hepExtended2021Reco_cff') 
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.EventContent.EventContent_cff")
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load("SimG4CMS.Calo.hcalSimHitAnalysis_cfi")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag 
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase1_2021_realistic', '')

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring('file:singlePion_dd4hep.root')
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string('hcalHitdd4hep.root'),
                                   closeFileFast = cms.untracked.bool(True)
                                   )

process.p = cms.Path(process.hcalSimHitAnalysis)
