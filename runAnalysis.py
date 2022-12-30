import FWCore.ParameterSet.Config as cms
import os
import sys
import subprocess

process = cms.Process("PracaMg")

# MessageLogger & co.
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1)
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(False))

dataDir='/eos/user/m/mglazews/Data/'
lsCommand='ls -1 '+dataDir+'|grep root'
print ('command: ',lsCommand)
dir=subprocess.Popen(lsCommand, stdout=subprocess.PIPE,shell=True,text=True)
lsOutput=dir.communicate()[0]
files=[]
for f in lsOutput.split():
  print ('file:'+dataDir+f)
  files.append('file:'+dataDir+f)

print (files)
print (len(files))
  


# input files (up to 255 files accepted)
process.source = cms.Source('PoolSource',
fileNames = cms.untracked.vstring( 
'file:FEVTSIM.root'
),
skipEvents =  cms.untracked.uint32(0),
inputCommands=cms.untracked.vstring( #Added to avoid problematic library for FEVTSIM_1.root file
                             "keep *",
                             "drop l1tTrackerMuons_L1TkMuonsGmt__fullsim"
)
)
#process.source.fileNames = files
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1)) #remember to change to -1! If 1, I was testing something

process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('Configuration.Geometry.GeometryDB_cff')
#process.load('Configuration.Geometry.GeometryExtended2023D17_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
#process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
#from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, '103X_upgrade2023_realistic_v2', '')
#process.load('FWCore.MessageService.MessageLogger_cfi')
#process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1)
#process.MessageLogger.suppressWarning  = cms.untracked.vstring('Geometry', 'AfterSource','L1T','L1GlobalTriggerRawToDigi','omtfTree')
#process.options = cms.untracked.PSet( wantSummary=cms.untracked.bool(False))
 
# PostLS1 geometry used
#process.load('Configuration.Geometry.GeometryExtended2016Reco_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

process.praca= cms.EDAnalyzer("Projekt")
process.MyPath = cms.Path(process.praca)
process.schedule = cms.Schedule(process.MyPath)
