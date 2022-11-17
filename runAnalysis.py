import FWCore.ParameterSet.Config as cms
import os
import sys
import subprocess

process = cms.Process("PracaMg")

# MessageLogger & co.
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1)
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(False))

#Setting data directory
dataDir='/scratch/CMSSW_12_3_0/src/Projekt/Magisterka/Data'
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
#'/store/mc/Phase2HLTTDRSummer20ReRECOMiniAOD/GluGluHToGG_M125_14TeV_powheg_pythia8_TuneCP5//GEN-SIM-DIGI-RAW-MINIAOD/NoPU_111X_mcRun4_realistic_T15_v1-v1/270000/073371ED-128E-E741-A5D4-6C3C87E14820.root',
'file:FEVTSIM.root' 
),
skipEvents =  cms.untracked.uint32(0)
)
#process.source.fileNames = files
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1))

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

process.praca= cms.EDAnalyzer("Projekt")

process.MyPath = cms.Path(process.praca)
process.schedule = cms.Schedule(process.MyPath)
