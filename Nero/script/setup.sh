#!/bin/bash

# Instruct builder to use a particular CMSSW release
# [CMSSW] CMSSW_7_6_4
# [Options] isData=False is25ns=True is50ns=False
# [fileList] /store/mc/RunIIFall15MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/70000/FA0A72D5-C7B8-E511-8B1D-901B0E6459E0.root 
# [MaxEvents] 5000

function CMSSW_7_6_4 {
	git cms-init
	git cms-addpkg CommonTools/PileupAlgos
	git cms-addpkg CommonTools/Utils
	git cms-addpkg JetMETCorrections/Configuration
	git cms-addpkg JetMETCorrections/Modules
	git cms-addpkg JetMETCorrections/Type1MET
	git cms-addpkg PhysicsTools/PatAlgos
	git cms-addpkg PhysicsTools/PatUtils
	git cms-addpkg RecoMET/METAlgorithms
	git cms-addpkg RecoMET/METProducers
	git cms-addpkg EgammaAnalysis/ElectronTools
	git cms-addpkg RecoJets/JetProducers
	git cms-merge-topic jbran:pileupJetId76X
	git cms-merge-topic amarini:egcorrection76x
	git cms-merge-topic cms-met:metTool76X
	git remote add blinkseb https://github.com/blinkseb/cmssw.git
	git fetch blinkseb
	git cherry-pick 4cca4688ae368bbbef2102e9bdc5bb00f6df959e
	git cms-merge-topic amarini:topic_met
	#git clone git@github.com:zdemirag/NeroProducer.git ## TO REMOVE
	cd RecoJets/JetProducers/data/
	wget https://github.com/jbrands/RecoJets-JetProducers/raw/3dad903ed25d025f68be94d6f781ca957d6f86ac/pileupJetId_76x_Eta0to2p5_BDT.weights.xml.gz
	wget https://github.com/jbrands/RecoJets-JetProducers/raw/3dad903ed25d025f68be94d6f781ca957d6f86ac/pileupJetId_76x_Eta2p5to2p75_BDT.weights.xml.gz
	wget https://github.com/jbrands/RecoJets-JetProducers/raw/3dad903ed25d025f68be94d6f781ca957d6f86ac/pileupJetId_76x_Eta2p75to3_BDT.weights.xml.gz
	wget https://github.com/jbrands/RecoJets-JetProducers/raw/3dad903ed25d025f68be94d6f781ca957d6f86ac/pileupJetId_76x_Eta3to5_BDT.weights.xml.gz
	cd $CMSSW_BASE/src
}

[ "X$1" == "X" ] && $1=$CMSSW_VERSION

$1
