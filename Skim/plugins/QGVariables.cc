#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"


// GEN
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"


#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"


//fastjet
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
using namespace fastjet;
using namespace std;

/**
 * EDProducer class to produced the  QG related variables with subclustering
 * Authors: andrea.carlo.marini@cern.ch
 */

class QGVariables : public edm::EDProducer{
   public:
      explicit QGVariables(const edm::ParameterSet&);
      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void produce(edm::Event&, const edm::EventSetup&);

      //std::tuple<int, int, int, float, float, float, float> calcVariables(const reco::Jet *jet, edm::Handle<reco::VertexCollection>& vC, float dR=0.0001);
      map<string,float> calcVariables(float jeta=0.0, float jphi=0.0, float dR=0.00001); // calc Variables from  inputParticles
      map<string,float> calcVariablesPat(const reco::Jet *jet,float dR=0.00001); // calc Variables from  inputParticles
      map<string,float> calcVariablesGen(const reco::GenJet *jet,float dR=0.00001); // calc Variables from  inputParticles

      template <typename T,typename J> void putInEvent(std::string, const edm::Handle< J >&, std::vector<T>*, edm::Event&);

      edm::EDGetTokenT<edm::View<reco::Jet>> 	jetsToken;
      edm::EDGetTokenT<reco::GenJetCollection> 	genjetsToken;

      vector<PseudoJet> input_particles;
      vector<float> dRToProduce{0.01,0.02,0.03,0.04,0.05,0.075,0.10,0.15,0.20,0.25};
      
};

/// Function to put product into event
template <typename T,typename J> void QGVariables::putInEvent(std::string name, const edm::Handle< J >& jets, std::vector<T>* product, edm::Event& iEvent){
  std::auto_ptr<edm::ValueMap<T>> out(new edm::ValueMap<T>());
  typename edm::ValueMap<T>::Filler filler(*out);
  filler.insert(jets, product->begin(), product->end());
  filler.fill();
  iEvent.put(out, name);
  delete product;
}

QGVariables::QGVariables(const edm::ParameterSet& iConfig) :
  jetsToken( 		consumes<edm::View<reco::Jet>>(		iConfig.getParameter<edm::InputTag>("srcJets"))),//slimmedJets
  genjetsToken( 	consumes<reco::GenJetCollection>(		iConfig.getParameter<edm::InputTag>("srcGenJets"))) //slimmedGenJets
{
  for(const auto dR : dRToProduce){
        produces<edm::ValueMap<float>>(Form("axis2_dR_0p%.0f",dR*1000));
        produces<edm::ValueMap<float>>(Form("axis1_dR_0p%.0f",dR*1000));
        produces<edm::ValueMap<int>>  (Form("mult_dR_0p%.0f",dR*1000));
        produces<edm::ValueMap<float>>(Form("ptD_dR_0p%.0f",dR*1000));
        produces<edm::ValueMap<float>>(Form("ptDrLog_dR_0p%.0f",dR*1000));
    }
}


/// Produce qgLikelihood using {mult, ptD, -log(axis2)}
void QGVariables::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //clear
  input_particles.clear();
  edm::Handle<edm::View<reco::Jet>> jets;
  iEvent.getByToken(jetsToken, 		jets);

  edm::Handle<reco::GenJetCollection> genjets;
  iEvent.getByToken(genjetsToken, 		genjets);
  //
  for(const auto dR : dRToProduce)
  {
    std::vector<float>* axis2Product 		= new std::vector<float>;
    std::vector<float>* axis1Product 		= new std::vector<float>;
    std::vector<int>*   multProduct 		= new std::vector<int>;
    std::vector<float>* ptDProduct 		    = new std::vector<float>;
    std::vector<float>* pt_dr_logProduct 	= new std::vector<float>;



    for(auto jet = jets->begin(); jet != jets->end(); ++jet)
    {
        map<string,float> vars = calcVariablesPat(&*jet,dR);

        axis1Product->push_back(vars["axis1"]);
        axis2Product->push_back(vars["axis2"]);
        multProduct ->push_back(vars["mult"]);
        ptDProduct  ->push_back(vars["ptD"]);
        pt_dr_logProduct->push_back(vars["pt_dr_log"]);
    }

    putInEvent(Form("axis2_dR_0p%.0f",dR*1000),        jets, axis2Product, iEvent);
    putInEvent(Form("axis1_dR_0p%.0f",dR*1000),        jets, axis1Product, iEvent);
    putInEvent(Form("mult_dR_0p%.0f",dR*1000),         jets, multProduct,  iEvent);
    putInEvent(Form("ptD_dR_0p%.0f",dR*1000),          jets, ptDProduct,   iEvent);
    putInEvent(Form("ptDrLog_dR_0p%.0f",dR*1000),    jets, pt_dr_logProduct,iEvent);
  }

  // gen loop
  for(const auto dR : dRToProduce)
  {
    std::vector<float>* axis2Product 		= new std::vector<float>;
    std::vector<float>* axis1Product 		= new std::vector<float>;
    std::vector<int>*   multProduct 		= new std::vector<int>;
    std::vector<float>* ptDProduct 		    = new std::vector<float>;
    std::vector<float>* pt_dr_logProduct 	= new std::vector<float>;



    for(auto jet = genjets->begin(); jet != genjets->end(); ++jet)
    {

        map<string,float> vars = calcVariablesGen(&*jet,dR);

        axis1Product->push_back(vars["axis1"]);
        axis2Product->push_back(vars["axis2"]);
        multProduct ->push_back(vars["mult"]);
        ptDProduct  ->push_back(vars["ptD"]);
        pt_dr_logProduct->push_back(vars["pt_dr_log"]);
    }

    putInEvent(Form("Gen_axis2_dR_0p%.0f",dR*1000),      genjets, axis2Product, iEvent);
    putInEvent(Form("Gen_axis1_dR_0p%.0f",dR*1000),      genjets, axis1Product, iEvent);
    putInEvent(Form("Gen_mult_dR_0p%.0f",dR*1000),       genjets, multProduct,  iEvent);
    putInEvent(Form("Gen_ptD_dR_0p%.0f",dR*1000),        genjets, ptDProduct,   iEvent);
    putInEvent(Form("Gen_ptDrLog_dR_0p%.0f",dR*1000),    genjets, pt_dr_logProduct,iEvent);
  }
}


/// Calculation of axis2, mult and ptD
// reco
std::map<string,float> QGVariables::calcVariablesPat(const reco::Jet *jet,float dR){
    //Loop over the jet constituents
    input_particles.clear();
    for(auto daughter : jet->getJetConstituentsQuick()){
            auto part = static_cast<const pat::PackedCandidate*>(daughter);

            if(part->charge()){
                if(!(part->fromPV() > 1 && part->trackHighPurity())) continue;
                //++mult;
                //++cmult;
            } else {
                if(part->pt() < 1.0) continue;
                //++mult;
                //++nmult;
            }

           input_particles.push_back(PseudoJet (part->px(),part->py(),part->pz(),part->energy()) );
    }
    return calcVariables(jet->eta(),jet->phi(),dR);
}

std::map<string,float> QGVariables::calcVariablesGen(const reco::GenJet *jet,float dR){
    //Loop over the jet constituents
    input_particles.clear();
    for(auto daughter : jet->getGenConstituents () ){
           if ( abs(daughter->pdgId() ) == 12 or abs(daughter->pdgId() ) == 14 or abs(daughter->pdgId() ) == 16 ) continue;
           input_particles.push_back(PseudoJet (daughter->px(),daughter->py(),daughter->pz(),daughter->energy()) );
    }
    return calcVariables(jet->eta(),jet->phi(),dR);
}

// gen

map<string,float> QGVariables::calcVariables(float jeta, float jphi,float dR){
    map<string,float> vars;

    JetDefinition ak_def(antikt_algorithm, dR);
    ClusterSequence seq(input_particles, ak_def);
    auto inclusive_jets = sorted_by_pt(seq.inclusive_jets(0.0)); // ptmin
    
    vars["mult"] = inclusive_jets.size();


    float sum_weight = 0., sum_deta = 0., sum_dphi = 0., sum_deta2 = 0., sum_dphi2 = 0., sum_detadphi = 0., sum_pt = 0.;
    float pt_dr_log = 0;

    for(size_t i=0;i<inclusive_jets.size() ;++i)
    {
        float pt= inclusive_jets.at(i).perp();
        float eta=inclusive_jets.at(i).eta();
        float phi=inclusive_jets.at(i).phi();

        float deta = eta - jeta;
        float dphi = reco::deltaPhi(phi, jphi);
        float weight = pt*pt;

        sum_weight += weight;
        sum_pt += pt;
        sum_deta += deta*weight;
        sum_dphi += dphi*weight;
        sum_deta2 += deta*deta*weight;
        sum_detadphi += deta*dphi*weight;
        sum_dphi2 += dphi*dphi*weight;


        float dr = std::sqrt(reco::deltaPhi(phi, jphi)*reco::deltaPhi(phi, jphi) + (eta-jeta)*(eta-jeta)) ;
        if (dr>0.0001) pt_dr_log += std::log(pt/dr);
        else pt_dr_log += std::log(pt/0.0001);
    }

    //Calculate axis2 and ptD
    float a = 0., b = 0., c = 0.;
    float ave_deta = 0., ave_dphi = 0., ave_deta2 = 0., ave_dphi2 = 0.;
    if(sum_weight > 0){
        ave_deta = sum_deta/sum_weight;
        ave_dphi = sum_dphi/sum_weight;
        ave_deta2 = sum_deta2/sum_weight;
        ave_dphi2 = sum_dphi2/sum_weight;
        a = ave_deta2 - ave_deta*ave_deta;                          
        b = ave_dphi2 - ave_dphi*ave_dphi;                          
        c = -(sum_detadphi/sum_weight - ave_deta*ave_dphi);                
    }
    float delta = sqrt(fabs((a-b)*(a-b)+4*c*c));
    float axis2 = (a+b-delta > 0 ?  sqrt(0.5*(a+b-delta)) : 0);
    float axis1 = (a+b+delta > 0 ?  sqrt(0.5*(a+b+delta)) : 0);
    float ptD   = (sum_weight > 0 ? sqrt(sum_weight)/sum_pt : 0);

    if (sum_pt>0 ) pt_dr_log = pt_dr_log / sum_pt;
    else pt_dr_log = -999;

    vars["ptD"] = ptD;
    vars["axis1"] = axis1;
    vars["axis2"] = axis2;
    vars["pt_dr_log"] = pt_dr_log;

    return vars;
}


/// Descriptions method
void QGVariables::fillDescriptions(edm::ConfigurationDescriptions& descriptions){
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("srcJets");
  desc.add<edm::InputTag>("srcRho");
  desc.add<std::string>("jetsLabel");
  descriptions.add("QGVariables",  desc);
}


//define this as a plug-in
DEFINE_FWK_MODULE(QGVariables);
