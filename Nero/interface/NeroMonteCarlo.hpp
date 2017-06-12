#ifndef NERO_MONTECARLO_H
#define NERO_MONTECARLO_H

#include "NeroProducer/Nero/interface/NeroCollection.hpp"
#include "NeroProducer/Core/interface/BareMonteCarlo.hpp"
#include "NeroProducer/Nero/interface/NeroRunLumi.hpp"

// we should use this gup for the Gen Xsec
#include "GeneratorInterface/Core/interface/GenXSecAnalyzer.h"
#include <memory>


class NeroMonteCarlo : virtual public NeroCollection,
    virtual public BareMonteCarlo,
    virtual public NeroRun,
    virtual public NeroLumi
{
    public:
        NeroMonteCarlo();
        ~NeroMonteCarlo();
        int analyze(const edm::Event& iEvent);
        //void defineBranches(TTree *t);
        virtual inline string name(){return "NeroMonteCarlo";};

        // --- specific fuctions
        int crossSection(edm::Run const & iRun, TH1F* h);
        // --- Handle
        edm::Handle<edm::View<pat::PackedGenParticle>> packed_handle;	
        edm::Handle<edm::View<reco::GenParticle> > pruned_handle;
        edm::Handle<GenEventInfoProduct> info_handle;
        edm::Handle<LHEEventProduct> lhe_handle;
        edm::Handle< std::vector<PileupSummaryInfo> > pu_handle;
        edm::Handle<reco::GenJetCollection> jet_handle;
        edm::Handle<GenRunInfoProduct> runinfo_handle; 
        edm::Handle<std::vector<int> > genBHadFlavour_handle;
        edm::Handle<std::vector<int> > genCHadJetIndex_handle;
        edm::Handle<std::vector<int> > genCHadBHadronId_handle;
        edm::Handle<int > genTtbarId_handle;

        // --- Token
        edm::EDGetTokenT<edm::View<pat::PackedGenParticle> > packed_token;
        edm::EDGetTokenT<edm::View<reco::GenParticle> > pruned_token;
        edm::EDGetTokenT<GenEventInfoProduct> info_token;
        edm::EDGetTokenT<LHEEventProduct> lhe_token;
        edm::EDGetTokenT< std::vector<PileupSummaryInfo> > pu_token;
        edm::EDGetTokenT<reco::GenJetCollection> jet_token;
        edm::EDGetTokenT<GenRunInfoProduct> runinfo_token;
        // bhadrons tokens
        edm::EDGetTokenT<std::vector<int> > genBHadFlavour_token;
        edm::EDGetTokenT<std::vector<int> > genCHadJetIndex_token;
        edm::EDGetTokenT<std::vector<int> > genCHadBHadronId_token;
        edm::EDGetTokenT<int > genTtbarId_token;

        // --- configuration
        float mMinGenParticlePt;
        float mMinGenJetPt;
        bool mParticleGun; 
        int isRealData;
       
        // ---  
        template<class T> // T is supposed to be: reco::GenParticles or Packed Gen Particles
        unsigned  ComputeFlags( T & p ) ;


        std::unique_ptr<GenXSecAnalyzer> gen_;
        int analyzeLumi(const edm::LuminosityBlock &iLumi ,edm::EventSetup const&iSetup,TTree*) ;
        int beginRun( edm::Run const &iRun,edm::EventSetup const&iSetup) override; 
        inline int analyzeRun(edm::Run const & iRun,edm::EventSetup const &iSetup, TH1F* h)override;

};

// template specification declaration
template<>
unsigned NeroMonteCarlo::ComputeFlags<const pat::PackedGenParticle>(const pat::PackedGenParticle &p );

// code implementation of templated functions
template<class T>
unsigned NeroMonteCarlo::ComputeFlags(T &p)
{ // this is called in the two loops
    unsigned flag=0;
    if (p.isPromptFinalState() ) flag |= PromptFinalState;
    if (p.isPromptDecayed() ) flag |= PromptDecayed;
    if (p.isDirectPromptTauDecayProductFinalState() ) flag |= DirectPromptTauDecayProductFinalState;
    if (p.isHardProcess() ) flag |= HardProcess;
    if (p.fromHardProcessBeforeFSR() ) flag |= HardProcessBeforeFSR;
    if (p.fromHardProcessDecayed() ) flag |= HardProcessDecayed;
    if (p.isLastCopy() ) flag |= LastCopy;
    if (p.isLastCopyBeforeFSR() ) flag |= LastCopyBeforeFSR;
    return flag;
}


#endif
// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
