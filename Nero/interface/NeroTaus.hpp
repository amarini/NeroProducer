#ifndef NERO_TAUS_H
#define NERO_TAUS_H

#include "NeroProducer/Nero/interface/NeroCollection.hpp"
#include "NeroProducer/Core/interface/BareTaus.hpp"
#include "NeroProducer/Nero/interface/NeroJets.hpp"
#include "NeroProducer/Nero/interface/NeroPF.hpp"

class NeroTaus : virtual public NeroCollection,
    virtual public BareTaus
{
    public:
        NeroTaus();
        ~NeroTaus();
        int analyze(const edm::Event &)  ;
        virtual inline string name(){return "NeroTaus";};


        // Token
        edm::EDGetTokenT<pat::TauCollection> token ;	
        // Handle
        edm::Handle<pat::TauCollection> handle;

        // --- configuration
        float mMinPt;
        int   mMinNtaus;
        float mMinEta;
        string mMinId;
        float mMaxIso;
        
        // ---- RC
        NeroJets *jets_;
        NeroPF  *pf_;
        inline float mean(vector<float>&a){ float S=0; for(auto& x : a) S+=x; S/=a.size() ; return S;}
        inline float rms(vector<float>&a){ float S=0; float m=mean(a); for(auto& x: a) S+= (x-m)*(x-m); S/=(a.size()-1); S=TMath::Sqrt(S) ;return S;}
        inline float median(vector<float>&a){ sort(a.begin(),a.end() );return a[ (a.size()-1)/2 ] ; }
};

#endif

// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
