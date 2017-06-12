#ifndef NERO_LUMI_RUN_H
#define NERO_LUMI_RUN_H

#include "NeroProducer/Nero/interface/Includes.hpp" 
#include "NeroProducer/Core/interface/BareCollection.hpp"

class NeroLumi: virtual public BareCollection{
    public:
        NeroLumi() : BareCollection(){}
        virtual ~NeroLumi(){}
        // --- Just virtual
        virtual int  beginLumi(const edm::LuminosityBlock &,edm::EventSetup const&iSetup) {return 0;}
        virtual int  analyzeLumi(const edm::LuminosityBlock &,edm::EventSetup const&iSetup,TTree*) {return 0;}
        //virtual int  analyzeLumi(const edm::LuminosityBlock &,TH1F*) {return 0;}
        virtual inline string name(){return "NeroLumi";}
};

class NeroRun: virtual public BareCollection{
    public:
        NeroRun() : BareCollection(){}
        virtual ~NeroRun(){}
        // --- Just virtual
        //virtual inline int  analyzeRun(const edm::Run &,TTree*) {return 0;}
        virtual inline int  beginRun(const edm::Run &,edm::EventSetup const&iSetup) {return 0;} 
        virtual inline int  analyzeRun(const edm::Run &,edm::EventSetup const&iSetup,TH1F*) {return 0;} 
        virtual inline string name(){return "NeroRun";}
};

#endif
// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
