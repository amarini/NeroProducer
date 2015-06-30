#ifndef NERO_PF_H
#define NERO_PF_H

#include "NeroProducer/Nero/interface/NeroCollection.hpp"
#include "NeroProducer/Core/interface/BareCollection.hpp"


class NeroPF : virtual public NeroCollection, virtual public BareCollection
{
    public:
        // this class just takes the handles so the other classes can share it
        NeroPF();
        ~NeroPF();
        // Bare Collection Virtual
        virtual void clear(){};
        virtual void defineBranches(TTree*) {} ;
        virtual void setBranchAddresses(TTree*) {};

        // Nero Collection Virtual
        int analyze(const edm::Event& iEvent);
        virtual inline string name(){return "NeroPF";};

        // --- Handle
        edm::Handle<pat::PackedCandidateCollection> handle;	

        // --- Token
        edm::EDGetTokenT<pat::PackedCandidateCollection> token;

};


#endif
// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
