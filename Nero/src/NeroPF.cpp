#include "NeroProducer/Nero/interface/NeroPF.hpp"
#include "NeroProducer/Nero/interface/Nero.hpp"

NeroPF::NeroPF():  BareCollection(),NeroCollection(){}

NeroPF::~NeroPF(){}

int NeroPF::analyze(const edm::Event& iEvent)
{
	if (mOnlyMc) return 0;
	iEvent.getByToken(token, handle);
	return 0;
}
