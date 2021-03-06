#include "NeroProducer/Bambu/interface/NeroMod.h"
#include "NeroProducer/Bambu/interface/TriggerFiller.h"
#include "NeroProducer/Bambu/interface/LeptonsFiller.h"
#include "NeroProducer/Bambu/interface/AllFiller.h"

#include "MitAna/TreeMod/interface/HLTFwkMod.h"
#include "MitAna/DataTree/interface/EventHeaderCol.h"
#include "MitAna/DataTree/interface/MCRunInfo.h"

#include <exception>
#include <cstring>

ClassImp(mithep::NeroMod)

mithep::NeroMod::NeroMod(char const* _name/* = "NeroMod"*/, char const* _title/* = "Nero-Bambu interface"*/) :
  BaseMod(_name, _title)
{
}

mithep::NeroMod::~NeroMod()
{
}

TObjArray*
mithep::NeroMod::GetFillers() const
{
  auto* arr = new TObjArray();
  for (unsigned iC(0); iC != nero::BaseFiller::nCollections; ++iC)
    arr->Add(filler_[iC]);

  return arr;
}

void
mithep::NeroMod::SlaveBegin()
{
  auto* outputFile = TFile::Open(fileName_, "recreate");
  auto* dir = outputFile->mkdir("nero");
  dir->cd();
  eventsTree_ = new TTree("events", "events");
  allTree_ = new TTree("all", "all");
  hXsec_ = new TH1F("xSec", "xSec", 20, -0.5, 19.5);
  hXsec_->Sumw2();

  tag_.Write();
  head_.Write();
  info_.Write();

  nero::BaseFiller::ProductGetter getter([this](char const* _name)->TObject const* {
      if (std::strcmp(_name, mithep::Names::gkMCRunInfoBrn) == 0)
        return this->GetMCRunInfo();

      return this->GetObject<TObject>(_name);
    });

  for (unsigned iC(0); iC != nero::BaseFiller::nCollections; ++iC) {
    if (!filler_[iC])
      continue;

    filler_[iC]->setOutputFile(outputFile);

    filler_[iC]->setProductGetter(getter);

    if (printLevel_ > 0)
      std::cout << "creating branches for " << filler_[iC]->getObject()->name() << std::endl;

    if (iC < nero::BaseFiller::nEventObjects)
      filler_[iC]->defineBranches(eventsTree_);
    else
      filler_[iC]->defineBranches(allTree_);

    if (printLevel_ > 0)
      std::cout << "initializing " << filler_[iC]->getObject()->name() << std::endl;

    filler_[iC]->setCrossRef(filler_);
    filler_[iC]->initialize();
  }
}

void
mithep::NeroMod::BeginRun()
{
  // if other fillers need similar switching, might make sense to define a conditions function on the filler side
  if (filler_[nero::BaseFiller::kTrigger] && (!HasHLTInfo() || !GetHltFwkMod()->HasData()))
    filler_[nero::BaseFiller::kTrigger]->disactivate();

  for (unsigned iC(0); iC != nero::BaseFiller::nCollections; ++iC) {
    if (!filler_[iC])
      continue;

    if (printLevel_ > 0)
      std::cout << "begin run for " << filler_[iC]->getObject()->name() << std::endl;

    filler_[iC]->callBegin();
  }
}

void
mithep::NeroMod::SlaveTerminate()
{
  for (unsigned iC(0); iC != nero::BaseFiller::nCollections; ++iC) {
    if (!filler_[iC])
      continue;

    if (printLevel_ > 0)
      std::cout << "finalizing " << filler_[iC]->getObject()->name() << std::endl;

    filler_[iC]->finalize();
  }

  TFile* outputFile = eventsTree_->GetCurrentFile();
  outputFile->cd("nero");
  eventsTree_->Write();
  allTree_->Write();
  delete outputFile;

  eventsTree_ = allTree_ = 0;
  hXsec_ = 0;
}

Bool_t
mithep::NeroMod::Notify()
{
  for (unsigned iC(0); iC != nero::BaseFiller::nCollections; ++iC) {
    if (!filler_[iC])
      continue;

    filler_[iC]->notify();
  }

  return true;
}

void
mithep::NeroMod::Process()
{
  // always fill allevents
  nero::AllFiller* allFiller = 0;

  if (filler_[nero::BaseFiller::kAll]) {
    allFiller = static_cast<nero::AllFiller*>(filler_[nero::BaseFiller::kAll]);
    if (printLevel_ > 1)
      std::cout << "filling " << allFiller->getObject()->name() << std::endl;

    allFiller->getObject()->clear();
    try {
      allFiller->callFill();
    }
    catch (std::exception& ex) {
      std::cerr << ex.what() << std::endl;
      AbortAnalysis();
    }

    if (printLevel_ > 1)
      std::cout << "filled " << allFiller->getObject()->name() << std::endl;
  }

  allTree_->Fill();

  // skip event if condition module is aborted
  if (condition_ && !condition_->IsActive())
    return;

  for (auto* filler : filler_) {
    if (filler == allFiller)
      continue;

    if (filler) {
      if (printLevel_ > 1)
        std::cout << "filling " << filler->getObject()->name() << std::endl;

      filler->getObject()->clear();
      try {
        filler->callFill();
      }
      catch (std::exception& ex) {
        std::cerr << ex.what() << std::endl;
        AbortAnalysis();
      }

      if (printLevel_ > 1)
        std::cout << "filled " << filler->getObject()->name() << std::endl;
    }
  }

  eventsTree_->Fill();
}
