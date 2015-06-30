#include "NeroProducer/Nero/interface/NeroTaus.hpp"
#include "NeroProducer/Nero/interface/Nero.hpp"

NeroTaus::NeroTaus(): BareTaus(){

    mMinPt = 18;
    mMinNtaus = 0;
    mMinEta = 2.3;
    mMinId = "decayModeFinding";
    mMaxIso = -1;
}

NeroTaus::~NeroTaus(){
}

int NeroTaus::analyze(const edm::Event & iEvent)
{
    if ( mOnlyMc  ) return 0;

    iEvent.getByToken(token, handle);
    for (const pat::Tau &tau : *handle) {

        if (tau.pt() < 18 ) continue;	
        if (tau.pt() < mMinPt ) continue;	
        
        /// miniaod taus = decayModeFindingNewDMs
        if ( mMinId != "" and !(tau.tauID(mMinId)) ) continue; // minimum requirement to be saved.
        if ( mMaxIso >=0 and tau.tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits") >= mMaxIso ) continue;

        if ( fabs(tau.eta()) > mMinEta ) continue;
        
        // ------------ END SELECTION 

         float phoIso = 0.; for(auto cand : tau.isolationGammaCands() ) phoIso += cand->pt();//tau.isolationPFGammaCandsEtSum() ;
         float chIso  = 0.; for(auto cand : tau.isolationChargedHadrCands() ) chIso += cand->pt();//tau.isolationPFChargedHadrCandsPtSum();
         float nhIso  = 0.; for(auto cand : tau.isolationNeutrHadrCands() ) nhIso += cand->pt(); // PF Cands not exists in miniAOD
         float totIso = phoIso + chIso + nhIso;
        
        //FILL
        new ( (*p4)[p4->GetEntriesFast()]) TLorentzVector(tau.px(), tau.py(), tau.pz(), tau.energy());
        id -> push_back( tau.tauID("decayModeFinding"));
        Q -> push_back( tau.charge() );
        M -> push_back( tau.mass() );
        iso -> push_back( totIso ) ; 

        if (IsExtend() ){

            chargedIsoPtSum  -> push_back( tau.tauID("chargedIsoPtSum") );
            neutralIsoPtSum  -> push_back( tau.tauID("neutralIsoPtSum") );
            isoDeltaBetaCorr -> push_back( tau.tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits"));

            againstEleLoose  -> push_back( tau.tauID("againstElectronLooseMVA5") );
            againstEleMedium -> push_back( tau.tauID("againstElectronMediumMVA5") );
            
            againstMuLoose   -> push_back( tau.tauID("againstMuonLoose3"));
            againstMuTight   -> push_back( tau.tauID("againstMuonTight3"));

            // ---- RC Isolation
            vector<float>  rcIsoTot_v;
            vector<float>  rcIsoCh_v;
            vector<float>  rcIsoNh_v;
            for(int i=1;i<10;++i)
            {
                if (i==5 ) continue; // no back to back
                float phi_rot = 2*TMath::Pi() *i/10;
                float eta= tau.eta();
                float phi = tau.phi() + phi_rot;
                float pt = tau.pt();
                float e = tau.energy();
                TLorentzVector p4_rot; 
                p4_rot.SetPtEtaPhiE(pt,eta,phi,e);

                float rcPhoIso = 0.;
                float rcChIso = 0;
                float rcNhIso = 0;

                float MinPt = 20;
                float MinDrJet = 0.6;
                float Cone = 0.3;
                bool cont= false;
                for(auto j : *(jets_->handle)  ) 
                {
                    if (j.pt() < MinPt ) continue;
                    TLorentzVector jP4(j.px(),j.py(),j.pz(),j.energy() );
                    if (jP4.DeltaR(p4_rot)<MinDrJet )  cont=true;
                    if(cont) break;
                }

                if (cont) continue;

                for(auto cand :  *(pf_->handle) ) 
                {
                    TLorentzVector cP4 ( cand.px(),cand.py(),cand.pz(),cand.energy() );
                    if (cP4.DeltaR(p4_rot) < Cone )
                    {
                     if( cand.charge() != 0 and 
                             abs(cand.pdgId())>20 and 
                             abs(cand.dz())<=0.1 and 
                             cand.fromPV()>1 and 
                             cand.trackHighPurity() 
                             )
                     {rcChIso += cand.pt() ;}
                     if (cand.pdgId() == 22  ) rcPhoIso += cand.pt();
                     if( cand.charge() == 0 and abs(cand.pdgId())>23 ) rcNhIso += cand.pt();

                    }
                }
                float rcTotIso = rcPhoIso + rcChIso + rcNhIso;
                rcIsoTot_v.push_back(rcTotIso);
                rcIsoCh_v.push_back(rcChIso);
                rcIsoNh_v.push_back(rcNhIso+ rcPhoIso);
            } // end loop over angles

            rcIsoTot->push_back(median(rcIsoTot_v));
            rcIsoCh ->push_back(median(rcIsoCh_v));
            rcIsoNh ->push_back(median(rcIsoNh_v));

        } // end Is Extended


    }
    if( int(id->size()) < mMinNtaus) return 1;
    return 0;
}


// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
