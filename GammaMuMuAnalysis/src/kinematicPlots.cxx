/**
 * @file kinematicPlots.cxx
 * @brief Analysis of Gamma->MuMu event kinematics
 * @author Andrew Whitbeck, fermilab
 */

#include "Framework/EventProcessor.h"
#include "Framework/ParameterSet.h"
#include <iostream>
#include "TH1.h"
#include "Event/Event.h"
#include "Event/SimParticle.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"

namespace ldmx {

    /**
     * @class kinematicPlots
     * @brief analyzer which isolates gamma->mu+mu- sim particles and computes kinematics for validation
     */
    class kinematicPlots : public ldmx::Analyzer {
        public:

            kinematicPlots(const std::string& name, ldmx::Process& process) : ldmx::Analyzer(name, process) {}

            virtual void configure(const ldmx::ParameterSet& ps) {
                particleCol_=ps.getString("simParticleCollection");
            }

            virtual void analyze(const ldmx::Event& event) {
                std::cout << "kinematicPlots: Analyzing an event!" << std::endl;
                const TClonesArray* tca=event.getCollection(particleCol_);
                ldmx::SimParticle* gamma=0,*mu_plus=0,*mu_minus=0;
                for (size_t i=0; i<tca->GetEntriesFast(); i++) {
                    ldmx::SimParticle* p=(ldmx::SimParticle*)(tca->At(i));
                    if( abs(p->getPdgID()) == 22 && p->getDaughterCount() == 2 && abs(p->getDaughter(0)->getPdgID()) == 13 && abs(p->getDaughter(1)->getPdgID()) == 13 ){
                        gamma = p;
                        if( p->getDaughter(0)->getCharge() == 1 ){
                            mu_plus = p->getDaughter(0);
                            mu_minus = p->getDaughter(1);
                        }else{
                            mu_plus = p->getDaughter(1);
                            mu_minus = p->getDaughter(0);
                        }// end if/else block to determine muon charges
                        break;
                    }// end if gamma/mu+/mu-
                }// end for loop over sim particles

                if( gamma == 0 && mu_plus == 0 && mu_minus == 0 ) return;
                TLorentzVector gamma_p4(gamma->getEnergy(),gamma->getMomentum()[0],gamma->getMomentum()[1],gamma->getMomentum()[2]);
                TLorentzVector mu_plus_p4(mu_plus->getEnergy(),mu_plus->getMomentum()[0],mu_plus->getMomentum()[1],mu_plus->getMomentum()[2]);
                TLorentzVector mu_minus_p4(mu_minus->getEnergy(),mu_minus->getMomentum()[0],mu_minus->getMomentum()[1],mu_minus->getMomentum()[2]);
                h_muon_q2->Fill((mu_plus_p4+mu_minus_p4).P());
                h_muonPhoton_q2->Fill((mu_plus_p4+mu_minus_p4-gamma_p4).P());
            }
            virtual void onFileOpen() {
                std::cout << "kinematicPlots: Opening a file!" << std::endl;
            }
            virtual void onFileClose() {
                std::cout << "kinematicPlots: Closing a file!" << std::endl;
                h_muon_q2->Write();
                h_muonPhoton_q2->Write();
            }
            virtual void onProcessStart() {
                std::cout << "kinematicPlots: Starting processing!" << std::endl;
                getHistoDirectory();
                h_muon_q2=new TH1F("Energy","Energy",500,0.0,10.0);
                h_muonPhoton_q2=new TH1F("Energy","Energy",500,0.0,10.0);
            }
            virtual void onProcessEnd() {
                std::cout << "kinematicPlots: Finishing processing!" << std::endl;
            }
        private:
            TH1 *h_muon_q2,*h_muonPhoton_q2;
            std::string particleCol_;
    };
}

DECLARE_ANALYZER_NS(ldmx, kinematicPlots);
