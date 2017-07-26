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

namespace ldmx {

    /**
     * @class kinematicPlots
     * @brief A dummy Analyzer implementation that just prints some messages and makes a simple histogram of calorimeter energy
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
                for (size_t i=0; i<tca->GetEntriesFast(); i++) {
                    ldmx::SimParticle* p=(ldmx::SimParticle*)(tca->At(i));
                    std::cout << p->getEnergy() << std::endl;
                }
            }
            virtual void onFileOpen() {
                std::cout << "kinematicPlots: Opening a file!" << std::endl;
            }
            virtual void onFileClose() {
                std::cout << "kinematicPlots: Closing a file!" << std::endl;
            }
            virtual void onProcessStart() {
                std::cout << "kinematicPlots: Starting processing!" << std::endl;
                getHistoDirectory();
                h_muon_q2=new TH1F("Energy","Energy",500,0.0,10.0);
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
