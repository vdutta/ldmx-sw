/**
 * @file eichl008DummyAnalyzer.cxx
 * @brief Class that defines a dummy Analyzer implementation that just prints some messages. Using this file to practice adding a new analyzer to the framework.
 * @author Jeremy Mans, University of Minnesota 
 * @author Tom Eichlersmith, University of Minnesota
 */

#ifndef PLAYTEST_EICHL008DUMMYANALYZER_H
#define PLAYTEST_EICHL008DUMMYANALYZER_H

#include "Framework/EventProcessor.h"
#include "Framework/ParameterSet.h"
#include <iostream>
#include "TH1.h"
#include "Event/Event.h"
#include "Event/CalorimeterHit.h"
#include "Event/HcalHit.h"
#include "TClonesArray.h"

namespace ldmx {

    /**
     * @class eichl008DummyAnalyzer
     * @brief A dummy Analyzer implementation that just prints some messages and makes a simple histogram of calorimeter energy
     */
    class eichl008DummyAnalyzer : public ldmx::Analyzer {
        public:

            eichl008DummyAnalyzer(const std::string& name, ldmx::Process& process) : ldmx::Analyzer(name, process) {}

            virtual void configure(const ldmx::ParameterSet& ps);

            virtual void analyze(const ldmx::Event& event);

            virtual void onFileOpen() {
            }

            virtual void onFileClose() {
            }

            virtual void onProcessStart() {
                getHistoDirectory();
                h_pe = new TH1F("h_pe","PE Distribution",500,0.5,500.5);
                h_energyperevent = new TH1F("h_energyperevent", "Energy Per Event Distribution [MeV]",500,0.0,1000);
            }

            virtual void onProcessEnd() {
            }

        private:
            TH1* h_pe;
            TH1* h_energyperevent;
            std::string caloCol_;
            int dropMod_;
            int keepMod_;
    };
}

#endif /* PLAYTEST_EICHL008DUMMYANALYZER_H */

