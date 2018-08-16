/**
 * @file TriggerAnalyzer.h
 * @brief Header file for TriggerAnalyzer class
 * @author Tom Eichlersmith, University of Minnesota
 */

#ifndef HCAL_TRIGGERANALYZER_H
#define HCAL_TRIGGERANALYZER_H

//ROOT
#include "TH1.h" //histograms
#include "TClonesArray.h" //get objects from event bus

//LDMX Framework
#include "Event/Event.h"
#include "Framework/EventProcessor.h" //Needed to declare processor
#include "Framework/ParameterSet.h" // Needed to import parameters from configuration file
#include "Event/TriggerResult.h" //studying triggers
#include "Event/SimParticle.h" //verifying trigger
#include "Event/HcalHit.h" //checking number of hcal hits

namespace ldmx {
    
    /**
     * @class TriggerAnalyzer
     * @brief Creates histogram of the number of tracks found by trigger
     */
    class TriggerAnalyzer : public ldmx::Analyzer {
        public:

            TriggerAnalyzer(const std::string& name, ldmx::Process& process) : ldmx::Analyzer(name, process) {}

            virtual void configure(const ldmx::ParameterSet& ps);

            virtual void analyze(const ldmx::Event& event);

            virtual void onFileOpen() { }

            virtual void onFileClose() { }

            virtual void onProcessStart() { } 

            virtual void onProcessEnd();

        private:
            
            /** Name of Trigger object */
            std::string hcalTriggerObjectName_;

            /** Name of pass that created trigger object */
            std::string hcalTriggerPassName_;

            /** Counters for trigger */
            unsigned int numFalsePass_;
            unsigned int numTruePass_;
            unsigned int numFalseFail_;
            unsigned int numTrueFail_;

   };
}

#endif /* HCAL_TRIGGERANALYZER_H */
