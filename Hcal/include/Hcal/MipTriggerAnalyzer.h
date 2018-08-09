/**
 * @file MipTriggerAnalyzer.h
 * @brief Header file for MipTriggerAnalyzer class
 * @author Tom Eichlersmith, University of Minnesota
 */

#ifndef HCAL_MIPTRIGGERANALYZER_H
#define HCAL_MIPTRIGGERANALYZER_H

//ROOT
#include "TH1.h" //histograms
#include "TClonesArray.h" //get objects from event bus

//LDMX Framework
#include "Event/Event.h"
#include "Framework/EventProcessor.h" //Needed to declare processor
#include "Framework/ParameterSet.h" // Needed to import parameters from configuration file
#include "Event/TriggerResult.h" //studying triggers

namespace ldmx {
    
    /**
     * @class MipTriggerAnalyzer
     * @brief Creates histogram of the number of tracks found by trigger
     */
    class MipTriggerAnalyzer : public ldmx::Analyzer {
        public:

            MipTriggerAnalyzer(const std::string& name, ldmx::Process& process) : ldmx::Analyzer(name, process) {}

            virtual void configure(const ldmx::ParameterSet& ps);

            virtual void analyze(const ldmx::Event& event);

            virtual void onFileOpen() { }

            virtual void onFileClose() { }

            virtual void onProcessStart(); 

            virtual void onProcessEnd() { }

        private:
            
            /** Name of MipTrigger object */
            std::string hcalMipTriggerObjectName_;

            /** Name of pass that created mip trigger object */
            std::string hcalMipTriggerPassName_;

            /** Number of "tracks" found per event */
            TH1F *hTracksPerEvent_;
    };
}

#endif /* HCAL_MIPTRIGGERANALYZER_H */
