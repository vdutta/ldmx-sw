/**
 * @file MuonAnalyzer.h
 * @brief Header file for MuonAnalyzer class.
 * @author Tom Eichlersmith, University of Minnesota
 */

#ifndef HCAL_MUONANALYZER_H
#define HCAL_MUONANALYZER_H

//ROOT
#include "TH1.h"
#include "TH2.h"

//LDMX Framework
#include "Event/Event.h"
#include "Framework/EventProcessor.h" //Needed to declare processor
#include "Framework/ParameterSet.h" // Needed to import parameters from configuration file
#include "Event/TriggerResult.h" //testing trigger
#include "Event/HcalHit.h" //Need for path length uncertainty calculation

namespace ldmx {
    
    /**
     * @class MuonAnalyzer
     * @brief Class to study the behavior of muons passing through the Hcal and the trigger trying to 
     *  pick them out.
     */
    class MuonAnalyzer : public ldmx::Analyzer {
        public:

            MuonAnalyzer(const std::string& name, ldmx::Process& process) : ldmx::Analyzer(name, process) {}

            virtual void configure(const ldmx::ParameterSet& ps);

            virtual void analyze(const ldmx::Event& event);

            virtual void onProcessStart(); 

        private:
            
            /** Trigger Object Name */
            std::string triggerObjectName_;

            /** Pass that created trigger */
            std::string triggerPassName_;

            /** HcalHit collection name */
            std::string hcalHitCollName_;

            /** Pass that created HcalHits */
            std::string hcalHitPassName_;

            //One set of histograms for back and side
            TH1F *hNumConsecLayers_[5];
            TH1F *hNumConsecStrips_[5];
            TH2F *hConsecLayersVPathLengthUnc_[2];
            TH2F *hConsecStripsVPathLengthUnc_[2];

            //Rough uncertainty measurement for all muons
            TH1F *hPathLengthUnc_;
            
            //Rough uncertainty measurement for passed muons
            TH1F *hPathLengthUncPassed_;
            
    };
}

#endif /* HCAL_MUONANALYZER_H */
