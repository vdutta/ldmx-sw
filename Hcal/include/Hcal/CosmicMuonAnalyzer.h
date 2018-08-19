/**
 * @file CosmicMuonAnalyzer.h
 * @brief Header file for CosmicMuonAnalyzer class.
 * @author Tom Eichlersmith, University of Minnesota
 */

#ifndef HCAL_COSMICMUONANALYZER_H
#define HCAL_COSMICMUONANALYZER_H

//STL
#include <set>
#include <iterator> //prev

//ROOT
#include "TH1.h"

//LDMX Framework
#include "Event/Event.h"
#include "Framework/EventProcessor.h" //Needed to declare processor
#include "Framework/ParameterSet.h" // Needed to import parameters from configuration file
#include "Event/HcalHit.h"

namespace ldmx {
    
    /**
     * @class CosmicMuonAnalyzer
     * @brief Class to study the behavior of cosmic muons passing through the back Hcal
     */
    class CosmicMuonAnalyzer : public ldmx::Analyzer {
        public:

            CosmicMuonAnalyzer(const std::string& name, ldmx::Process& process) : ldmx::Analyzer(name, process) {}

            virtual void configure(const ldmx::ParameterSet& ps);

            virtual void analyze(const ldmx::Event& event);

            virtual void onProcessStart(); 

        private:
            /**
             * Counts consecutive numbers in list.
             * Returns the maximum number of consecutives.
             */
            int consecutiveCount( const std::set<int> &list ) const;
            
            //One set of histograms for each section
            TH1F *hMinStrip_[5];
            TH1F *hMaxStrip_[5];
            TH1F *hMinLayer_[5];
            TH1F *hMaxLayer_[5];
            TH1F *hNumConsecLayers_[5];
            TH1F *hNumConsecStrips_[5];
    };
}

#endif /* HCAL_COSMICMUONANALYZER_H */
