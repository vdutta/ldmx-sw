/**
 * @file HcalLayerAnalyzer.h
 * @brief Class implementation of analyzer to study layer differences in HCal
 * @author Tom Eichlersmith, University of Minnesota
 */

#ifndef PLAYTEST_HCALLAYERANALYZER_H
#define PLAYTEST_HCALLAYERANALYZER_H

//Standard Libraries
#include <iostream> //Checks to std::cout for development purposes
#include <vector> //Vector of hits per layer
#include <algorithm> //Sort hitlog vectors

//ROOT
#include "TH1.h" //One Dimensional Histograms
#include "TClonesArray.h" //Array to get list of hits from Event

//LDMX Framework
#include "Framework/EventProcessor.h" //Needed to declare analyzer
#include "Framework/ParameterSet.h" // Needed to import parameters from configuration file
#include "Event/Event.h" //Study event by event
#include "Event/HcalHit.h" //Get Hcal Specific information from hit
#include "DetDescr/HcalID.h" //Get HcalSection enum for filtering out Hcal Hits

namespace ldmx {

    /**
     * @class HcalLayerAnalyzer
     * @brief ldmx::Analyzer that constructs histograms studying how layers in the Hcal behave differently.
     *
     * @note Currently, the Hcal strip orientation is not specified (x or y), and simulation is only split
     *  along y (as of June 18, 2018). This means any orientation-related studies in this analyzer will be
     *  making external assumptions that SHOULD BE REMOVED if strip orientation is specifed in the future.
     */
    class HcalLayerAnalyzer : public ldmx::Analyzer {
        public:

            HcalLayerAnalyzer(const std::string& name, ldmx::Process& process) : ldmx::Analyzer(name, process) {}

            virtual void configure(const ldmx::ParameterSet& ps);

            virtual void analyze(const ldmx::Event& event);

            virtual void onFileOpen() {}

            virtual void onFileClose() {}

            virtual void onProcessStart(); 

            virtual void onProcessEnd();

        private:
            std::string caloCol_; //* Obtain the calorimeter hit information
            
            int nNotIncluded_; //* Number of hits not in the back Hcal or above minPE_ threshold

            float minPE_; //* Minimum number of PEs to not be considered noise

            int nStripsPerLayer_;
            int nLayers_;

            TH1F* h_includedhits; //* PE distribution of included hits
    };
}

#endif /* PLAYTEST_HCALLAYERANALYZER_H */

