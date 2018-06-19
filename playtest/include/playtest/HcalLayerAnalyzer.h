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
#include <map> //std::map for storage tree in HitLog class

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
     * @class HitLog
     * @brief Helper class that logs the hits for a certain event in a useful way.
     *
     * @note The implementation of this class will be variable. Not currently confident in any particular
     *  method of storing the hits, so it may go through several versions.
     */
    class HitLog {
        public:
            
            /**
             * Default Constructor.
             *
             * Sets layermod_ to 1000 and sets maps comparison function to dictorder.
             */
            HitLog();

            /**
             * Preferred Constructor.
             *
             * Sets layermod_ to input and sets comparison function to dictorder.
             */
            HitLog(const int layermod);

        private:
            
            //* storage tree for logging hits and indices
            std::vector< int , const ldmx::HcalHit* > log_;

            //* integer modulus for layer, should be larger than the number of layers
            const int layermod_;

            /**
             * Function to determine key for a hit.
             *
             * @param hit hit that key needs to be generated for
             * @return integer key
             */
            int keygen( const ldmx::HcalHit* ) const;

            /**
             * Dictionary order to determine ordering in std::map.
             *
             * @param lhs left hand side of inequality
             * @param rhs right hand side of inequality
             * @return true if lhs goes before rhs
             */
            bool dictorder( const int lhs , const int rhs ) const;
    };

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

            int nStripsPerLayer_; //* number of strips in each layer of back hcal
            int nLayers_; //* number of layers in back hcal

            TH1F* h_includedhits; //* PE distribution of included hits
    };
}

#endif /* PLAYTEST_HCALLAYERANALYZER_H */

