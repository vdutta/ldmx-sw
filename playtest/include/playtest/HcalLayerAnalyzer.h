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
#include <map> //std::map for storage tree in ClusterLog class
#include <utility> //std::pair for storage tree in ClusterLog class
#include <iterator> //std::next and std::prev for search through map
#include <cmath> //floor and ceil for stripbounds calculations

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
    
    //* type that will be used to reference hits
    typedef const ldmx::HcalHit* HitPtr;

    //* type that will be used for logging hits by layer,strip
    typedef std::map< int , HitPtr > HitLog;

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

            int layermod_; //* layer modulus to be used for key gen

            int nStrips_; //* number of strips per layer of Hcal

            int nEcalThickness_; //* thickness of Ecal in number of strips (rounded up)

            float origin_; //* origin of Ecal in strips
            float lowside_; //* low side of Ecal in strips
            float upside_; //* upper side of Ecal in strips

            TH1F* h_includedhits; //* PE distribution of included hits over all layers and events

            /**
             * Function to generate key for a given hit.
             * Relies on layermod_.
             *
             * @param hit pointer to HcalHit instance that a key is needed for
             * @return integer key value
             */
            int keygen( HitPtr hit ) const;

            /**
             * Function to search a specific range of a HitLog for a hit.
             * Returns a pair of nullptrs if failed to find an isolated hit.
             *
             * @param log HitLog to be searched
             * @param lowkey lower bound key
             * @param upkey upper bound key
             * @return HitPtr pair to isolated hit (single strip or two strip combo)
             */
            std::pair< HitPtr , HitPtr > search( const HitLog log , const int lowkey , const int upkey ) const;

            /**
             * Function to find strip bounds for input layer given seed layer.
             *
             * @param seedlayer layer number of seed
             * @param seedstrip strip number of seed
             * @param layer layer number of interest
             * @param lowstrip strip number of lower bound
             * @param upstrip strip number of upper bound
             * @return true if projected range intersects the input layer
             */
            bool stripbounds( const int seedlayer , const int seedstrip , const int layer , int &lowstrip , int &upstrip ) const;

            /**
             * Function to find seedstrip from seedlayer by finding an isolated hit.
             * Will change seedlayer if unable to find isolated hit in that layer.
             *
             * @param log HitLog to be searched
             * @param seedlayer layer number of seed
             * @param seedstrip strip number of seed
             * @return true if search for seed was successful
             */
            bool findseed( const HitLog log , int &seedlayer , int &seedstrip ) const;

    };
}

#endif /* PLAYTEST_HCALLAYERANALYZER_H */

