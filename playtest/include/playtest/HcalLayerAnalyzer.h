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
#include <utilities> //std::pair for storage tree in ClusterLog class

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
     * New type just to save on typing
     * Clusters hits togther if mip passes through two neighbor strips in same layer.
     */
    typedef std::pair< const ldmx::HcalHit* , const ldmx::HcalHit* > Cluster;
    
    /**
     * constant layer modulus for cluster keys.
     */
    const int HCAL_LAYER_MOD = 1000;

    /**
     * @class Cluster
     * @brief Grouping class to store pointers to hits and cluster key
     *
     */
    class Cluster {
        public:
            
    };

    /**
     * @class ClusterLog
     * @brief Helper class that logs the hits for a certain event in a useful way.
     *
     * @note The implementation of this class will be variable. Not currently confident in any particular
     *  method of storing the hits, so it may go through several versions.
     */
    class ClusterLog {
        public:
            
            /**
             * Default Constructor.
             */
            ClusterLog();

            /**
             * Input new element function.
             *
             * @param cluster instance of Cluster that stores a pair of hits
             */
            void Insert( Cluster cluster );

            /**
             * Find cluster in range between lowstrip and highstrip in layer.
             * Will return FIRST cluster it finds, does NOT check for multiple clusters.
             *
             * @param lowstrip lower bound on strip number
             * @param highstrip upper bound on strip number
             * @param layer layer to search
             * @return Cluster that contains findings (two nullptrs if nothing)
             */
            Cluster Find( const int layer , const int lowstrip , const int highstrip ) const;

            /**
             * Removes a vector of keys from log.
             *
             * @param track vector of keys
             */
            void RemoveTrack( std::vector< int > track );

        private:
            
            //* storage tree for logging hits and indices
            std::map< int , Cluster > log_;
            
            /**
             * Key generator for a specific layer and strip pair.
             *
             * @param layer integer layer number
             * @param strip integer strip number
             * @return integer key value
             */
            int KeyGen( const int layer , const int strip ) const;
            
            /**
             * Key generator for new cluster.
             *
             * @param cluster Cluster instance for which key will be generated.
             * @return integer key value
             */
            int KeyGen( Cluster cluster ) const;

             
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

