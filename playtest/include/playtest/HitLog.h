/**
 * @file HitLog.h
 * @brief Class implementation of hit log to search and reconstruct.
 * @author Tom Eichlersmith, University of Minnesota
 */

#ifndef PLAYTEST_HITLOG_H
#define PLAYTEST_HITLOG_H

//Standard Libraries
#include <iostream> //Checks to std::cout for development purposes
#include <vector> //Vector of hits per layer
#include <map> //std::map for storage tree
#include <stack> //std::stack for search cone
#include <utility> //std::pair for storage tree
#include <iterator> //std::next and std::prev for search through map
#include <cmath> //floor and ceil for stripbounds calculations

//LDMX Framework
#include "Event/HcalHit.h" //Get Hcal Specific information from hit


namespace ldmx {

    //* type that will be used to reference hits
    typedef const ldmx::HcalHit* HitPtr;

    /**
     * @class HitLog
     * @brief Stores HitPtrs in a std::map for easy searching and track reconstruction.
     */
    class HitLog {
        public:
            /**
             * Default Constructor.
             */
            HitLog();

            /**
             * Preferred Constructor.
             */
            HitLog( const float minPE , const int conedepth , const int coneangle , const float origin ,
                    const float lowside , const float upside , const int nstrips , const int nlayers );

            /**
             * Add hit to log
             */
            void AddHit( HitPtr hit );


        private:
            
            std::map< int , HitPtr > log_; //* map that will be used to store the hits

            std::vector< bool > layercheck_; //* vector of layers to see if searched or not

            std::stack< std::pair< int , int > > cone_; //* search cone in keys around seed
            
            int nlayers_; //* number of layers in detector
            int nstrips_; //* number of strips per layer

            float minPE_; //* Minimum number of PEs to not be considered noise

            int conedepth_; //* depth of search cone around seed in layers

            int coneangle_; //* angular opening of cone around seed in strips across the first layer

            float origin_; //* center of Ecal in strips
            float lowside_; //* low side of Ecal in strips
            float upside_; //* upper side of Ecal in strips

            /**
             * Function to generate key for a given hit.
             * Relies on layermod_.
             *
             * @param hit pointer to HcalHit instance that a key is needed for
             * @return integer key value
             */
            int KeyGen( HitPtr hit ) const;

            /**
             * Constructs search cone around seed
             *
             */
            void SetSearchCone( const int seedlayer , const int seedstrip );

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
             * This function assumes that the possible track originates from the Ecal.
             * Within a few layers of the seed, this function will do a good job shrinking the parameter space.
             * Farther away from the seed (in layers), this function will only work if the track did originate
             *  from the Ecal.
             *
             * @param seedlayer layer number of seed
             * @param seedstrip strip number of seed
             * @param layer layer number of interest
             * @param lowstrip strip number of lower bound
             * @param upstrip strip number of upper bound
             * @return true if projected range intersects the input layer
             */
            bool stripbounds( const int seedlayer , const int seedstrip , const int layer , int &lowstrip , int &upstrip ) const;

    };

}

#endif /* PLAYTEST_HITLOG_H */
