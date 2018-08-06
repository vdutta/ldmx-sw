/**
 * @file HcalMipTrackProducer.h
 * @brief Header file for HcalMipTrackProducer class
 * @author Tom Eichlersmith, University of Minnesota
 */

#ifndef HCAL_HCALMIPTRACKPRODUCER_H
#define HCAL_HCALMIPTRACKPRODUCER_H

//STL
#include <string> //names of collections
#include <map> //cluster and hit logs
#include <iterator> //std::prev
#include <cmath> //sqrt and abs

//LDMX Framework
#include "Event/Event.h"
#include "Framework/EventProcessor.h" //Needed to declare processor
#include "Framework/ParameterSet.h" // Needed to import parameters from configuration file
#include "Event/EventConstants.h" //for HcalMipTrack string
#include "Event/HcalMipTrack.h" //mip track container
#include "Hcal/HcalDetectorGeometry.h" //calculating real space coordinates

namespace ldmx {
    
    /**
     * @class HcalMipTrackProducer
     * @brief Producer that reconstructs MIP tracks through the hcal.
     */
    class HcalMipTrackProducer : public ldmx::Producer {
        public:

            HcalMipTrackProducer(const std::string& name, ldmx::Process& process) : ldmx::Producer(name, process) { }

            virtual void configure(const ldmx::ParameterSet& ps);

            virtual void produce(ldmx::Event& event);

            virtual void onFileOpen() { }

            virtual void onFileClose() { }

            virtual void onProcessStart() { } 

            virtual void onProcessEnd() { }

        private:
            
            /**
             * Check if a single HcalHit is not noise.
             *
             * @param hit HcalHit to check
             * @return true if not noise
             */
            bool isNotNoise( const HcalHit* hit ) const;

            /**
             * Check if cluster is a mip.
             *
             * @param cluster MipCluster to check
             * @bool true if a mip
             */
            bool isMip( const MipCluster* cluster ) const;

            /**
             * Cluster HcalHits into MipClusters.
             *
             * @note Currently, this only clusters hits in the same section,layer. One can update this
             *  to cluster across layers.
             */
            void clusterHits();

            /**
             * Determine if line (defined by origin and direction) hits the box defined by minimum and maximum corners.
             * Assumes all inputs are formatted correctly.
             * This algorithm is the Fast Ray-Box Intersection Algorithm developed by Andrew Woo
             *  from "Graphics Gems", Academic Press, 1990.
             * This algorithm returns false if ray points away from box (even if the extended line would intersect it),
             *  so it is suggested to project the ray back to a plane one one side of all boxes and then use an origin there.
             *
             * @param origin originating point of ray
             * @param dir direction vector of ray
             * @param minBox minimum corner of box
             * @param maxBox maximum corner of box
             * @return true if ray hits box
             */
            bool rayHitBox( const std::vector<double> origin , const std::vector<double> dir , 
                             const std::vector<double> minBox , const std::vector<double> maxBox ) const;
            
            /**
             * Compares two mip tracks.
             *
             * @param track1 considered on the "less than" side
             * @param track2 considered on the "greather than" side
             * @return true if track1 is "worse" than track2
             */
            bool compMipTracks( const HcalMipTrack &track1 , const HcalMipTrack &track2 ) const;

            /**
             * Finds best track out of the clusters that are left.
             * 
             * @param track_mipids vector of cluster ids in track that is found
             * @return true if track was found
             */
            bool buildTrack( std::vector< unsigned int > &track_mipids );

            /** Name of collection of HcalHits */
            std::string hcalHitCollName_;

            /** Name of pass to get HcalHits collection from */
            std::string hcalHitPassName_;

            /** Container for reconstructed tracks */
            TClonesArray *hcalMipTracks_;

            /** Collection name for event bus */
            std::string hcalMipTracksCollName_;

            /** Minimum number of PE to be considered not noise */
            double minPE_;

            /** Maximum energy of a cluster to be considered a mip */
            double maxEnergy_;

            /** Minimum number of clusters in track to be considered for ranking */
            int minNumClusters_;
            
            /** Geometry class instance to calculate transformation between detector id and real space */
            static HcalDetectorGeometry hdg_;

            /** Log of HcalHits, sorted by section,layer,strip information */
            std::map< unsigned int , HcalHit* > hcalHitLog_;

            /** Log of clusters still being considered for track */
            std::map< unsigned int , MipCluster > clusterLog_;

            
    };
}

#endif /* HCAL_HCALMIPTRACKPRODUCER_H */
