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
#include <set> //set of bad seed ids
#include <iterator> //std::prev
#include <cmath> //sqrt, abs, asin
#include <iostream> //cerr and cout
#include <ctime> //timing produce function

//LDMX Framework
#include "Event/Event.h"
#include "Framework/EventProcessor.h" //Needed to declare processor
#include "Framework/ParameterSet.h" // Needed to import parameters from configuration file
#include "Event/EventConstants.h" //for HcalMipTrack string
#include "Event/HcalMipTrack.h" //mip track container
#include "Hcal/MipCluster.h" //cluster object
#include "Tools/HitBox.h" //real space representation of MipCluster

namespace ldmx {
    
    /**
     * @class HcalMipTrackProducer
     * @brief Producer that reconstructs MIP tracks through the hcal.
     *
     */
    class HcalMipTrackProducer : public ldmx::Producer {
        public:

            HcalMipTrackProducer(const std::string& name, ldmx::Process& process) : ldmx::Producer(name, process) { }

            virtual void configure(const ldmx::ParameterSet& ps);

            virtual void produce(ldmx::Event& event);

            virtual void onFileOpen() { }

            virtual void onFileClose() { }

            virtual void onProcessStart() { } 
            
            /**
             * Prints performance trackers
             */
            virtual void onProcessEnd(); 

        private:
            
            /**
             * Check if a single HcalHit is not noise.
             *
             * @param hit HcalHit to check
             * @return true if not noise
             */
            bool isNotNoise( HcalHit* hit ) const;

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
             * Finds a seed to construct a track from.
             * Choose a good seed cluster that has the lowest z coordinate.
             * Unless using the median is specified.
             *
             * @return true if successfully found a seed
             */
            bool findSeed( const bool useMedian );

            /**
             * Finds best track out of the clusters that are left.
             * 
             * @param track_mipids vector of cluster ids in track that is found
             * @return true if track was found
             */
            bool buildTrack( std::vector< unsigned int > &track_mipids );

            /**
             * Determine if line (defined by origin and direction) hits the axis-oriented box
             *  defined by minimum and maximum corners.
             * Assumes all inputs are formatted correctly.
             * This algorithm is the Fast Ray-Box Intersection Algorithm developed by Andrew Woo
             *  from "Graphics Gems", Academic Press, 1990.
             * This algorithm returns false if ray points away from box
             *  (even if the extended line would intersect it),
             *  so it is suggested to check both directions from origin.
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
             * Determine if track from the list of mip ids is acceptable.
             *
             * @return true if acceptable
             */
            bool isAcceptableTrack( const std::vector< unsigned int > &track_mipids ) const;

            /**
             * Checks if the two tracks should be merged together or not.
             * If the two tracks direction vectors have an angle between them that is 
             *  smaller than the input parameter, the two tracks are merged.
             *
             * Uses cross product to calculate angle between tracks (theta):
             *  |first x second| = |first||second|sin(theta)
             *  theta = asin( | first x second | / |first||second| )
             *
             * Compares abs(theta) to input parameter.
             *
             * @return true if should be merged
             */
            bool shouldMergeTracks( HcalMipTrack *first , HcalMipTrack *second ) const;

            /**
             * Calculates closest distance between a line and a point
             * Assumes inputs are correctly formatted and the direction vector of the line
             *  is normalized (unit length);
             */
            double distPt2Line( const std::vector<double> &start, const std::vector<double> &dir,
                                const std::vector<double> &point ) const;

            /** Name of collection of HcalHits */
            std::string hcalHitCollName_;

            /** Name of pass to get HcalHits collection from */
            std::string hcalHitPassName_;

            /** Container for reconstructed tracks */
            TClonesArray *hcalMipTracks_;

            /** Collection name for event bus */
            std::string hcalMipTracksCollName_;
            
            /** Maximum number of tracks to be allowed, prevents infinite looping */
            int maxTrackCount_;

            /** Minimum number of PE to be considered not noise */
            double minPE_;

            /** Maximum energy of a cluster to be considered a mip */
            double maxEnergy_;

            /** Fraction of total number of clusters required to find a seed */
            double fracTotalClusters_;

            /** Fraction of clusters currently in log to consider track acceptable */
            double fracClustersLeft_;

            /** Maximum difference in slope angles of two tracks to merge */
            double maxSlopeAngleDiff_;

            /** Maximum distance between two tracks to merge */
            double maxTrackDist_;

            /** Log of HcalHits, sorted by section,layer,strip information */
            std::map< unsigned int , HcalHit* > hcalHitLog_;

            /** Log of clusters still being considered for track */
            std::map< unsigned int , MipCluster > clusterLog_;

            /** Minimum number of clusters to allow to search for more tracks */
            int minNumClusters_;

            /** Seed cluster box */
            HitBox seedBox_;
            
            /** Seed Mip ID */
            unsigned int seedID_;

            // PERFORMANCE TRACKERS

            /** Record the number of events that have a certain number of tracks */
            std::map< int , int > numTracksPerEvent_;

            /** Record the number of clusters in each track per event (separated by numTracksPerEvent) */
            std::map< int , double > meanClustersPerTrack_;

            /** Record the mean time it takes to produce the tracks */
            double meanTime_produce_;
            
            /** Number of accesses to both logs in current event */
            unsigned long int numTouchLogs_;

            /** Number of accesses to both logs */
            double meanNumTouchLogs_;

            /** Number of clusters not included in track */
            double meanClustersIgnored_;
    };
}

#endif /* HCAL_HCALMIPTRACKPRODUCER_H */
