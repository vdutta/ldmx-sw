/**
 * @file HcalTrackProducer.h
 * @brief Class implementation of hit log to search and reconstruct.
 * @author Tom Eichlersmith, University of Minnesota
 */

#ifndef HCAL_HCALTRACKPRODUCER_H
#define HCAL_HCALTRACKPROCUCER_H

//Standard Libraries
#include <iostream> //Checks to std::cout for development purposes
#include <vector> //Vector of hits per layer
#include <map> //std::map for storage tree
#include <list> //std::queue for search cone and layer list
#include <set> //std::set for layers that haven't been checked yet
#include <utility> //std::pair for storage tree
#include <iterator> //std::next and std::prev for search through map
#include <cmath> //floor and ceil for stripbounds calculations
#include <algorithm> //std::sort for layerlist

//ROOT
#include "TClonesArray.h" //Add new array of tracks to event bus

//LDMX Framework
#include "Event/Event.h" //add new TConesArray to event bus
#include "Framework/EventProcessor.h" //inherit from Producer class
#include "Framework/ParameterSet.h" //get parameters from config script
#include "Event/HcalTrack.h" //store tracks generated
#include "DetDescr/HcalID.h" //for HcalSection enum
#include "Hcal/MipHit.h" //for MipHit class

namespace ldmx {

    /**
     * @class HcalTrackProducer
     * @brief Stores HitPtrs in a std::map for easy searching and track reconstruction.
     *
     * @note Currently, alternating bar/strip orientation is not implemented in the HCAL simulation.
     *  Therefore, all of the functions in this class assume that all the layers in the BACK HCAL have the same orientation.
     *
     * @note The idea to incorporate the different sections and the alternating orientations in the Back Hcal is to resolve
     *  tracks in each section separately and then see if any can be joined. The horizontal layers and the vertical layers
     *  in the Back Hcal will be treated as separte sections.
     */
    class HcalTrackProducer : public Producer {
        public:
            /**
             * Default Constructor.
             */
            HcalTrackProducer( const std::string& name , Process& process ) : Producer( name , process ) { }
            
            /** Destructor */
            ~HcalTrackProducer();

            virtual void configure( const ParameterSet& ps );

            virtual void produce( Event& event );

            virtual void onFileOpen() { }

            virtual void onFileClose() { }

            virtual void onProcessStart() { }

            virtual void onProcessEnd() { }

        private:
            
            /**
             * Attempt to find reconstruct tracks in the input section.
             *
             * @param section HcalSection to search
             * @return true if at least one track was found
             */
            bool TrackSearch( HcalSection section );
   
            /**
             * Function to generate key from section,layer,hit information.
             *
             * @param section HcalSection enum specify the section
             * @param layer layer number
             * @param strip strip number
             * @return integer key value
             */
            int KeyGen( const int section , const int layer , const int strip ) const;

            /**
             * Function to generate key for a given hit.
             * Relies on layer and section moduli.
             *
             * @param hit pointer to HcalHit instance that a key is needed for
             * @return integer key value
             */
            int KeyGen( HitPtr hit ) const;

            /**
             * Function to generate key for a given mip.
             *
             * @param mip pointer to MipHit instance that a key is needed for
             * @return integer key value
             */
            int KeyGen( MipHitPtr mip ) const;

            /**
             * Function to produce end points for a given section.
             *
             * @param section HcalSection to search
             * @param endpts std::pair of MipHitPtrs that will be the end points
             * @return true if endponts were found
             */
            bool findEndPts( HcalSection section , std::pair< MipHitPtr > &endpts );
            
            /**
             * Connect Given End Points into a possible track.
             *
             * @param section HcalSection to search through
             * @param endpts std::pair of MipHitPtrs that are the end points to use
             * @param possibletrack std::vector of MipHitPtrs that contain the possible track
             * @return true if possibletrack is acceptable as a track
             */
            bool connectTrack( HcalSection section , const std::pair< MipHitPtr > &endpts , std::vector< MipHitPtr > possibletrack );

            /**
             * Check if plausible track is acceptable.
             *
             * @param track vector of MipHitPtrs to check
             * @return true if acceptable
             */
            bool isAcceptableTrack( const std::vector< MipHitPtr > &track ) const;

            /**
             * Function to determine whether a group of hits can be considered a mip.
             * checks on number of individual hits, and total energy of group.
             * Assumes group is isolated (defined as all hits in the group are adjacent strips
             *  and the strips on either side of the group are empty).
             *
             * @param group vector of HitPtr that contains the hits in a group
             * @return true if considered a mip
             */
            bool isMIP( const std::vector< HitPtr > &group ) const;

            std::string hitCollName_; //* name of collection of hits
            std::string hitPassName_; //* name of pass that made hit collection

            int nLayersBack_; //* number of layers in Back HCAL
            int nStripsBack_; //* number of strips per layer in Back HCAL
            int nLayersTopBot_; //* number of layers in Top and Bottom HCAL
            int nStripsTopBot_; //* number of strips per layer in Top and Bottom HCAL
            int nLayersLeftRight_; //* number of layers in Left and Right HCAL
            int nStripsLeftRight_; //* number of strips per layer in Left and Right HCAL
            
            int layerMod_; //* layer modulus to use for hit keys
            int sectionMod_; //* sectio modulus to use for hit keys

            float minPE_; //* Minimum number of PEs to not be considered noise
            float maxEnergy_; //* Maximum energy of a hit to not be considered a non-mip
            
            int trackWidth_; //* width of extended track to search in number of strips
            
            int minTrackLayHits_; //* minimum number of layers hit in a full track for it to be accepted
            
            int maxTrackCnt_; //* maximum number of tracks that can be found before exiting

            std::string hcalTracksName_; //* name of track collection to be put into event bus
            TClonesArray* hcalTracks_; //* array of HcalTracks that are found in a given event

            std::map< int , HitPtr > nonoiseLog_; //* map that will be used to store the hits that aren't noise (above minPE)

            std::map< int , MipHitPtr > mipLog_ //* map that stores the pre-processed mips

            std::map< int , HcalTrackPtr > trackLog_; //* map that stores the tracks found in each section

            std::set< int > badEndPts_; //* set of mip keys that do not produce good tracks

    };

}

#endif /* HCAL_HCALTRACKPRODUCER_H */
