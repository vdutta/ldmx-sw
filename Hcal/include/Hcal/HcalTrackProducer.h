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
#include "Hcal/MIPHit.h" //for MIPHit class

namespace ldmx {

    /**
     * @class HcalTrackProducer
     * @brief Stores HitPtrs in a std::map for easy searching and track reconstruction.
     *
     * @note Currently, alternating bar/strip orientation is not implemented in the HCAL simulation.
     *  Therefore, all of the functions in this class assume that all the layers in the BACK HCAL have the same orientation.
     */
    class HcalTrackProducer : public Producer {
        public:
            /**
             * Default Constructor.
             */
            HcalTrackProducer( const std::string& name , Process& process ) : Producer( name , process ) {}
            
            virtual void configure( const ParameterSet& ps );

            virtual void produce( Event& event );

            virtual void onFileOpen() { }

            virtual void onFileClose() { }

            virtual void onProcessStart() { }

            virtual void onProcessEnd() { }

        private:
            
            /**
             * Remove list of hits for a track.
             *
             * @param track HcalTrack to be removed.
             */
            void RemoveTrack( const HcalTrack *track );

            /**
             * Attempt to reconstruct a track from a seed layer.
             *
             * @param track plausible track - should be empty
             * @return true if a track was found
             */
            bool TrackSearch( HcalTrack *track );
   
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
             * Corrects for strip numbers outside of real range.
             * Negative strip numbers are set to zero and numbers greater
             *  than the number of strips are set to nstrips_.
             *
             * @param strip the strip number to correct
             */
            void CorrectStrip( int &strip) const;
            
            /**
             * Finds a seed strip given a seed layer.
             * Recursively looks through more different seed layers if doesn't find one in input layer.
             * Will return false if it exhausts its layer options for a seed.
             *
             * @return true if found a seed (seedlayer_ and seedstrip_ are its position)
             */
            bool FindSeed();
            
            /**
             * Constructs search cone around seed and list of layers that aren't in cone or seed.
             */
            void SetSearchCone();
            
            /**
             * Begins partial track by searching through cone around seed.
             *
             * @param track HcalTrack that stores beginning of track
             * @return true if successfully started track
             */
            bool BeginPartialTrack( HcalTrack *track );

            /**
             * Search through layers for mips to add to partial track.
             * Assumes track has AT LEAST two hits in it.
             *
             * @param track HcalTrack to be extended
             * @return true if acceptable track was created
             */
            bool ExtendTrack( HcalTrack *track );

            /**
             * Check if plausible track is acceptable.
             *
             * @param track HcalTrack to check
             * @return true if acceptable
             */
            bool isAcceptableTrack( const HcalTrack *track ) const;

            /**
             * Function to search a specific range of a log for a hit.
             * Will add hit(s) to track that it considers to be the preferred mip.
             * If more than one mip is found, the one closest to prefkey is the one added.
             *
             * @param lowkey lower bound key
             * @param upkey upper bound key
             * @param track partial track to add hit to if found
             * @param prefkey hit key that is given preference if multiple hits are found
             * @return true if successfully found a hit in the given key range 
             */
            bool SearchByKey( const int lowkey , const int upkey , HcalTrack *track , const float prefkey = -1.0 );

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

            std::map< HcalSection , std::map< int , MIPHit > > mipLog_ //* map that stores the pre-processed mips

    };

}

#endif /* HCAL_HCALTRACKPRODUCER_H */
