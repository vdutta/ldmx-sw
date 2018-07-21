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
#include "TGraph.h" //Linearly extrapolate strip to new layer
#include "TF1.h" //Linearly extrapolate strip to new layer

//LDMX Framework
#include "Event/Event.h" //add new TConesArray to event bus
#include "Framework/EventProcessor.h" //inherit from Producer class
#include "Framework/ParameterSet.h" //get parameters from config script
#include "Event/HcalTrack.h" //store tracks generated
#include "DetDescr/HcalID.h" //for HcalSection enum

namespace ldmx {

    /**
     * @class HcalTrackProducer
     * @brief Stores HitPtrs in a std::map for easy searching and track reconstruction.
     *
     * @note Currently, alternating bar/strip orientation is not implemented in the HCAL simulation.
     *  Therefore, all of the functions in this class assume that all the layers have the same orientation.
     *  In order to switch to alternating bar orientation, replace "true ) { //" with the empty string.
     *  In vim and in command mode type ":%s/true\ )\ {\ \/\///g" when in the source file.
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
             * Add hit to log_.
             *
             * @param hit pointer to hit that will be added to log_.
             */
            void AddHit( HitPtr hit );
            
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
             * Search for next mip given layer and partial track.
             * Will add to track if found a mip hit.
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
            
            std::string hitcollname_; //* name of collection of hits
            std::string hitpassname_; //* name of pass that made hit collection

            int nlayers_; //* number of layers in detector
            int nstrips_; //* number of strips per layer
            
            int layermod_; //* layer modulus to use for hit keys
            int sectionmod_; //* section modulus to use for hit keys

            float minPE_; //* Minimum number of PEs to not be considered noise
            float maxEnergy_; //* Maximum energy of a hit to not be considered a non-mip
            
            int firstseedlayer_; //* first seed layer to try

            int conedepth_; //* depth of search cone around seed in layers
            int coneangle_; //* angular opening of cone around seed in strips across the first layer
            int minconehits_; //* minimum number of hits in cone around seed to allow for seed to be accepted

            int trackwidth_; //* width of extended track to search in number of strips
            
            int mintracklayhits_; //* minimum number of layers hit in a full track for it to be accepted
            
            int maxtrackcnt_; //* maximum number of tracks that can be found before exiting

            int seedlayer_; //* current layer number used as seed for current track
            int seedstrip_; //* current strip number used as seed for current track

            std::string hcaltracksname_; //* name of track collection to be put into event bus
            TClonesArray* hcaltracks_; //* array of HcalTracks that are found in a given event

            std::map< int , HitPtr > log_; //* map that will be used to store the hits

            std::set< int > layercheck_; //* set of layers that haven't been checked exhaustively yet 

            std::list< std::pair< int , int > > cone_; //* search cone in keys around seed
            std::list< int > layerlist_; //* list of layers to go through after partial track is begun
            std::set< int > badseeds_; //* set of seedkeys that end up not being able to start a track
            
    };

}

#endif /* HCAL_HCALTRACKPRODUCER_H */
