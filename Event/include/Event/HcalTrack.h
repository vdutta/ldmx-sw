/**
 * @file HcalTrack.h
 * @brief Class declaration of a track through the Hcal
 * @author Tom Eichlersmith, University of Minnesota
 */

#ifndef PLAYTEST_HCALTRACK_H
#define PLAYTEST_HCALTRACK_H

//Standard Libraries
#include <iostream> //Checks to std::cout for development purposes

//ROOT
#include "TRefArray.h" //store pointers to hits
#include "TObject.h" //inherit from TObject

//LDMX Framework
#include "Event/HcalHit.h" //Get Hcal Specific information from hit

namespace ldmx {

    //* type that will be used to reference hits
    typedef ldmx::HcalHit* HitPtr;

    /**
     * @class HcalTrack
     * @brief Storage object for a track through the Hcal.
     */
    class HcalTrack : public TObject {
        public:
            /**
             * Default Constructor
             */
            HcalTrack() 
                : TObject(), hits_(new TRefArray()), nlayhits_(0),
                  seedlayer_(0), seedstrip_(0) { }
            
            /**
             * Destructor
             * Clears the TRefArray and Frees the memory
             */
            ~HcalTrack(); 
            
            /**
             * Assignment Operator
             * Explicity copies over every member variable including the information in the TRefArray (using its copy constructor).
             */
            HcalTrack& operator= ( const HcalTrack &track );  

            /**
             * Clear the track
             * Sets member variables to zero and empties the TRefArray (does NONT free the memory).
             */
            void Clear(Option_t *opt = ""); 

            /**
             * Add a hit to the track.
             */
            void addHit( HitPtr hit );

            /**
             * Increment the number of layers hit by one.
             */
            void incLayHit();

            /**
             * Set seed information
             */
            void setSeed(int seedlayer , int seedstrip);

            /**
             * Add a group of hits to the track.
             */
            void addGroup( const std::vector<HitPtr> group ); 

            /**
             * Get number of hits in track
             */
            int getNHits() const; 

            /**
             * Get number of layers hit in track
             */
            int getNLayHits() const; 
            
            /**
             * Get seed layer information
             */
            int getSeedLayer() const;
            
            /**
             * Get seed strip information
             */
            int getSeedStrip() const; 

            /**
             * Get hit at a certain index in track.
             */
            HitPtr getHit( int i ) const; 

        private:
            
            TRefArray *hits_; //* references to hits in the track
            int nlayhits_; //* number of layers hit in the track
            
            int seedlayer_; //* layer of seed for this track
            int seedstrip_; //* strip of seed for this track

            /**
             * ROOT Class Definition
             */
            ClassDef( HcalTrack , 17 );

    };

}

#endif /* PLAYTEST_HCALTRACK_H */
