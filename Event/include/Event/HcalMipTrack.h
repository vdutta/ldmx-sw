/**
 * @file HcalMipTrack.h
 * @brief Class declaration of a track through the Hcal
 * @author Tom Eichlersmith, University of Minnesota
 */

#ifndef PLAYTEST_HCALMIPTRACK_H
#define PLAYTEST_HCALMIPTRACK_H

//Standard Libraries
#include <iostream> //Checks to std::cout for development purposes

//ROOT
#include "TRefArray.h" //store pointers to hits
#include "TObject.h" //inherit from TObject
#include "TF1.h" //fit function
#include "TGraph.h" //store points for fit

//LDMX Framework
#include "Event/HcalHit.h" //Get Hcal Specific information from hit

namespace ldmx {

    /**
     * @class HcalMipTrack
     * @brief Storage object for a track through the Hcal.
     */
    class HcalMipTrack : public TObject {
        public:
            /**
             * Default Constructor
             */
            HcalMipTrack(); 
            
            /**
             * Destructor
             * Clears the TRefArray and Frees the memory
             * Deltes the TGraphs
             */
            ~HcalMipTrack(); 
            
            /**
             * Assignment Operator
             * Explicity copies over every member variable including the information
             *  in the TRefArray (using its copy constructor).
             */
            HcalMipTrack& operator= ( const HcalMipTrack &track );  

            /**
             * Clear the track
             * Sets member variables to zero and empties the TRefArray (does NOT free the memory).
             * Deletes the TGraphs and creates new (empty) ones.
             */
            void Clear(Option_t *opt = ""); 

            /**
             * Add a hit to the track.
             * This function also adds the layer,strip information to the TGraphs for fitting
             */
            void addHit( HcalHit* hit );
            
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
            HcalHit* getHit( int i ) const; 
            
            /**
             * Check to see if HcalMipTrack is empty.
             * Checks size of TRefArray.
             */
            bool isEmpty() const;

            /**
             * Check if HcalMipTrack is broken.
             * There is a hit that is a nullptr.
             */
            bool isBroken() const;

        private:
            
            TRefArray *hcalHits_; //* references to hits in the track
            
            /**
             * ROOT Class Definition
             */
            ClassDef( HcalMipTrack , 1 );

    };

}

#endif /* PLAYTEST_HCALMIPTRACK_H */
