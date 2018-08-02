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
             * Add a MipCluster to the track.
             * This function takes the clusters points and errors for fitting,
             *  and puts the hit in the TRefArray.
             */
            void addCluster( const MipCluster &cluster );
            
            /**
             * Get number of hits in track
             */
            int getNHits() const; 

            /**
             * Get hit at a certain index in track.
             */
            const HcalHit* getHit( int i ) const; 

            /**
             * Fit the graphs linearly and assign their values
             *  to the corresponding inputs.
             *
             * @param z z coordinate to evaluate at
             * @param x x coordinate to find
             * @param y y coordinate to find
             */
            void evalFit( const float z , float &x , float &y );
            
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
            
            /** references to hits in the track */
            TRefArray *hcalHits_; 
            
            /** Graph relating z and x coordinates */
            TGraphErrors zxGr_;

            /** Graph relating z and y coordinates */
            TGraphErrors zyGr_;
            
            /**
             * ROOT Class Definition
             */
            ClassDef( HcalMipTrack , 1 );

    };

}

#endif /* PLAYTEST_HCALMIPTRACK_H */
