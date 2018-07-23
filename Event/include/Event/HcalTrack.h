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
#include "TF1.h" //fit function
#include "TGraph.h" //store points for fit

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
            HcalTrack(); 
            
            /**
             * Destructor
             * Clears the TRefArray and Frees the memory
             * Deltes the TGraphs and TF1s
             */
            ~HcalTrack(); 
            
            /**
             * Assignment Operator
             * Explicity copies over every member variable including the information in the TRefArray (using its copy constructor).
             */
            HcalTrack& operator= ( const HcalTrack &track );  

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
             * Evaluate the fit at a given layer.
             */
            float evalFit( const int layer );

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
            
            /**
             * Check to see if HcalTrack is empty.
             * Checks size of TRefArray.
             */
            bool isEmpty() const;

            /**
             * Check if HcalTrack is broken.
             * There is a hit that is a nullptr.
             */
            bool isBroken() const;

        private:
            
            TRefArray *hits_; //* references to hits in the track
            int nlayhits_; //* number of layers hit in the track
            
            int seedlayer_; //* layer of seed for this track
            int seedstrip_; //* strip of seed for this track

            TGraph* oddgr_; //* most up-to-date data points for odd layers
            TGraph* evengr_; //* most up-to-date data points for even layers
            
            TF1* fitres_; //* resulting function when fit is performed on TGraphs.

            /**
             * ROOT Class Definition
             */
            ClassDef( HcalTrack , 23 );

    };

}

#endif /* PLAYTEST_HCALTRACK_H */
