/**
 * @file HcalTrack.h
 * @brief Class implementation of a track through the Hcal
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
            HcalTrack() : nhits_(0) { }
            
            /**
             * Copy Constructor
             */
            HcalTrack( const HcalTrack& tocopy ) : hits_(tocopy.hits_),nhits_(tocopy.nhits_)  {
                
            }

            /**
             * Clear the track
             */
            void Clear() {
                TObject::Clear();

                hits_.Clear();
                nhits_ = 0;

                return;
            }

            /**
             * Add a hit to the track.
             */
            void addHit( HitPtr hit ) {
                hits_.Add( hit );
                nhits_++;
                return;
            }

            /**
             * Get number of hits in track
             */
            int getNHits() const {
                return nhits_;
            }

            /**
             * Get hit at a certain index in track.
             */
            HitPtr getHit( int i ) const {
                return (HitPtr)(hits_.At(i));
            }

        private:
            
            TRefArray hits_; //* references to hits in the track
            int nhits_; //* number of hits in the track
            
            /**
             * ROOT Class Definition
             */
            ClassDef( HcalTrack , 1 );

    };

}

#endif /* PLAYTEST_HCALTRACK_H */
