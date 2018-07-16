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
            HcalTrack() : TObject(), nhits_(0), nlayhits_(0), seedlayer_(1000), seedstrip_(1000) { }
            
            /**
             * Destructor
             */
            ~HcalTrack() {
                Clear();
            }
            
            /**
             * Clear the track
             */
            void Clear() {

                TObject::Clear();

                hits_.Clear();
                nhits_ = 0;
                nlayhits_ = 0;

                seedlayer_ = 1000;
                seedstrip_ = 1000;

                return;
            }

            /**
             * Add a hit to the track.
             */
            void addHit( HitPtr hit ) {
                hits_.AddAtFree( hit );
                nhits_++;
                return;
            }

            /**
             * Increment the number of layers hit by one.
             */
            void incLayHit() {
                nlayhits_++;
                return;
            }

            /**
             * Set seed information
             */
            void setSeed(int seedlayer , int seedstrip) {
                seedlayer_ = seedlayer;
                seedstrip_ = seedstrip;
                return;
            }

            /**
             * Add a group of hits to the track.
             */
            void addGroup( const std::vector<HitPtr> group ) {
                for( auto it = group.begin(); it != group.end(); ++it ) {
                    addHit( *it );
                }
                return;
            }

            /**
             * Get number of hits in track
             */
            int getNHits() const {
                return nhits_;
            }

            /**
             * Get number of layers hit in track
             */
            int getNLayHits() const {
                return nlayhits_;
            }
            
            /**
             * Get seed information
             */
            int getSeedLayer() const {
                return seedlayer_;
            }

            int getSeedStrip() const {
                return seedstrip_;
            }

            /**
             * Get hit at a certain index in track.
             */
            HitPtr getHit( int i ) const {
                return (HitPtr)(hits_.At(i));
            }

            /**
             * Get full TRefArray.
             */
            const TRefArray &getTrack() {
                return hits_;
            }

        private:
            
            TRefArray hits_; //* references to hits in the track
            int nhits_; //* number of hits in the track
            int nlayhits_; //* number of layers hit in the track
            
            int seedlayer_; //* layer of seed for this track
            int seedstrip_; //* strip of seed for this track

            /**
             * ROOT Class Definition
             */
            ClassDef( HcalTrack , 1 );

    };

}

#endif /* PLAYTEST_HCALTRACK_H */
