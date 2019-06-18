/**
 * @file Trigscinthit.h
 * @brief Class that stores Stores reconstructed hit information from the HCAL
 * @author Andrew Whitbeck, Texas Tech University
 */

#ifndef EVENT_TRIGSCINTHIT_H_
#define EVENT_TRIGSCINTHIT_H_

// LDMX
#include "Event/HcalHit.h"

namespace ldmx {

    /**
     * @class Trigscinthit
     * @brief Stores reconstructed hit information from the HCAL
     *
     * @note This class represents the reconstructed hit information
     * from the trigger scintillator. 
     */
    class TrigScintHit : public HcalHit {

        public:

            /**
             * Class constructor.
             */
            TrigScintHit() {
            }

            /**
             * Class destructor.
             */
            virtual ~TrigScintHit() {
            }

            /**
             * Clear the data in the object.
             */
            void Clear(Option_t *option = "");

            /**
             * Print out the object.
             */
            void Print(Option_t *option = "") const;

            /** Decode the section associated with the hit from the ID. */
            virtual int getSection() const override{
	        return (getID() & 0x7000) >> 12;
            }

            /** Decode the strip associated with the hit from the ID. */
            virtual int getStrip() const override{
                return (getID() & 0x7F8000) >> 15;
            }

	    /** Set time of hit. */
	    void setTime(float t){
	      time_ = t;
	    }

	    /** Set beam energy fraction of hit. */
	    void setBeamEfrac(float e){
	      beamEfrac_ = e;
	    }
	  
        private:

            /** The time estimated for this hit. */
            float time_{0};
	    
	    /** The fraction of energy associated with beam electrons. */
	    float beamEfrac_{0};

            /**
             * The ROOT class definition.
             */
            ClassDef(TrigScintHit, 1);
    };

}

#endif /* EVENT_TRIGSCINTHIT_H_ */
