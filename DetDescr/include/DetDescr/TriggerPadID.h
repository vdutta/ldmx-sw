/**
 * @file TriggerPadID.h
 * @brief Class that defines a TriggerPad sensitive detector
 * @author Jeremy McCormick, SLAC National Accelerator Laboratory
 * @author Lene Kristian Bryngemark, Stanford University
 */

#ifndef DETDESCR_TRIGGERPADID_H_
#define DETDESCR_TRIGGERPADID_H_

// LDMX
#include "DetDescr/DefaultDetectorID.h"

namespace ldmx {

    /**
     * Encodes the pad of the TriggerPad based on the 'pad' field value.
     */
    enum TriggerPadNb {
        TAGGER = 0,
        UP = 1,
        DOWN = 2,
    };

    /**
     * @class TriggerPadID
     * @brief Implements sensitive detector for Trigger scintillator pad subdetector
     */
    class TriggerPadID : public DefaultDetectorID {

        public:

            TriggerPadID() {
                this->getFieldList()->push_back(new IDField("pad", 1, 12, 13));
                this->getFieldList()->push_back(new IDField("layer", 2, 14, 14));
                this->getFieldList()->push_back(new IDField("strip", 3, 15, 21));
                init();
            }

            /**
             * Get the value of the 'pad' field from the ID.
             * @return The value of the 'pad' field.
             */
            int getPad() {
                return this->getFieldValue(1);
            }

            /**
             * Get the value of the 'layer' field from the ID.
             * @return The value of the 'layer' field.
             */
            int getLayer() {
                return this->getFieldValue(2);
            }

            /**
             * Get the value of the 'strip' field from the ID.
             * @return The value of 'strip' field.
             */
            int getStrip() {
                return this->getFieldValue(3);
            }
    };
}

#endif
