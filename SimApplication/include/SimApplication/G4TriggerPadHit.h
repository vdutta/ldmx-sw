/**
 * @file G4TriggerPadHit.h
 * @brief Class defining a triggerpad hit which is used to create output SimTriggerPadHit collection
 * @author Jeremy McCormick, SLAC National Accelerator Laboratory
 * @author Lene Kristian Bryngemark, Stanford University
 */

#ifndef SIMAPPLICATION_G4SIMTRIGGERPADHIT_H_
#define SIMAPPLICATION_G4SIMTRIGGERPADHIT_H_

// Geant4
#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"

// LDMX
#include "Event/SimTriggerPadHit.h"

// STL
#include <ostream>
#include <vector>

using std::vector;

namespace ldmx {

    /**
     * @class G4TriggerPadHit
     * @brief Track hit which is used to create output SimTriggerPadHit collection
     *
     * @note
     * One of these is created for every step in a TriggerPadSD.  These are basically
     * copied verbatim into SimTriggerPadHit objects by the RootPersistencyManager.
     */
    class G4TriggerPadHit : public G4VHit {

        public:

            /**
             * Class constructor.
             */
            G4TriggerPadHit() {
            }

            /**
             * Class destructor.
             */
            virtual ~G4TriggerPadHit() {
            }

            /**
             * Draw the hit in the Geant4 runtime.
             */
            void Draw();

            /**
             * Print out hit information.
             */
            void Print();

            /**
             * Print hit information to an output stream.
             * @param os The output stream.
             * @return The same output stream.
             */
            std::ostream& print(std::ostream& os);

            /**
             * Create a new hit object.
             * @param s The size of the hit object.
             */
            inline void *operator new(size_t s);

            /**
             * Delete a hit object.
             * @param aHit The hit to delete.
             */
            inline void operator delete(void *aHit);

            /**
             * Set the track ID.
             * @param trackID  The track ID.
             */
            void setTrackID(int trackID) {
                this->trackID_ = trackID;
            }

            /**
             * Get the track ID.
             * @return The track ID.
             */
            int getTrackID() {
                return this->trackID_;
            }

            /**
             * Get the detector ID.
             * @return The detector ID.
             */
            int getID() {
                return id_;
            }

            /**
             * Set the detector ID.
             * @param id The detector ID.
             */
            void setID(int id) {
                this->id_ = id;
            }

            /**
             * Get the PDG ID.
             * @return The detector ID.
             */
            int getPdgID() {
                return pdgid_;
            }

            /**
             * Set the PDG ID.
             * @param id The detector ID.
             */
            void setPdgID(int pdgid) {
                this->pdgid_ = pdgid;
            }            

            /**
             * Get the layer ID.
             * @return The layer ID.
             */
            int getLayerID() {
                return layerID_;
            }

            /**
             * Set the layer ID.
             * @param layerID The layer ID.
             */
            void setLayerID(int layerID) {
                this->layerID_ = layerID;
            }

            /**
             * Get the strip ID.
             * @return The strip ID.
             */
            int getStripID() {
                return stripID_;
            }

            /**
             * Set the strip ID.
             * @param stripID The strip ID.
             */
            void setStripID(int stripID) {
                this->stripID_ = stripID;
            }
            
            /** 
             * Get the pad ID associated with a hit.  This is used to 
             * uniquely identify a sensor within a layer.
             * @return The pad ID associated with a hit.
             */
            int getPadID() const { return padID_; };
            
            /** 
             * Set the pad ID associated with a hit.  
             * @return padID The pad ID associated with a hit.
             */
            void setPadID(const int padID) { this->padID_ = padID; };

            /**
             * Get the energy deposition.
             * @return The energy deposition.
             */
            float getEdep() {
                return edep_;
            }

            /**
             * Get the energy .
             * @return The energy .
             */
            float getEnergy() {
                return energy_;
            }            

            /**
             * Set the energy deposition.
             * @param edep The energy deposition.
             */
            void setEdep(float edep) {
                this->edep_ = edep;
            }

            /**
             * Get the global time [ns].
             * @return The global time [ns].
             */
            float getTime() {
                return time_;
            }

            /**
             * Set the global time.
             * @param time The global time [ns].
             */
            void setTime(float time) {
                this->time_ = time;
            }

            /**
             * Get the XYZ momentum.
             * @return The XYZ momentum.
             */
            const G4ThreeVector& getMomentum() {
                return momentum_;
            }

            /**
             * Set the momentum.
             * @param px The X momentum.
             * @param py The Y momentum.
             * @param pz The Z momentum.
             */
            void setMomentum(float px, float py, float pz) {
                momentum_.setX(px);
                momentum_.setY(py);
                momentum_.setZ(pz);
            }

            /**
             * Get the XYZ hit position [mm].
             * @return The hit position.
             */
            const G4ThreeVector& getPosition() {
                return position_;
            }

            /**
             * Set the hit position [mm].
             * @param x The X position.
             * @param y The Y position.
             * 2param z The Z position.
             */
            void setPosition(float x, float y, float z) {
                position_.setX(x);
                position_.setY(y);
                position_.setZ(z);
            }

            /**
             * Set the energy.
             * @param edep The energy.
             */
            void setEnergy(float energy) {
                this->energy_ = energy;
            }            

            /**
             * Get the path length from the pre and post step points [mm].
             * @return The path length.
             */
            float getPathLength() {
                return pathLength_;
            }

            /**
             * Set the path length [mm].
             * @pathLength The path length.
             */
            void setPathLength(float pathLength) {
                this->pathLength_ = pathLength;
            }

        private:

            /**
             * The track ID.
             */
            G4int trackID_ {0};

            /**
             * The detector ID.
             */
            int id_ {0};

            /**
             * The detector ID.
             */
            int pdgid_ {0};
            /**
             * The layer ID.
             */
            int layerID_ {0};

            /** The pad ID. */
            int padID_{0}; 
            
            /**
             * The strip ID.
             */
            int stripID_ {0};

            /**
             * The energy deposition.
             */
            float edep_ {0};

            /**
             * The global time.
             */
            float time_ {0};

            /**
             * The XYZ momentum.
             */
            G4ThreeVector momentum_;

            /**
             * The XYZ position.
             */
            G4ThreeVector position_;

            /**
             * The energy .
             */
            float energy_ {0};

            /**
             * The path length.
             */
            float pathLength_ {0};
    };

    /**
     * Template instantiation of G4 hits collection class.
     */
    typedef G4THitsCollection<G4TriggerPadHit> G4TriggerPadHitsCollection;

    /**
     * Memory allocator for objects of this class.
     */
    extern G4Allocator<G4TriggerPadHit> G4TriggerPadHitAllocator;

    /**
     * Implementation of custom new operator.
     */
    inline void* G4TriggerPadHit::operator new(size_t) {
        void* aHit;
        aHit = (void*) G4TriggerPadHitAllocator.MallocSingle();
        return aHit;
    }

    /**
     * Implementation of custom delete operator.
     */
    inline void G4TriggerPadHit::operator delete(void *aHit) {
        G4TriggerPadHitAllocator.FreeSingle((G4TriggerPadHit*) aHit);
    }

}

#endif
