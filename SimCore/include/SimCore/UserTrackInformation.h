/**
 * @file UserTrackInformation.h
 * @brief Class providing extra information associated to a Geant4 track
 * @author Jeremy McCormick
 * @author Omar Moreno, SLAC National Accelerator Laboratory
 */

#ifndef _SIMCORE_USER_TRACK_INFORMATION_H_
#define _SIMCORE_USER_TRACK_INFORMATION_H_

//------------//
//   Geant4   //
//------------//
#include "G4VUserTrackInformation.hh"
#include "G4ThreeVector.hh"

namespace ldmx {

    /**
     * @class UserTrackInformation
     * @note Provides extra information associated to a Geant4 track
     */
    class UserTrackInformation : public G4VUserTrackInformation {

        public:

            /** Overload pure virtual method (we don't implement it!). */
            void Print() const {;}

            /**
             * Get the flag which indicates whether this track should be saved
             * as a Trajectory.
             * @return The save flag.
             */
            inline bool getSaveFlag() const { return saveFlag_; }

            /**
             * Set the save flag so the associated track will be persisted
             * as a Trajectory.
             * @param saveFlag True to save the associated track.
             */
            inline void setSaveFlag(const bool& saveFlag) { saveFlag_ = saveFlag; }

            /**
             * Get the initial momentum 3-vector of the track [MeV].
             * @return The initial momentum of the track.
             */
            inline const G4ThreeVector& getInitialMomentum() { return initialMomentum_; } 
           
            /**
             * Set the initial momentum of the associated track.
             * @param p The initial momentum of the track.
             */ 
            inline void setInitialMomentum(const G4ThreeVector& p) { initialMomentum_.set(p.x(), p.y(), p.z()); } 
        
            /** 
             * Set the flag used by TargetBremFilter to tag photon tracks 
             * produced in the target via the eBrem process by electron 
             * satisfying an energy threshold requirement.
             * @param isBremCandidate True if the track is a brem whose parent
             *                        is an electron below threshold.
             */
            inline void tagBremCandidate(const bool& isBremCandidate = true) { isBremCandidate_ = isBremCandidate; }

            /**
             * @return true if the track is a brem candidate. 
             */
            inline bool isBremCandidate() const { return isBremCandidate_; } 


        private:

            /** Flag for saving the track as a Trajectory. */
            bool saveFlag_{false};

            /** Flag used to tag a track as a brem candidate. */
            bool isBremCandidate_{false}; 

            /** The initial momentum of the track. */
            G4ThreeVector initialMomentum_;

    };  // UserTrackInformation

} // ldmx

#endif // _SIMCORE_USER_TRACK_INFORMATION_H_
