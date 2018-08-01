/**
 * @file MipCluster.h
 * @author Tom Eichlersmith, University of Minnesota
 * @brief Header file for class MipCluster
 */

#ifndef HCAL_MIPCLUSTER_H
#define HCAL_MIPCLUSTER_H

//Standard Libraries
#include <vector> //for storage vector of HcalHit pointers

//ROOT Libs

//LDMX Framework
#include "DetDescr/HcalID.h" //for the HcalSection enum

namespace ldmx {
    
    /**
     * @class MipCluster
     * @brief Stores pointers to HcalHits that are considered a single MIP hit (usually due to proximity).
     */
    class MipCluster {
        public:
            /**
             * Default Constructor
             */
            MipCluster();

            /**
             * Add an HcalHit to the MipCluster
             */
            void addHit( ldmx::HcalHit* hit );
            
            /**
             * Get the total energy of the MipCluster
             */
            float getEnergy() const { return totalEnergy_; }

            /**
             * Get the number of HcalHits in this MipCluster
             */
            int getNumHits() const { return hcalHits_.size(); }

            /**
             * Get vector of HcalHits in this MipCluster
             */
            const std::vector< ldmx::HcalHit* > &getHcalHits() const { return hcalHits_; }

        private:
            
            /**
             * Calculate the total energy of the MipCluster.
             */
            void setTotalEnergy();
            
            /** The total energy of the MipCluster */
            float totalEnergy_;
 
            /** Storage vector of pointers to HcalHits */
            std::vector< ldmx::HcalHit* > hcalHits_;
   };

}

#endif /* HCAL_MIPCLUSTER_H */
