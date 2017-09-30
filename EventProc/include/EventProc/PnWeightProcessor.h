/**
 * @file pnWeightProcessor.h
 * @brief Processor that calculates pnWeight based on photonNuclear track 
 *        properties.
 * @author Alex Patterson, UCSB
 * @author Omar Moreno, SLAC National Accelerator Laboratory
 *
 * @note
 * pnWeightProcessor calculates an event weight which is added to the 
 * collection as a pnWeight object. This weight is based on simParticles
 * arising from photonNuclear reactions, and is intended to correct
 * the simulation in the case of high-momentum, backwards-going nucleons
 * arising from those reactions.
 *   fit variable W_p = 0.5*(p_tot + K)*(1.12-0.5*(p_z/p))
 *   where p_tot = sqrt(K^2 + 2*K*m)
 *           K = kinetic energy of nucleon at PN vertex
 *           p, p_z = momentum, z-component of nucleon at PN vertex
 */

#ifndef EVENTPROC_PNWEIGHTPROCESSOR_H_
#define EVENTPROC_PNWEIGHTPROCESSOR_H_

//----------------//
//   C++ StdLib   //
//----------------//
#include <iostream>
#include <cmath>
#include <algorithm>

//----------//
//   LDMX   //
//----------//
#include "Event/PnWeightResult.h"
#include "Event/SimParticle.h"
#include "Framework/EventProcessor.h"

//----------//
//   ROOT   //
//----------//
#include "TTree.h"
#include "TFile.h"
#include "TH1F.h"
#include "TClonesArray.h"

namespace ldmx {

    class PnWeightProcessor : public Producer {

        public:
            
            //---------------//
            //   Constants   //
            //---------------//

            /** Proton PDG ID */
            static const int PROTON_PDGID; 

            /** Neutron PDG ID */
            static const int NEUTRON_PDGID; 

            /** Constructor */
            PnWeightProcessor(const std::string& name, Process& process);

            /** Destructor */
            ~PnWeightProcessor();

            /**
             *  Read in user-specified parameters
             */
            void configure(const ParameterSet& pSet);

            /**
             *  Run the weight calculation and create a pnWeightResult
             */
            void produce(Event& event);

            double calculateWeight(double w); 

            /**
             * Calculate the measured W defined as
             *     W(measured) = 0.5*(p_tot + K)*(1.12-0.5*(p_z/p)) 
             * where 
             *     p is the total momentum of the particle
             *     K is its kinetic energy
             *     p_z is the z component of the momentum
             * all defined at the hardest PN vertex.
             *
             * @param particle SimParticle used to calculate W.
             * @return W
             */
            double calculateW(SimParticle* particle, double delta = 0.5);

        private:
        
            /** Threshold after which to apply W reweighting. */
            double wThreshold_{1400 /* MeV */};

            /**
             * Minimum angle for backwards-going hadron
             */
            double thetaThreshold_{100 /* degrees */};

            PnWeightResult result_;
    };
}

#endif
