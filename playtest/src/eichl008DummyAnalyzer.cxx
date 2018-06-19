/**
 * @file eichl008DummyAnalyzer.cxx
 * @brief Class that defines a dummy Analyzer implementation that just prints some messages. Using this file to practice adding a new analyzer to the framework.
 * @author Jeremy Mans, University of Minnesota 
 * @author Tom Eichlersmith, University of Minnesota
 */

#include "playtest/eichl008DummyAnalyzer.h" 

namespace ldmx {

    void eichl008DummyAnalyzer::configure(const ldmx::ParameterSet& ps) {
        caloCol_=ps.getString("caloHitCollection");
        keepMod_=ps.getInteger("keepEventModulus",0);
        dropMod_=ps.getInteger("dropEventModulus",0);
    }

    void eichl008DummyAnalyzer::analyze(const ldmx::Event& event) {
        const TClonesArray* tca=event.getCollection(caloCol_);

	    float total_energy = 0.0; //total energy from all hcal hits in this event
	    for (size_t i=0; i<tca->GetEntriesFast(); i++) {
            const ldmx::HcalHit* chit=(const ldmx::HcalHit*)(tca->At(i));
            total_energy += chit->getEnergy();
            h_pe->Fill(chit->getPE());
	    }
        h_energyperevent->Fill(total_energy);
            
            
	   	int ievent=event.getEventHeader()->getEventNumber();
        if (keepMod_>0 && !(ievent%keepMod_)) setStorageHint(hint_shouldKeep);
        if (dropMod_>0 && !(ievent%dropMod_)) setStorageHint(hint_shouldDrop);
         
        /* Failing attempt at hinting to drop events
	    if ( total_energy > 0.1 ) //total energy arbitrary cutoff
	   	    setStorageHint(ldmx::StorageControlHint::hint_mustDrop);
	    	
	    if (true) {	
	        setStorageHint(hint_shouldDrop);
	    } else {
	    	setStorageHint(hint_shouldKeep);
	    }
        */
    }
}

DECLARE_ANALYZER_NS(ldmx, eichl008DummyAnalyzer);
