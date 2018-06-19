/**
 * @file HcalLayerAnalyzer.cxx
 * @brief Implementation file for HcalLayerAnalyzer
 * @author Tom Eichlersmith, University of Minnesota
 */

#include "playtest/HcalLayerAnalyzer.h"

namespace ldmx {

    void HcalLayerAnalyzer::configure(const ldmx::ParameterSet& ps) {
        caloCol_=ps.getString("caloHitCollection");
        minPE_ = static_cast<float>(ps.getDouble("minPE"));
    }

    void HcalLayerAnalyzer::analyze(const ldmx::Event& event) {
        const TClonesArray* tca=event.getCollection(caloCol_);
        

	    for (size_t i = 0; i < tca->GetEntriesFast(); i++) {
            
            //get current hit
            const ldmx::HcalHit* chit=(const ldmx::HcalHit*)(tca->At(i));
            float curr_PE = chit->getPE();

            //only process hits in the back hcal
            if ( chit->getSection() == ldmx::HcalSection::BACK ) {
                
                h_includedhits->Fill( curr_PE );

            } //only process non-noise hits in the back hcal
            else {
                nNotIncluded++;
            }
	    
        } //loop through all entries in calorimeter hits for current event (i)
        
        //Now hitlog has 

    }

    void HcalLayerAnalyzer::onProcessStart() {
        getHistoDirectory();
        h_includedhits = new TH1F("h_includedhits","PE Distribution of included hits",500,0.5,500.5);
        //Declare a histogram per layer of back hcal
        /*
        for ( ) {
            TH1F* h_newhistogram = new TH1F("name","title",500,0.0,1.0); //declare new histogram
            h_hitsperlayer_.pushback( h_newhistogram ); //put new histogram into vector
        } //loop through all layers in back hcal
        */

        nNotIncluded_ = 0;
        nStripsPerLayer_ = 1000; //NEED REAL NUMBER
        nLayers_ = 81; //NEED REAL NUMBER
    }

    void HcalLayerAnalyzer::onProcessEnd() {
        std::cout << "Number Hits NOT included in analysis: " << nNotIncluded_ << std::endl;
    }

}

DECLARE_ANALYZER_NS(ldmx, HcalLayerAnalyzer);
