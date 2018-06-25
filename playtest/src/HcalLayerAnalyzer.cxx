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
        nStrips_= ps.getInteger("nStrips");
        nEcalThickness_ = ps.getInteger("nEcalThickness");
        conedepth_ = ps.getInteger("conedepth");
        coneangle_ = ps.getInteger("coneangle");

        origin_ = static_cast<float>(nStrips_)/2;
        lowside_ = origin_ - static_cast<float>(nEcalThickness_)/2;
        upside_ = origin_ + static_cast<float>(nEcalThickness_)/2;
    }

    void HcalLayerAnalyzer::analyze(const ldmx::Event& event) {
        const TClonesArray* tca=event.getCollection(caloCol_);

	    for (size_t i = 0; i < tca->GetEntriesFast(); i++) {
            
            //get current hit
            const ldmx::HcalHit* chit=(const ldmx::HcalHit*)(tca->At(i));
            float curr_PE = chit->getPE();

            //only process non-noise hits in the back hcal
            if ( curr_PE > minPE_ and chit->getSection() == ldmx::HcalSection::BACK ) {
                
                h_includedhits->Fill( curr_PE );

            } //only process non-noise hits in the back hcal
            else {
                nNotIncluded_++;
            }
	    
        } //loop through all entries in calorimeter hits for current event (i)
        
        //Now log has the non-noise hits in it


    }

    void HcalLayerAnalyzer::onProcessStart() {
        getHistoDirectory();
        h_includedhits = new TH1F("h_includedhits","PE Distribution of included hits",500,0.5,500.5);

        nNotIncluded_ = 0;
        layermod_ = 1000;
    }

    void HcalLayerAnalyzer::onProcessEnd() {
        std::cout << "Number Hits NOT included in analysis: " << nNotIncluded_ << std::endl;
    }

}

DECLARE_ANALYZER_NS(ldmx, HcalLayerAnalyzer);
