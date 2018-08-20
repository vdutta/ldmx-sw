/**
 * @file MuonAnalyzer.cxx
 * @brief Implementation file for MuonAnalyzer class.
 * @author Tom Eichlersmith, University of Minnesota
 */

#include "Hcal/MuonAnalyzer.h"

namespace ldmx {

    void MuonAnalyzer::configure(const ldmx::ParameterSet& ps) {
        
        triggerObjectName_ = ps.getString( "TriggerObjectName" );

        triggerPassName_ = ps.getString( "TriggerPassName" );

        return;
    }

    void MuonAnalyzer::analyze(const ldmx::Event& event) {
        
        //get trigger array
        const TClonesArray *triggers = event.getCollection( "Trigger" , triggerPassName_ );

        const TriggerResult *muonTrigger;
        int nTriggers = triggers->GetEntriesFast();
        int iT;
        for ( iT = 0; iT < nTriggers; iT++ ) {
            muonTrigger = (const TriggerResult *)(triggers->At(iT));
            if ( muonTrigger->getName() == triggerObjectName_ )
                break;
        }

        if ( iT >= nTriggers ) {
            std::cerr << "WARNING [ MuonAnalyzer ] : Could not find trigger object '" << triggerObjectName_ << "'." << std::endl;
        } else {
    
            //fill histograms for each section
            for ( int s = 0; s < 5; s++ ) {
                hNumConsecLayers_[s]->Fill( muonTrigger->getAlgoVar( 4*s+2 ) );
                hNumConsecStrips_[s]->Fill( muonTrigger->getAlgoVar( 4*s+3 ) );
            }

        } //check if found trigger

        return;
    }

    void MuonAnalyzer::onProcessStart() {
        
        getHistoDirectory();
        
        for ( int s = 0; s < 5; s++ ) {

            std::string sect = std::to_string(s);

            hNumConsecStrips_[s] = new TH1F( ("hNumConsecStrips_"+triggerObjectName_+sect).c_str() , ("Num Consecutive Strips in Hcal Section "+sect).c_str() ,
                40 , 0.0 , 40.0 );
            hNumConsecLayers_[s] = new TH1F( ("hNumConsecLayers_"+triggerObjectName_+sect).c_str() , ("Num Consecutive Layers in Hcal Section "+sect).c_str() ,
                150 , 0.0 , 150.0 );
        }
        
        
        return;
    }
}

DECLARE_ANALYZER_NS(ldmx, MuonAnalyzer);
