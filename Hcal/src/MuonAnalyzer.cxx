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
            int maxconsecLayers(0), maxconsecStrips(0);
            for ( int s = 0; s < 5; s++ ) {
                int consecLayers = muonTrigger->getAlgoVar( 4*s+2 );
                int consecStrips = muonTrigger->getAlgoVar( 4*s+3 );
                double pathUnc = muonTrigger->getAlgoVar( 4*5 );
                
                hNumConsecLayers_[s]->Fill( consecLayers );
                hNumConsecStrips_[s]->Fill( consecStrips );
                
                if ( consecLayers > maxconsecLayers )
                    maxconsecLayers = consecLayers;

                if ( consecStrips > maxconsecStrips )
                    maxconsecStrips = consecStrips;
            }
            double pathUnc = muonTrigger->getAlgoVar( 20 );

            if ( pathUnc > 0.0 ) {
                if ( pathUnc > 10.0 )
                    pathUnc = 10.0;
                hUncertainPathLength_->Fill( pathUnc );
                hUncertainVConsecLayers_->Fill( maxconsecLayers , pathUnc );
                hUncertainVConsecStrips_->Fill( maxconsecStrips , pathUnc );
            }//skip uncalculated path uncs

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

        hUncertainPathLength_ = new TH1F( "hUncertainPathLength_" , "Uncertainty in Path Length" ,
            100 , 0.0 , 10.0 );

        hUncertainVConsecLayers_ = new TH2F( "hUncertainVConsecLayers_" , "Uncertainty in Path Length vs N Consecutive Layers" ,
            150 , 0.0 , 150.0 , 
            100 , 0.0 , 10.0 );

        hUncertainVConsecStrips_ = new TH2F( "hUncertainVConsecStrips_" , "Uncertainty in Path Length vs N Consecutive Strips" ,
            150 , 0.0 , 70.0 , 
            100 , 0.0 , 10.0 );
       
        return;
    }
}

DECLARE_ANALYZER_NS(ldmx, MuonAnalyzer);