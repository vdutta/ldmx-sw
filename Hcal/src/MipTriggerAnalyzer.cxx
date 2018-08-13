/**
 * @file MipTriggerAnalyzer.cxx
 * @brief Implementation file for MipTriggerAnalyzer class
 * @author Tom Eichlersmith, University of Minnesota
 */

#include "Hcal/MipTriggerAnalyzer.h"

namespace ldmx {

    void MipTriggerAnalyzer::configure(const ldmx::ParameterSet& ps) {
        
        hcalMipTriggerObjectName_ = ps.getString( "HcalMipTriggerObjectName" );

        hcalMipTriggerPassName_ = ps.getString( "HcalMipTriggerPassName" );

        numFalsePass_ = 0;
        numTruePass_ = 0;
        numFalseFail_ = 0;
        numTrueFail_ = 0;

        return;
    }

    void MipTriggerAnalyzer::analyze(const ldmx::Event& event) {
        
        //get list of triggers
        const TClonesArray *triggers = event.getCollection( "Trigger" , hcalMipTriggerPassName_ );
        //get list of actual sim particles
        const TClonesArray *simparticles = event.getCollection( "SimParticles" , "sim" );
        
        //get hcal mip trigger
        int ntriggers = triggers->GetEntriesFast();
        const TriggerResult *hcalMipTrigger;
        for ( int i = 0; i < ntriggers; i++ ) {
            hcalMipTrigger = (const TriggerResult *)(triggers->At(i));
            if ( hcalMipTrigger->getName() == hcalMipTriggerObjectName_ )
                break;
        }
        
        if ( hcalMipTrigger ) {
            hTracksPerEvent_->Fill( hcalMipTrigger->getAlgoVar4() );
        } else {
            std::cerr << hcalMipTriggerObjectName_ << " was not found in Trigger Collection in pass ";
            std::cerr << hcalMipTriggerPassName_ << std::endl;
            return;
        }

        //count number of actual muons
        int nmuons = 0;
        for ( int iP = 0; iP < simparticles->GetEntriesFast(); iP++ ) {
            SimParticle *simp = (SimParticle *)(simparticles->At(iP));
            if ( simp->getPdgID() == 13 or simp->getPdgID() == -13 ) {
                nmuons++;
            } //if a muon
        }

        //check accuracy of mip trigger
        bool triggerpass = hcalMipTrigger->passed();
        bool realpass = ( nmuons > 0 );
        if ( triggerpass and realpass)
            numTruePass_++;
        else if ( triggerpass and !realpass )
            numFalsePass_++;
        else if ( !triggerpass and realpass )
            numFalseFail_++;
        else
            numTrueFail_++;

        return;
    }

    void MipTriggerAnalyzer::onProcessStart() {
        
        getHistoDirectory();

        hTracksPerEvent_ = new TH1F( "hTracksPerEvent_" , "Tracks Found Per Event" ,
            11 , -0.5 , 10.5 );

        return;
    }

    void MipTriggerAnalyzer::onProcessEnd() {
        
        unsigned int numEvents = numTruePass_+numTrueFail_+numFalsePass_+numTrueFail_;
        double triggeraccuracy = (numTruePass_ + numTrueFail_)/(double)(numEvents);
        printf( "\n" );
        printf( " ===============================\n" );
        printf( " | Mip Trigger Confusion Table |\n" );
        printf( " | Mip     ||    Sim Particle  |\n" );
        printf( " | Trigger ||   Pass | Fail    |\n" );
        printf( " |    Pass ||%7d | %-7d |\n" , numTruePass_ , numFalsePass_ );
        printf( " |    Fail ||%7d | %-7d |\n" , numFalseFail_ , numTrueFail_ );
        printf( " |=============================|\n" );
        printf( " | Accuracy | %-16f |\n" , triggeraccuracy );
        printf( " ===============================\n" );

        return;
    }

}

DECLARE_ANALYZER_NS(ldmx, MipTriggerAnalyzer);
