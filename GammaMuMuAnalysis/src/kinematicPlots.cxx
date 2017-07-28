/**
 * @file kinematicPlots.cxx
 * @brief Analysis of Gamma->MuMu event kinematics
 * @author Andrew Whitbeck, fermilab
 */

#include "Framework/EventProcessor.h"
#include "Framework/ParameterSet.h"

#include "TH1.h"
#include "TH2.h"
#include "Event/Event.h"
#include "Event/EcalHit.h"
#include "Event/SimParticle.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"

#include <iostream>
#include <cmath>

namespace ldmx {

    /**
     * @class kinematicPlots
     * @brief analyzer which isolates gamma->mu+mu- sim particles and computes kinematics for validation
     */
    class kinematicPlots : public ldmx::Analyzer {
        public:

            kinematicPlots(const std::string& name, ldmx::Process& process) : ldmx::Analyzer(name, process) {}

            virtual void configure(const ldmx::ParameterSet& ps) {
                particleCol_=ps.getString("simParticleCollection");
            }

            virtual void analyze(const ldmx::Event& event) {
                const TClonesArray* tca = event.getCollection("ecalDigis");
                double cal_e=0.;
                double digi_e=0.;
                int digi_layer=0;
                for(size_t h=0; h<tca->GetEntriesFast(); h++){
                    ldmx::EcalHit* hits = (ldmx::EcalHit*)(tca->At(h));
                    //std::cout << "ecal energy: " << hits->getEnergy() << std::endl;
                    digi_e = hits->getEnergy();
                    digi_layer = hits->getLayer();
                    cal_e+=digi_e/0.130*calib_weights[digi_layer]+digi_e*0.948;
                }
                //std::cout << "calibrated energy: " << cal_e << std::endl;
                h_calibrated_energy->Fill(cal_e);

                tca=event.getCollection(particleCol_);
                ldmx::SimParticle* gamma=0,*mu_plus=0,*mu_minus=0;
                for (size_t i=0; i<tca->GetEntriesFast(); i++) {
                    ldmx::SimParticle* p=(ldmx::SimParticle*)(tca->At(i));
                    if( abs(p->getPdgID()) == 22 && p->getDaughterCount() == 2 && abs(p->getDaughter(0)->getPdgID()) == 13 && abs(p->getDaughter(1)->getPdgID()) == 13 ){
                        gamma = p;
                        if( p->getDaughter(0)->getCharge() == 1 ){
                            mu_plus = p->getDaughter(0);
                            mu_minus = p->getDaughter(1);
                        }else{
                            mu_plus = p->getDaughter(1);
                            mu_minus = p->getDaughter(0);
                        }// end if/else block to determine muon charges
                        break;
                    }// end if gamma/mu+/mu-
                }// end for loop over sim particles

                if( gamma == 0 && mu_plus == 0 && mu_minus == 0 ) return;
                TLorentzVector gamma_p4(gamma->getEnergy(),gamma->getMomentum()[0],gamma->getMomentum()[1],gamma->getMomentum()[2]);
                TLorentzVector mu_plus_p4(mu_plus->getEnergy(),mu_plus->getMomentum()[0],mu_plus->getMomentum()[1],mu_plus->getMomentum()[2]);
                TLorentzVector mu_minus_p4(mu_minus->getEnergy(),mu_minus->getMomentum()[0],mu_minus->getMomentum()[1],mu_minus->getMomentum()[2]);
                h_muon_q2->Fill(log10(2)+log10((mu_plus_p4+mu_minus_p4).P()));
                h_muonPhoton_q2->Fill(log10(2)+log10((mu_plus_p4+mu_minus_p4-gamma_p4).P()));
                h_muonPhoton_q2_v_photonEnergy->Fill(log10(2)+log10((mu_plus_p4+mu_minus_p4-gamma_p4).P()),gamma_p4.E());
                h_muon_theta->Fill(mu_plus_p4.Theta());  
                h_muon_theta->Fill(mu_minus_p4.Theta());
                h_muon_phi->Fill(mu_plus_p4.Phi());
                h_muon_phi->Fill(mu_minus_p4.Phi());
                h_muon_pt->Fill(mu_plus_p4.Pt());
                h_muon_pt->Fill(mu_minus_p4.Pt());
                h_muon_e->Fill(mu_plus_p4.E());
                h_muon_e->Fill(mu_minus_p4.E());
                h_muon_pz->Fill(mu_plus_p4.Pz());
                h_muon_pz->Fill(mu_minus_p4.Pz());
                
                h_xP->Fill(mu_plus_p4.E()/gamma_p4.E());
                h_xM->Fill(mu_minus_p4.E()/gamma_p4.E());
                double uP = mu_plus_p4.Gamma()*mu_plus_p4.Vect().Angle(gamma_p4.Vect());
                double uM = mu_minus_p4.Gamma()*mu_minus_p4.Vect().Angle(gamma_p4.Vect());
                double u = (uP+uM)/2.;
                double xi = uP-uM;
                double beta = u*mu_plus_p4.DeltaPhi(mu_minus_p4);
                h_t->Fill(1./(1.+u*u));
                h_psi->Fill(TMath::ATan2(beta,xi));
                h_rho->Fill(sqrt(xi*xi+beta*beta));
                h_beta->Fill(beta);
                h_xi->Fill(xi);
            }
            virtual void onFileOpen() {
                std::cout << "kinematicPlots: Opening a file!" << std::endl;
            }
            virtual void onFileClose() {
                std::cout << "kinematicPlots: Closing a file!" << std::endl;
                h_muon_q2->Write();
                h_muonPhoton_q2->Write();
            }
            virtual void onProcessStart() {
                std::cout << "kinematicPlots: Starting processing!" << std::endl;

                calib_weights = {1.641, 3.526, 5.184, 6.841, 8.222, 8.775, 8.775, 8.775, 8.775, 8.775, 8.775, 8.775, 8.775, 8.775, 8.775, 8.775, 8.775, 8.775, 8.775, 8.775, 8.775, 8.775, 12.642, 16.51, 16.51, 16.51, 16.51, 16.51, 16.51, 16.51, 16.51, 16.51, 16.51, 16.51, 16.51, 16.51, 16.51, 16.51, 16.51, 8.45};

                getHistoDirectory();
                h_muon_theta = new TH1F("h_muon_theta",";#theta_{#mu};Events",100,0.,2*3.1415);
                h_muon_phi = new TH1F("h_muon_phi",";#Phi_{#mu};Events",100,0.,2*3.1415);
                h_muon_pt = new TH1F("h_muon_pt",";p_{T,#mu};Events",100,0.,4000.);
                h_muon_e = new TH1F("h_muon_e",";#E_{#mu};Events",100,0.,4000.);
                h_muon_pz = new TH1F("h_muon_pz",";#p_{z,#mu};Events",100,0.,4000.);

                h_muon_q2=new TH1F("muon_q2",";log_{10}(Q^{2});Events",500,0.0,10.0);
                h_muonPhoton_q2=new TH1F("muonPhoton_q2",";log_{10}(Q^{2});Events",500,0.0,10.0);
                h_muonPhoton_q2_v_photonEnergy=new TH2F("muonPhoton_q2_v_photonEnergy",";log_{10}(Q^{2});E_{#gamma}",500,0.,10.,500,0.0,4000.0);
                
                h_calibrated_energy = new TH1F("h_calibrated_energy",";Energy [Mev];Events",200,0.,5000.);
                h_xP = new TH1F("h_xP",";x+;count",50,0,1);
                h_xM = new TH1F("h_xM",";x-;count",50,0,1);
                h_t = new TH1F("h_t",";t;count",50,0,2);
                h_psi = new TH1F("h_psi",";#psi;count",50,-2,2);
                h_rho = new TH1F("h_rho",";#rho;count",50,0,2);
                h_beta = new TH1F("h_beta",";#beta;count",50,0,2);
                h_xi = new TH1F("h_xi",";#xi;count",50,0,2);
           }
            virtual void onProcessEnd() {
                std::cout << "kinematicPlots: Finishing processing!" << std::endl;
            }
        private:
            TH1 *h_muon_theta, *h_muon_phi, *h_muon_pt, *h_muon_e, *h_muon_pz;
            TH1 *h_xP, *h_xM, *h_t, *h_psi, *h_rho,*h_beta,*h_xi;
            TH1 *h_muon_q2,*h_muonPhoton_q2;
            TH2 *h_muonPhoton_q2_v_photonEnergy;
        TH1 *h_calibrated_energy;
            std::string particleCol_;
            std::vector<double> calib_weights;

    };
}

DECLARE_ANALYZER_NS(ldmx, kinematicPlots);
