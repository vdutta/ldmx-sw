/**
 * @file HcalMipTrackProducer.cxx
 * @brief Implementation file for HcalMipTrackProducer class.
 * @author Tom Eichlersmith, University of Minnesota
 */

#include "Hcal/HcalMipTrackProducer.h"

namespace ldmx {

    void HcalMipTrackProducer::configure(const ldmx::ParameterSet& ps) {

        hcalHitCollName_ = ps.getString( "HcalHitCollectionName" );

        hcalHitPassName_ = ps.getString( "HcalHitPassName" );

        maxTrackCount_ = 100;

        hcalMipTracks_ = new TClonesArray( "ldmx::HcalMipTrack" , maxTrackCount_ );

        hcalMipTracksCollName_ = ps.getString( "HcalMipTrackCollectionName" );

        minPE_ = ps.getDouble( "MinimumPE" );

        maxEnergy_ = ps.getDouble( "MaximumEnergy" );

        fracTotalClusters_ = ps.getDouble( "FractionTotalClusters" );

        if ( fracTotalClusters_ < 0.0 or fracTotalClusters_ > 1.0 ) {
            std::cerr << "[ HcalMipTrackProducer::configure ] : FractionTotalClusters is out of viable range!\n";
            std::cerr << "                                      Must be set within [0,1]\n";
            fracTotalClusters_ = 0.2;
        }
        
        fracClustersLeft_ = ps.getDouble( "FractionClustersLeft" );
        
        if ( fracClustersLeft_ < 0.0 or fracClustersLeft_ > 1.0 ) {
            std::cerr << "[ HcalMipTrackProducer::configure ] : FractionClustersLeft is out of viable range!\n";
            std::cerr << "                                      Must be set within [0,1]\n";
            fracClustersLeft_ = 0.8;
        }

        maxSlopeAngleDiff_ = ps.getDouble( "MaximumSlopeAngleDifference" );

        meanTime_produce_ = 0.0;
        meanNumTouchLogs_ = 0.0;
        meanClustersIgnored_ = 0.0;

        return;
    }

    void HcalMipTrackProducer::produce(ldmx::Event& event) {
        
        numTouchLogs_ = 0;
        std::clock_t start_produce;
        start_produce = std::clock();

        //Clear event containers
        hcalHitLog_.clear();
        clusterLog_.clear();

        const TClonesArray *rawhits = event.getCollection( hcalHitCollName_ , hcalHitPassName_ );

        //go through raw hits and ignore noise hits
        int nhits = rawhits->GetEntriesFast();
        for ( int iH = 0; iH < nhits; iH++ ) {
            numTouchLogs_++;
            HcalHit* chit = dynamic_cast<HcalHit*>(rawhits->At( iH ));

            if ( isNotNoise( chit ) ) {
                int section = chit->getSection();
                int layer = chit->getLayer();
                int strip = chit->getStrip();

                unsigned int key = section*1000*100 + layer*100 + strip;

                hcalHitLog_[ key ] = chit;
            } //if not noise hit

        } //iterate through rawhits (iH)
        
        clusterHits();

        std::vector< unsigned int > track_mipids;
        int trackcnt = 0;
        int nclustersintracks = 0;
        while ( findSeed( false ) and trackcnt < maxTrackCount_ ) {
            
            if ( buildTrack( track_mipids ) ) {

                //able to build track from seed 
                HcalMipTrack *track = new HcalMipTrack(); 
                for ( unsigned int mipid : track_mipids ) {
                    numTouchLogs_++;
                    MipCluster* cmip = &clusterLog_.at( mipid );
                     
                    //add HcalHits to track
                    for ( int i = 0; i < cmip->getNumHits(); i++ ) {
                        track->addHit( cmip->getHcalHit( i ) );
                    }//iterate through hits in cluster
                    
                    //add real point to track
                    HitBox cmip_box = cmip->getBox();
                    track->addPoint( cmip_box.getMin() , cmip_box.getOrigin() , cmip_box.getMax() );
                    
                    //erase mipid from log
                    clusterLog_.erase( mipid );

                } //add clusters with mipids to track
                
                //construct fit so we can compare to previous tracks for merging
                track->setFit();
                
                //check if track should be merged with a previous track
                int iT;
                for ( iT = 0; iT < hcalMipTracks_->GetEntriesFast(); iT++ ) {
                    HcalMipTrack *listedtrack = (HcalMipTrack *)(hcalMipTracks_->At( iT ));
                    if ( shouldMergeTracks( listedtrack , track ) ) {
                        //tracks look like two pieces of one track
                        // Use current iT for merging
                        break;
                    }
                }
                
                //put newly built track into tclonesarray
                // this may point to a new HcalMipTrack created at end of TClonesArray or an old one
                // iT == trackcnt if new
                HcalMipTrack *listedtrack = (HcalMipTrack *)(hcalMipTracks_->ConstructedAt( iT ));
                listedtrack->merge( track );
                delete track;
                
                nclustersintracks += track_mipids.size();

                if ( iT == trackcnt ) {
                    trackcnt++;
                } else if ( iT > trackcnt ) {
                    std::cerr << "[ HcalMipTrackProducer::produce ] : Iterated past the number of tracks!" << std::endl;
                }

            } else {
                //unable to build track, mark as bad seed     
                clusterLog_.at( seedID_ ).wasBadSeed();

            } //if able to build track

        } //repeat track construction until no more viable seeds
        
        //store collection in event bus
        event.add( hcalMipTracksCollName_ , hcalMipTracks_ ); 

        numTracksPerEvent_[ trackcnt ] ++;
        
        if ( meanClustersPerTrack_.count( trackcnt ) < 1 )
            meanClustersPerTrack_[ trackcnt ] = 0.0;
        
        double meanclusters = 0.0;
        if ( trackcnt > 0 )
            meanclusters = (double)(nclustersintracks)/(double)(trackcnt);
        double prev_mean = meanClustersPerTrack_.at(trackcnt);
        int nevents = numTracksPerEvent_.at( trackcnt );
        meanClustersPerTrack_[ trackcnt ] = ((double)(nevents-1)/(double)(nevents))*prev_mean +
            meanclusters/nevents;
        
        int ievent = event.getEventHeader()->getEventNumber();
        double time_produce = (std::clock() - start_produce)/(double)(CLOCKS_PER_SEC / 1000 ); //ms
        meanTime_produce_ = ((double)(ievent)/(double)(ievent+1))*meanTime_produce_ +
            time_produce/(double)(ievent+1);
        meanNumTouchLogs_ = ((double)(ievent)/(double)(ievent+1))*meanNumTouchLogs_ +
            numTouchLogs_/(double)(ievent+1);
        meanClustersIgnored_ = ((double)(ievent)/(double)(ievent+1))*meanClustersIgnored_ +
            clusterLog_.size()/(double)(ievent+1);

        return;
    }

    void HcalMipTrackProducer::onProcessEnd() {
                
        printf( "\n" );
        printf( " ==========================================\n" );
        printf( " |HcalMipTrackProducer - Performance Stats|\n");
        printf( " |========================================|\n");
        printf( " |                Stat : Mean             |\n");
        printf( " |        Time produce : %-10.8fs      |\n", meanTime_produce_/1000 );
        printf( " |         Log Touches : %-10.2f       |\n" , meanNumTouchLogs_ );
        printf( " |    Clusters Ignored : %-10.2f       |\n" , meanClustersIgnored_ );
        printf( " |========================================|\n" );
        printf( " | N Tracks : N Events : Mean N Clusters  |\n" );
        for ( auto keyval : numTracksPerEvent_ ) {
            printf( " |%9d : %-8d : %-16.2f |\n" , 
                keyval.first , keyval.second , meanClustersPerTrack_.at(keyval.first) );
        }
        printf( " ==========================================\n" );

        return;
    }

    bool HcalMipTrackProducer::isNotNoise( HcalHit* hit ) const {
        return ( !hit->getNoise() and hit->getPE() > minPE_ );
    }

    bool HcalMipTrackProducer::isMip( const MipCluster* cluster ) const {
        return ( cluster->getEnergy() < maxEnergy_ );
    }

    void HcalMipTrackProducer::clusterHits() {
        
        // cluster hits in same layer
        std::map< unsigned int , HcalHit* >::iterator itH;
        std::map< unsigned int , HcalHit* >::iterator prev_itH = hcalHitLog_.begin();
        MipCluster *current_cluster = new MipCluster();
        for ( itH = hcalHitLog_.begin(); itH != hcalHitLog_.end(); ++itH ) {
            numTouchLogs_++;
            //itH and prev_itH will both point to hcalHitLog_.begin() on first loop
            unsigned int keydif = itH->first - prev_itH->first;

            if ( keydif > 1 ) {
                //current hit is in different cluster
                // add previous cluster to log and clear temporary cluster
                current_cluster->setUID( prev_itH->first );
                current_cluster->set(); //calculate real space point
                if ( isMip( current_cluster ) ) {
                    clusterLog_[ current_cluster->getUID() ] = *current_cluster;
                }//check if mip
                delete current_cluster;
                current_cluster = new MipCluster();
            }//check if separate cluster

            current_cluster->addHit( itH->second );

            prev_itH = itH;

        }//iterate through sorted hcalHitLog (itH)

        //clean up at end of hit log
        if ( current_cluster ) {
            current_cluster->setUID( prev_itH->first );
            current_cluster->set();
            if ( isMip( current_cluster ) ) {
                clusterLog_[ current_cluster->getUID() ] = *current_cluster;
            }
            delete current_cluster;
            current_cluster = nullptr;
        }
        
        //set the total number of clusters
        minNumClusters_ = static_cast<int>( clusterLog_.size()*fracTotalClusters_ );

        return; 
    }

    bool HcalMipTrackProducer::findSeed( const bool useMedian ) {
        
        seedID_ = 0;

        if ( clusterLog_.size() > minNumClusters_ ) {

            std::map< const double , unsigned int > zpos_id;
            
            HitBox mipbox;
            std::vector<double> point;
            for ( auto keyclust : clusterLog_ ) {
                numTouchLogs_++;
                if ( (keyclust.second).isGoodSeed() ) {
                    mipbox = (keyclust.second).getBox();
                    point = mipbox.getOrigin();
                    zpos_id[ point[2] ] = keyclust.first;
                } //if not been a seed before
            } //go through clusterLog_
            
            std::map< const double , unsigned int >::iterator seed_it = zpos_id.begin();
            
            if ( seed_it != zpos_id.end() ) {
            
                if ( useMedian ) {
                    //move seed_it to rough media
                    int mid_ind = zpos_id.size()/2;
                    for ( int i = 0; i < mid_ind; i++ )
                        ++seed_it;
                }
    
                seedID_ = seed_it->second;
                seedBox_ = clusterLog_.at( seedID_ ).getBox();
                numTouchLogs_++;

            } //zpos_id is not empty

        } //enough clusters in log
        
        return (seedID_ != 0);
    }
 
    bool HcalMipTrackProducer::buildTrack( std::vector< unsigned int > &track_mipids ) {

        track_mipids.clear();
        
        std::vector< double > seedpoint = seedBox_.getOrigin();
        std::vector< double > seedmin = seedBox_.getMin();
        std::vector< double > seedmax = seedBox_.getMax();

        //iterate through all pairs of points
        HcalMipTrack best_track;
        std::map< unsigned int , MipCluster >::iterator itEnd, itC; //iterators for map
        for ( itEnd = clusterLog_.begin(); itEnd != clusterLog_.end(); ++itEnd ) {
            
            if ( itEnd->first != seedID_ ) {
                
                //construct track in cylinder
                HitBox endBox = (itEnd->second).getBox();
                std::vector< double > endpoint = endBox.getOrigin();
                std::vector< double > endmin = endBox.getMin();
                std::vector< double > endmax = endBox.getMax();

                //calculate line properties 
                //  (direction, negative direction, smudging factor)
                std::vector< double > direction( 3 , 0.0 ), negdirection( 3 , 0.0 );
                std::vector< double > linesmudge( 3 );
                for ( unsigned int iC = 0; iC < 3; iC++ ) {
                    direction[ iC ] = endpoint[iC] - seedpoint[iC];
                    negdirection[ iC ] = -1*direction[iC];
                    
                    //determine line smudge
                    std::vector< double > err_opts;
                    err_opts.push_back( abs( endpoint[iC] - endmin[iC] ) );
                    err_opts.push_back( abs( endpoint[iC] - endmax[iC] ) );
                    err_opts.push_back( abs( seedpoint[iC] - seedmin[iC] ) );
                    err_opts.push_back( abs( seedpoint[iC] - seedmax[iC] ) );
                    linesmudge[iC] = err_opts[0];
                    for ( unsigned int ie = 0; ie < 4; ie++ ) {
                        if ( err_opts[ie] < linesmudge[iC] )
                            linesmudge[iC] = err_opts[ie];
                    }

                } //space coordinates (iC)
                
                std::vector< unsigned int > ctrack_mipids; //ids of mip clusters in track
                //iterate through all clusters to see if they are in track
                for ( itC = clusterLog_.begin(); itC != clusterLog_.end(); ++itC ) {
                    numTouchLogs_++;
                    
                    HitBox mip_box = (itC->second).getBox();
                    
                    //see if ray hits box on either side of origin along line
                    if ( rayHitBox( seedpoint , direction    , mip_box.getMin() , mip_box.getMax() ) or
                         rayHitBox( seedpoint , negdirection , mip_box.getMin() , mip_box.getMax() ) ) {
                        ctrack_mipids.push_back( itC->first );
                    } 
    
                } //iterate through all clusters to see if they are in track (itC)
                
                //check if current track is acceptable and an improvement
                if ( isAcceptableTrack( ctrack_mipids ) and 
                     ctrack_mipids.size() > track_mipids.size() ) {
                    track_mipids = ctrack_mipids;
                }//ctrack is a acceptable and includes more clusters than other track
            
            }//make sure origin and end aren't the same

        } //go through remaining hits as end point (itEnd)
        
        return (!track_mipids.empty());
    }
   
    bool HcalMipTrackProducer::rayHitBox( const std::vector<double> origin , const std::vector<double> dir , 
                    const std::vector<double> minBox , const std::vector<double> maxBox ) const {
        
        bool originInside = true;
        bool originBetween[3];

        //Determine planes that are on the "front" of the box w.r.t. the origin of the ray
        std::vector<double> candidatePlane( 3 , 0.0 );
        for ( unsigned int iC = 0; iC < 3; iC++ ) {
            
            if ( origin.at(iC) < minBox.at(iC) ) {
                originBetween[iC] = false;
                candidatePlane[iC] = minBox.at(iC);
                originInside = false;
            } else if ( origin.at(iC) > maxBox.at(iC) ) {
                originBetween[iC] = false;
                candidatePlane[iC] = maxBox.at(iC);
                originInside = false;
            } else {
                originBetween[iC] = true;
            } //where origin is w.r.t. box
                
        } //iterate through coordinates (iC)
        
        //Origin Inside Box ==> Ray Intersects Box
        if ( originInside ) {
            return true;
        }

        //Calculate maximum T distances to candidatePlanes
        std::vector<double> maxT( 3 , 0.0 );
        for ( unsigned int iC = 0; iC < 3; iC++ ) {
            
            if ( !originBetween[iC] and dir.at(iC) != 0.0 ) {
                maxT[ iC ] = ( candidatePlane[iC] - origin.at(iC) ) / dir.at(iC);        
            } else {
                maxT[ iC ] = -1.0;
            }

        } //iterate through coordinates (iC)

        //Get largest of maxTs for the final choice of intersection
        unsigned int iMax = 0;
        for ( unsigned int iC = 0; iC < 3; iC++ ) {
            if ( maxT[ iMax ] < maxT[ iC ] ) {
                iMax = iC;
            }
        } //iterate through coordinates (iC)
        
        //Check if final candidate is inside box
        if ( maxT[ iMax ] < 0.0 ) {
            return false;
        }

        for ( unsigned int iC = 0; iC < 3; iC++ ) {
            
            if ( iMax != iC ) {
                double coordinate = origin.at(iC) + maxT[iC]*dir.at(iC);
                if ( coordinate < minBox.at(iC) or coordinate > maxBox.at(iC) ) {
                    //coordinate outside box
                    return false;
                }
            } //if coordinate is not maximum T plane
        } //iterate through coordinates (iC)

        return true;
    }

    bool HcalMipTrackProducer::isAcceptableTrack( const std::vector< unsigned int > &track_mipids ) const {
        return ( track_mipids.size() > fracClustersLeft_*clusterLog_.size() );
    }

    bool HcalMipTrackProducer::shouldMergeTracks( HcalMipTrack *first , HcalMipTrack *secon ) const {
       
        bool yes = false;

        std::vector<double> firstStart = first->getStart();
        std::vector<double> firstEnd = first->getEnd();
        std::vector<double> firstSlope( 3 , 0.0 );

        std::vector<double> seconStart = secon->getStart();
        std::vector<double> seconEnd = secon->getEnd();
        std::vector<double> seconSlope( 3 , 0.0 );

        double firstSlopeMag(0.0), seconSlopeMag(0.0);
        for ( int i = 0; i < 3; i++ ) {
            
            //slope of first and magnitude
            firstSlope[i] = firstEnd[i] - firstStart[i];
            firstSlopeMag += firstSlope[i]*firstSlope[i];
            
            //slope of secon and magnitude
            seconSlope[i] = seconEnd[i] - seconStart[i];
            seconSlopeMag += seconSlope[i]*seconSlope[i];
        }
        
        double crossProdMag = 0.0;
        for ( int i = 0; i < 3; i++ ) {
            
            int hi =  (i+1) % 3;

            double coordinate_val = (firstSlope[i]*seconSlope[hi]) +
                (firstSlope[hi]*seconSlope[i]);
            
            crossProdMag += coordinate_val*coordinate_val;
        }
        crossProdMag = sqrt(crossProdMag/(firstSlopeMag*seconSlopeMag));
        
        //check if angle between direction vectors is less than
        // parameter input by user
        double angledif = abs(asin( crossProdMag ));
        if ( angledif < maxSlopeAngleDiff_ ) {
            //merge tracks 
            yes = true;
        } 

        return yes;
    }
}

DECLARE_PRODUCER_NS( ldmx , HcalMipTrackProducer );
