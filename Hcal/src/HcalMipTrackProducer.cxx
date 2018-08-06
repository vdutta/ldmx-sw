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

        TClonesArray *hcalMipTracks_ = new TClonesArray( "ldmx::HcalMipTrack" , maxTrackCount_ );

        hcalMipTracksCollName_ = ps.getString( "HcalMipTrackCollectionName" );

        minPE_ = ps.getDouble( "MinimumPE" );

        maxEnergy_ = ps.getDouble( "MaximumEnergy" );

        minNumClusters_ = ps.getInteger( "MinimumNumClusters" );
        
        //tracks are started with a pair of points, so
        //  the absolute minimum number of clusters in a track is 2
        if ( minNumClusters_ < 2 ) {
            minNumClusters_ = 2;
        } //minNumClusters_ validation check
        
        return;
    }

    void HcalMipTrackProducer::produce(ldmx::Event& event) {
        
        const TClonesArray *rawhits = event.getCollection( hcalHitCollName_ , hcalHitPassName_ );

        //go through raw hits and ignore noise hits
        int nhits = rawhits->GetEntriesFast();
        for ( int iH = 0; iH < nhits; iH++ ) {
            
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
        while ( buildTrack( track_mipids ) and trackcnt < maxTrackCount_ ) {
            //store best in collection and delete used mips from log
            std::cout << trackcnt << std::endl;
            HcalMipTrack *track = (HcalMipTrack *)(hcalMipTracks_->ConstructedAt( trackcnt ));
            
            for ( unsigned int mipid : track_mipids ) {
                
                MipCluster* cmip = &clusterLog_[ mipid ];

                //add HcalHits to track
                for ( int i = 0; i < cmip->getNumHits(); i++ ) {
                    track->addHit( cmip->getHcalHit( i ) );
                }//iterate through hits in cluster
                
                //add real point to track
                std::vector<double> point, errors;
                cmip->getPoint( point , errors );
                track->addPoint( point , errors );
                
                //erase mipid from log
                clusterLog_.erase( mipid );

            } //add clusters with mipids to track
            
            trackcnt++;

        } //repeat track construction until no more pairs of clusters
        
        event.add( hcalMipTracksCollName_ , hcalMipTracks_ ); 
        
        numTracksPerEvent_[ trackcnt ] ++;

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

        return; 
    }
 
    bool HcalMipTrackProducer::buildTrack( std::vector< unsigned int > &track_mipids ) {
        //for clusters:
        //  no suffix means that it isn't an endpoint
        //  suffice {1,2} means that it is one of the endpoints
        track_mipids.clear();

        //iterate through all pairs of points
        HcalMipTrack best_track;
        std::map< unsigned int , MipCluster >::iterator itC1, itC2, itC; //iterators for map
        for ( itC1 = clusterLog_.begin(); itC1 != clusterLog_.end(); ++itC1 ) {
            for ( itC2 = std::next(itC1); itC2 != clusterLog_.end(); ++itC2 ) {
                //construct track in cylinder
                
                std::vector< double > point1, point2 , errors1 , errors2;
                (itC1->second).getPoint( point1 , errors1 );
                (itC2->second).getPoint( point2 , errors2 );
                
                //calculate line properties 
                //  (origin, direction, negative direction, smudging factor)
                std::vector< double > origin( point1 ), direction( 3 , 0.0 );
                std::vector< double > negdirection( 3 , 0.0 ), linesmudge( errors1 );
                for ( unsigned int iC = 0; iC < 3; iC++ ) {
                    direction[ iC ] = point2[iC] - point1[iC];
                    negdirection[ iC ] = -1*direction[iC];

                    //determine line smudging
                    if ( errors2[iC] < linesmudge[iC] ) {
                        linesmudge[iC] = errors2[iC];
                    }
                } //space coordinates (iC)
                
                std::vector< unsigned int > ctrack_mipids; //ids of mip clusters in track
                //iterate through all clusters to see if they are in track
                for ( itC = clusterLog_.begin(); itC != clusterLog_.end(); ++itC ) {
                    
                    std::vector< double> point, errors;
                    (itC->second).getPoint( point , errors );
                    
                    //construct hit box
                    // could add a fudge factor controlled by user
                    std::vector< double > maxBox( 3 , 0.0 ), minBox( 3 , 0.0 );
                    for ( unsigned int iC = 0; iC < 3; iC++ ) {
                        maxBox[iC] = point[iC] + errors[iC] + linesmudge[iC];
                        minBox[iC] = point[iC] - errors[iC] - linesmudge[iC];
                    }
                    
                    //see if ray hits box on either side of origin along line
                    if ( rayHitBox( origin , direction    , minBox , maxBox ) or
                         rayHitBox( origin , negdirection , minBox , maxBox ) ) {
                        ctrack_mipids.push_back( itC->first );
                    } 

                } //iterate through all clusters to see if they are in track (itC)
                
                //check if plausible track
                if ( ctrack_mipids.size() > minNumClusters_ ) {
                    //create fit for ctrack
                    HcalMipTrack ctrack;
                    for ( std::vector< unsigned int >::iterator it = ctrack_mipids.begin();
                        it != ctrack_mipids.end(); ++it ) {
                        
                        MipCluster* cmip = &clusterLog_[ *it ];
                        for ( int i = 0; i < cmip->getNumHits(); i++ ) {
                            ctrack.addHit( cmip->getHcalHit( i ) );
                        }//iterate through hits in cluster

                        std::vector<double> point,errors;
                        cmip->getPoint( point , errors );
                        ctrack.addPoint( point , errors );

                    } //add clusters with mipids to ctrack
                    
                    if ( compMipTracks( best_track , ctrack ) ) {
                        best_track = ctrack;
                        track_mipids = ctrack_mipids;
                    }//ctrack is better than best_track

                }//ctrack is a plausible track

            } //go through remaining hits as second end point (itC2)
        } //go through all hits as first end point (itC1)
        
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

    bool HcalMipTrackProducer::compMipTracks( const HcalMipTrack &track1 , const HcalMipTrack &track2 ) const {
        
        bool better = false;
        if ( track1.isEmpty() ) {
            better = true;
        } else {
            //TEMPORARY
            better = ( track1.getNHits() < track2.getNHits() );
        } //if track1 is empty

        return better;
    }
}

DECLARE_PRODUCER_NS( ldmx , HcalMipTrackProducer );
