#include "redux/momfbd/constraints.hpp"

#include "redux/logging/logger.hpp"
#include "redux/momfbd/momfbdjob.hpp"
#include "redux/file/fileana.hpp"
#include "redux/math/linalg.hpp"
#include "redux/util/cache.hpp"
#include "redux/util/datautil.hpp"

#include <algorithm>
#include <iostream>
#include <numeric>

#include <boost/functional/hash.hpp>
#include <boost/iterator/counting_iterator.hpp>

using namespace redux::momfbd;
using namespace redux::file;
using namespace redux::logging;
using namespace redux::math;
using namespace redux::util;
using namespace redux;

using namespace std;


void Constraints::Constraint::replaceIndices( const map<int32_t, int32_t>& indexMap ) {
    map<int32_t, int8_t> newEntries;
    for( auto & entry : entries ) {
        newEntries[indexMap.at( entry.first )] = entry.second;
    }
    entries = std::move( newEntries );
}


bool Constraints::Constraint::operator>( const Constraint& rhs ) const {

    auto it = entries.begin();
    auto rhsit = rhs.entries.begin();
    int valueOrder(0);
    while( it != entries.end() ) {
        if( it->first != rhsit->first ) return ( it->first < rhsit->first ); // sort on leftmost entry
        else if( !valueOrder && it->second != rhsit->second ) {

            valueOrder = (it->second > rhsit->second)?1:-1;                     // if the indices are all the same, sort on largest (leftmost) value
        }
        it++;
        rhsit++;
    }
    if( entries.size() != rhs.entries.size() ) return ( entries.size() > rhs.entries.size() );
    if( valueOrder ) return (valueOrder>0);

    return false;
}


Constraints::NullSpace::NullSpace( const std::map<int32_t, int8_t>& e, int32_t np, int32_t nc, bool old_ns ) :
    c_entries(e), entriesHash(0), nParameters(np), nConstraints(nc) {
    
    entriesHash = boost::hash_range( c_entries.begin(), c_entries.end() );
    string path_name = "Nullspace_";
    if( old_ns ) path_name += "old_";
    path_name +=  to_string(np)+"x"+to_string(nc)+"_"+to_string(entriesHash);
    setPath( path_name );
    isLoaded = false;
    
}


#define NS_THRESHOLD 1E-10
void Constraints::NullSpace::mapNullspace( void ) {
    ns_entries.clear();
    int nCols = nParameters - nConstraints;
    int count = 0;
    for( auto& value: ns ) {
        if( abs(value) > NS_THRESHOLD ) {
            ns_entries.insert(make_pair(PointI(count/nCols,count%nCols),value));
        }
        count++;
    }
}
#undef NS_THRESHOLD


void Constraints::NullSpace::calculateNullspace( logging::Logger& logger, bool store, bool old_ns ) {

    int32_t nFreeParams = nParameters-nConstraints;
    if( nFreeParams < 1 ) {
        ns.clear();
        return;
    } else {
        if(!old_ns) {
            LOG << "Constructing optimal nullspace basis for (" << nConstraints << "x" << nParameters << ") group." << ende;
            ns.resize( nParameters, nFreeParams );
            if( nFreeParams == 1 ) {    // 1D nullspace (i.e. non-tilts)
                ns = 1.0/sqrt(nParameters);
            } else {
                Array<int32_t> C( nConstraints, nParameters );
                C.zero();
                ns.zero();
                shared_ptr<double*> tmpNS = ns.reshape( nFreeParams, nParameters ); // we will work in row-major form, then transpose at the end
                double** nsArr = tmpNS.get();
                
                for( auto & entry: c_entries ) {
                    C( entry.first/nParameters, entry.first%nParameters ) = entry.second;
                }
                int32_t row_index(0);
                int32_t max_element_count(0);
                for( int32_t i(0); i<nConstraints; ++i ) {
                    int32_t this_element_count(0);
                    for( int32_t j(0); j<nParameters; ++j ) this_element_count += (C( i, j ) != 0);
                    if( this_element_count > max_element_count ) {
                        max_element_count = this_element_count;
                        row_index = i;
                    }
                }
                // The constraint just selected will be the constraint that the sum of all tilts is 0, i.e.
                // it will be a row consisting of max_element_count x 1 and the rest 0
                
                int32_t block_length(1);
                size_t basis_count(0);
                while( block_length < max_element_count ) {
                    int32_t offset(0);
                    while( offset+block_length < max_element_count ) {
                        int32_t rem = max(min(max_element_count-offset-block_length,block_length),0);
                        for( int32_t i(0); i<block_length; ++i ) {
                            nsArr[basis_count][offset+i] = 1.0;
                            if( offset+block_length+i < max_element_count ) {
                                nsArr[basis_count][offset+block_length+i] = -block_length*1.0/rem;
                            }
                        }
                        basis_count++;
                        offset += 2*block_length;
                    }
                    block_length *= 2;
                }
                for( int32_t i(0); i<nFreeParams; ++i ) {    // loop over nullspace base vectors
                    for( int32_t j(0); j<max_element_count; ++j ) {
                        double ns_val = nsArr[i][j];
                        if( fabs(ns_val) > 0 ) {;
                            for( int32_t k(0); k<nConstraints; ++k ) {
                                int32_t v1 = C(k,j);
                                if( (k != row_index) && v1 ) {    // skip the constraint we already used to construct the nullspace
                                    // find the other non-zero value (there will only be 1)
                                    for( int32_t p(0); p<nParameters; ++p ) {
                                        if( p == j ) continue;
                                        int32_t v2 = C(k,p);
                                        if( v2 ) nsArr[i][p] = -v1*ns_val/v2;
                                    }
                                }
                            }
                        }
                    }
                }
                for( int32_t i(0); i<nFreeParams; ++i ) {    // loop over nullspace base vectors
                    double sum(0);
                    for( int32_t k=0; k<nParameters; ++k ) {
                        sum += nsArr[i][k]*nsArr[i][k];
                    }
                    for( int32_t j(i+1); j<nFreeParams; ++j ) {
                        double scalar_product(0);
                        for( int32_t k=0; k<nParameters; ++k ) {
                            scalar_product += nsArr[i][k]*nsArr[j][k];
                        }
                        if( fabs(scalar_product) > 0.0 ) {
                            scalar_product /= sum;
                            std::transform( nsArr[j], nsArr[j]+nParameters, nsArr[i], nsArr[j], [&]( double& a, double& b ){ return a-scalar_product*b; } );
                        }
                    }
                    if( sum > 0 && sum != 1.0 ) {
                        sum = 1.0/sqrt(sum);
                        std::transform( nsArr[i], nsArr[i]+nParameters, nsArr[i], [&]( double& a ){ return a*sum; } );
                    }
                }
                transpose( ns.get(), nFreeParams, nParameters );
            }
            
        } else {
            LOG << "Calculating QR-decomposition for (" << nConstraints << "x" << nParameters << ") group." << ende;
            Array<double> C,R;
            C.resize(nConstraints, nParameters);
            ns.resize(nParameters, nParameters);    // use ns rather than a temporary "Q" array
            R.resize(nParameters, nConstraints);
            C.zero();
            for( auto & entry: c_entries ) {
                C( entry.first/nParameters, entry.first%nParameters ) = entry.second;
            }
            // we need C^T to find the nullspace of C
            transpose(C.get(), nConstraints, nParameters);
            qr_decomp(C.get(), nParameters, nConstraints, ns.get(), R.get());
            ns.setLimits(0,nParameters-1,nConstraints,nParameters-1);
            ns.trim();
        }
    }
    isLoaded = true;
    if(store) cacheStore();
    
}


bool Constraints::NullSpace::verify( logging::Logger& logger, const std::map<int32_t, int8_t>& e, int32_t nP, int32_t nC ) {
    
    if(nP != nParameters || nC != nConstraints || c_entries != e) return false;
    
    if( !isLoaded )  {
        if( cacheLoad() ) {
            if( nP == nParameters && nC == nConstraints && c_entries == e ) {
                LOG << "Found matching constraint-group in file: " << itemPath << ende;
                return true;
            }
        } else {
            LOG << "No matching constraint-group found." << ende;
            ns.resize(0);
            ns_entries.clear();
        }
    }
        
    return true;
}


size_t Constraints::NullSpace::csize(void) const {
    size_t sz = sizeof(entriesHash) + sizeof(nParameters) + sizeof(nConstraints);
    sz += c_entries.size()*(sizeof(int32_t)+sizeof(int8_t)) + sizeof(uint64_t);
    sz += ns.size();
    return sz;
}


uint64_t Constraints::NullSpace::cpack(char* ptr) const {
    using redux::util::pack;
    uint64_t count = pack(ptr, entriesHash);
    count += pack(ptr+count,nParameters);
    count += pack(ptr+count,nConstraints);
    count += pack(ptr+count,c_entries);
    count += ns.pack(ptr+count);
    return count;
}


uint64_t Constraints::NullSpace::cunpack(const char* ptr, bool swap_endian) {
    using redux::util::unpack;
    uint64_t count = unpack(ptr, entriesHash, swap_endian);
    count += unpack(ptr+count,nParameters, swap_endian);
    count += unpack(ptr+count,nConstraints, swap_endian);
    count += unpack(ptr+count,c_entries, swap_endian);
    count += ns.unpack(ptr+count, swap_endian);
    return count;
}


void Constraints::NullSpace::cclear(void) {
    ns.resize(0);
    c_entries.clear();
}


bool Constraints::NullSpace::operator<(const NullSpace& rhs) const {
    if( entriesHash == rhs.entriesHash) {
        return (c_entries < rhs.c_entries);
    }
    return (entriesHash < rhs.entriesHash);
}


Constraints::Group::Group( logging::Logger& l, shared_ptr<Constraint>& con, bool ons ) : nParameters(0), entriesHash(0), logger(l), old_ns(ons) {
    add( con );
}


void Constraints::Group::add( const shared_ptr<Constraint>& con ) {
    if( find( constraints.begin(), constraints.end(), con ) == constraints.end() ) { // don't add the same row twice
        constraints.push_back( con );
        for( auto & entry : con->entries ) {
            indices.insert( entry.first );
        }
    }
}


void Constraints::Group::addConnectedConstraints( vector<shared_ptr<Constraint>>& cons ) {
    for( auto it = cons.begin(); it != cons.end(); ) {
        bool connected( false );
        for( auto& entry : ( *it )->entries ) {  // if one of the column-indices matches this group, add the constraint to group.
            if( indices.find( entry.first ) != indices.end() ) {
                connected = true;
            }
        }
        if( connected ) {
            add( *it );
            cons.erase( it );       // remove this row
            it = cons.begin();      // restart loop (since "indices" might contain new values)
        } else it++;
    }
}


/*
 *  Swap column order so that this group has all entries in consecutive columns
 *  This makes the constraint matrix block-diagonal.
 */
void Constraints::Group::blockify( int32_t* columnOrdering, int32_t& cOffset, int32_t& rOffset ) {

    groupOffset.y = rOffset;    // store the offset of this group (in the nullspace matrix).
    groupOffset.x = cOffset;
    
    nParameters = indices.size();
    int32_t nConstraints = constraints.size();
    
    vector<int32_t> newIndices(nParameters);
    std::iota( newIndices.begin(), newIndices.end(), rOffset );         // newIndices is a dense range starting at rOffset
    
    vector<int32_t> entry_count( nParameters, 0 );
    struct my_comp {
        bool operator() ( const PointI& lhs, const PointI& rhs ) const {
            if( lhs.y != rhs.y ) return lhs.y>rhs.y;
            return lhs.x<rhs.x;
        }
    };
    set<PointI,my_comp> entryCounts;
    for( auto& i: indices ) {
        PointI p(0, i);
        for( auto &c: constraints ) {
            if( c->has(i) ) p.y++;
        }
        entryCounts.insert(p);
    }

    map<int32_t, int32_t> indexMap;
    size_t cnt(0);
    for( auto& ec: entryCounts ) {
        columnOrdering[newIndices[cnt]] = ec.x;
        indexMap[ec.x] = newIndices[cnt++];
    }

    indices.clear();
    indices.insert( newIndices.begin(), newIndices.end() );
    
    std::map<int32_t, int8_t> tmpEntries;
    int32_t indexOffset = 0;
    for( auto & constraint : constraints ) {
        constraint->replaceIndices( indexMap );
        for( auto & ind : constraint->entries ) {
            tmpEntries.insert(make_pair(indexOffset-rOffset+ind.first, ind.second));
        }
        indexOffset += nParameters;     // next row
    }
    entries = std::move( tmpEntries );
    entriesHash = boost::hash_range(entries.begin(), entries.end());
    
    rOffset += nParameters;
    cOffset += (nParameters - nConstraints);

    mapNullspace();
    
}


void Constraints::Group::sortRows(void) {

    // sort rows by number of first entry -> nEntries -> value
    std::sort( constraints.begin(), constraints.end(),
            []( const shared_ptr<Constraint> &lhs, const shared_ptr<Constraint> &rhs ) {
                return ( (*lhs) > (*rhs) );
            }
        );
    
}


bool Constraints::Group::dense( void ) const {

    int32_t first = *indices.begin();
    int32_t last = *indices.rbegin();
    return ( ( last - first + 1 ) == (int32_t)indices.size() );

}


void Constraints::Group::mapNullspace(void) {
    
    int32_t nConstraints = constraints.size();
    ns_entries.clear();
    if( nParameters > nConstraints ) {
        
        shared_ptr<NullSpace> tmpNS( new NullSpace( entries, nParameters, nConstraints, old_ns ) );
        nullspace = redux::util::Cache::get( entriesHash, tmpNS );
        auto lock = nullspace->getLock();
        bool storeNS = true;
        
        if( !nullspace->verify(logger, entries,nParameters,nConstraints) ) {
            LOG_WARN << "NullSpace verification failed: possible hash-collision. Using local instance." << ende;
            LOG_DEBUG << printArray(entries,"\ngroupEntries") << printArray(nullspace->c_entries,"\nnsEntries") << ende;
            nullspace = tmpNS;
            storeNS = false;    // storing will overwrite the group we collided with.
            nullspace->ns.clear();
        }
        if( nullspace->ns.nDimensions() < 1 ) {
            nullspace->calculateNullspace( logger, storeNS, old_ns );
        }
        
        if( nullspace->ns_entries.empty() ) {
            LOG << "Mapping (" << nParameters << "x" << (nParameters-nConstraints) << ") nullspace." << ende;
            nullspace->mapNullspace();
        }

        for( const auto& entry: nullspace->ns_entries ) {
            ns_entries.insert(make_pair(entry.first + groupOffset,entry.second));
        }
        
    } else if (nParameters == nConstraints) {
        for(int32_t i=0; i<nParameters; ++i) ns_entries.insert(make_pair(PointI(groupOffset.y+i,-1),0.0));
    }
 
}


Constraints::Constraints( MomfbdJob& j ) : type(CT_UNDEF), job( j ), logger(j.logger),
    blockified(false), loaded(false), nParameters(0), nFreeParameters(0) {

}


Constraints::~Constraints() {

}


void Constraints::blockifyGroups( void ) {
    int32_t cOffset = 0;
    int32_t rOffset = 0;
    const int32_t* colMap = parameterOrder.get();
    std::set<int32_t> unused1(boost::counting_iterator<int32_t>(0), boost::counting_iterator<int32_t>(nParameters));            // Initializes set with values 0, ... , (nParameters-1)
    std::set<int32_t> unused2(boost::counting_iterator<int32_t>(0), boost::counting_iterator<int32_t>(nFreeParameters));
    for( auto & g : groups ) {
        g.blockify( parameterOrder.get(), cOffset, rOffset );
        for( auto& entry: g.ns_entries ) {
            auto index = entry.first;
            index.y = colMap[index.y];     // map parameter-index back to original order, so ns_entries can be applied directly to alpha.
            unused1.erase(index.y);
            unused2.erase(index.x);
            if( index.x < 0 ) {     // if this group has nParams = nConstriants just map it to 0 to prevent out-of-bounds, value is 0 anyway.
                index.x = 0;
            }
            ns_entries.insert(make_pair(index,entry.second));
        }
    }

    // These should always be the same size, corresponding to the number of "unconstrained" parameters (i.e. image-numbers existing only in 1 channel)
    if( unused1.size() == unused2.size() ) {
        if(unused1.size() > 0 ) {
            auto it1 = unused1.begin();
            for( auto & it2: unused2) {
                PointI index(*it1++,it2);
                ns_entries.insert(make_pair(index,1.0));
            }
        }
    } else {
        std::vector <int32_t> output;
        std::copy(unused1.begin(), unused1.end(), std::back_inserter(output));
        redux::file::Ana::write("unused1",output);
        output.clear();
        std::copy(unused2.begin(), unused2.end(), std::back_inserter(output));
        redux::file::Ana::write("unused2",output);
        LOG_ERR << "Size mismatch in Constraints::blockifyGroups(). This should *never* happen, so go looking for the bug. :-P";
    }
    
    blockified = true;
}


void Constraints::groupConnectedVariables( void ) {

    auto tmpC = constraints;
    while( tmpC.size() ) {
        Group g( logger, tmpC[0], (job.runFlags&RF_OLD_NS)  );
        g.addConnectedConstraints( tmpC );  // adds connected constraints to group and removes them from tmpC
        g.sortRows();
        groups.push_back( g );
    }

}


void Constraints::sortConstraints( bool blockwise ) {
    if(blockwise) {
        for( auto & g : groups ) {
            g.sortRows();
        }
    } else {
        sort( constraints.begin(), constraints.end(),
                []( const shared_ptr<Constraint> &lhs, const shared_ptr<Constraint> &rhs ) {
                    return ( ( *lhs ) > ( *rhs ) );
                }
            );
    }
}


uint64_t Constraints::size( void ) const {
    uint64_t sz = sizeof(nParameters) + sizeof(nFreeParameters);
    sz += sizeof(uint64_t) + ns_entries.size()*(PointI::size()+sizeof(double));
    return sz;
}


uint64_t Constraints::pack( char* ptr ) const {
    using redux::util::pack;
    uint64_t count = pack( ptr, nParameters );
    count += pack( ptr+count, nFreeParameters );
    count += pack( ptr+count, (uint64_t)ns_entries.size() );
    for(auto& entry: ns_entries) {
        count += entry.first.pack( ptr+count );
        count += pack( ptr+count, entry.second );
    }
    return count;
}


uint64_t Constraints::unpack( const char* ptr, bool swap_endian ) {
    using redux::util::unpack;
    uint64_t count = unpack( ptr, nParameters, swap_endian );
    count += unpack( ptr+count, nFreeParameters, swap_endian );
    uint64_t nEntries;
    count += unpack( ptr+count, nEntries, swap_endian );
    while(nEntries--) {
        pair<PointI,double> tmp;
        count += tmp.first.unpack( ptr+count, swap_endian);
        count += unpack( ptr+count, tmp.second );
        ns_entries.insert(tmp);
    }
    return count;
}


void Constraints::init( void ) {
    
    int32_t nModes = job.nModes;
    nParameters = job.nImages() * nModes;

    parameterOrder.reset( new int32_t[nParameters] );
    std::iota( parameterOrder.get(), parameterOrder.get()+nParameters, 0 );

    auto &objects = job.getObjects();
    if( nParameters && objects.size() ) {
        
        for( unsigned int modeIndex = 0; modeIndex < job.nModes; ++modeIndex ) {
            int32_t imageCount = 0;
            if( modeIndex == 0 || modeIndex == 1 ) {      // tilts
                shared_ptr<Constraint> thisConstraint( new Constraint( logger, imageCount * nModes + modeIndex, 1 ) );
                for( size_t count=objects[0]->getChannels()[0]->nImages(); count; --count ) {  // only for first object and channel
                    thisConstraint->addEntry( imageCount * nModes + modeIndex, 1.0 );
                    imageCount++;
                }
                constraints.push_back( thisConstraint );

                //int32_t kOffset = 0;
                map<int32_t, Constraint> wfCons;
                int32_t parameterOffset(0);
                for( auto& obj : objects ) {            // \alpha_{tkm} - \alpha_{t1m} - \alpha_{1km} + \alpha_{11m} = 0
                    for( auto& ch : obj->getChannels() ) {
                        for( auto& wf : ch->waveFrontList ) {
                            auto ret = wfCons.insert( make_pair( wf, Constraint( logger, parameterOffset + modeIndex, 1 ) ) );
                            if( !ret.second ) { // wavefront already existed in wfCons => constrain this image/mode coefficient
                                thisConstraint.reset( new Constraint( logger, ret.first->second ) );
                                thisConstraint->addEntry( parameterOffset + modeIndex, -1 );
                                constraints.push_back( thisConstraint );
                            }
                            parameterOffset += nModes;
                        }
//                         if( kOffset ) { // skip reference channel
//                             int32_t tOffset = 0;
//                             for( size_t count=ch->nImages(); count; --count ) {
// //                                 if( tOffset ) { // skip reference wavefront      TODO: select type with cfg option.
// //                                                                                  //TBD:  proper dealing with fixed differences.
// //                                     thisConstraint.reset( new Constraint( logger, modeIndex, 1 ) );          // \alpha_{11m}
// //                                     thisConstraint->addEntry( tOffset + modeIndex, -1 );                     // \alpha_{t1m}
// //                                     thisConstraint->addEntry( kOffset + tOffset + modeIndex, 1 );            // \alpha_{tkm}
// //                                     thisConstraint->addEntry( kOffset + modeIndex, -1 );                     // \alpha_{1km}
// //                                     constraints.push_back( thisConstraint );
// //                                 }
//                                 //if( tOffset ) { // skip reference wavefront
//                                    thisConstraint.reset( new Constraint( logger, tOffset + modeIndex, 1 ) ); // \alpha_{tkm}
//                                     thisConstraint->addEntry( kOffset + tOffset + modeIndex, -1 );          // \alpha_{t(k-1)m}
//                                     constraints.push_back( thisConstraint );
//                                 //}
//                                 tOffset += nModes;
//                             }
//                         }
//                         kOffset += ch->nImages() * nModes;
                    }
                }

            } else {                // for non-tilts, the mode coefficients are the same for co-temporal images (same wavefront)
                map<int32_t, Constraint> wfCons;
                for( auto& obj : objects ) {
                    for( auto& ch : obj->getChannels() ) {
                        for( auto& wf : ch->waveFrontList ) {
                            int32_t parameterOffset = imageCount * nModes;
                            auto ret = wfCons.insert( make_pair( wf, Constraint( logger, parameterOffset + modeIndex, 1 ) ) );
                            if( !ret.second ) { // wavefront already existed in wfCons => constrain this image/mode coefficient
                                shared_ptr<Constraint> thisConstraint( new Constraint( logger, ret.first->second ) );
                                thisConstraint->addEntry( parameterOffset + modeIndex, -1 );
                                constraints.push_back( thisConstraint );
                            }
                            imageCount++;
                        }
                    }
                }
            }
        }
        
        nFreeParameters = nParameters - constraints.size();

        groupConnectedVariables();
        
        blockifyGroups();
        
        loaded = true;

    }
}


Array<int16_t> Constraints::getMatrix( bool blocked ) const {

    int32_t nConstraints = constraints.size();
    if( nConstraints == 0 || nParameters == 0 ) {
        return Array<int16_t>();
    }

    Array<int16_t> c( nConstraints, nParameters );
    c.zero();
    const int32_t* colMap = parameterOrder.get();
    if( !blocked || groups.empty() ) {
        for( int i = 0; i < nConstraints; ++i ) {
            for( auto& entry : constraints[i]->entries ) {
                c( i, colMap[entry.first] ) = entry.second;
            }
        }
    } else {
        int32_t rowEnd = 0;
        int32_t colEnd = 0;
        //int cnt = 0;
        for( auto& group : groups ) {
            if( ! group.dense() ) {
                // TODO: Silent for now, but should throw or warn.
            }
            int32_t nConstr = group.constraints.size();
            int32_t nPar = group.indices.size();
            int32_t firstIndex = *group.indices.begin();
            int32_t rowBegin = rowEnd;
            int32_t colBegin = colEnd;
            rowEnd = rowBegin + nConstr - 1;
            colEnd = colBegin + nPar - 1;
            Array<int16_t> cc( c, rowBegin, rowEnd, colBegin, colEnd );  // sub-array of c
            for( int i = 0; i < nConstr; ++i ) {
                for( auto& entry : group.constraints[i]->entries ) {
                    if( entry.first - firstIndex >= nPar ) {
                        // TODO: Skip silently for now, but should throw or warn.
                        continue;
                    }
                    cc( i, entry.first - firstIndex ) = entry.second;
                }
            }
            rowEnd++;
            colEnd++;
        }
    }

    return c;

}


Array<double> Constraints::getNullMatrix( void ) const {

    if( nFreeParameters == 0 || nParameters == 0 ) {
        return Array<double>();
    }

    Array<double> ns( nParameters, nFreeParameters );
    ns.zero();
    
    for( auto& entry: ns_entries ) {
        ns(entry.first.y,entry.first.x) = entry.second;
    }
    
    return ns;

}


Array<int16_t> Constraints::getReducedMatrix( bool blocked ) const {

    int32_t nConstraints = constraints.size();
    if( nConstraints == 0 || nParameters == 0 ) {
        return Array<int16_t>();
    }

    Array<int16_t> c( nConstraints, nParameters );
    c.zero();
    const int32_t* colMap = parameterOrder.get();
    if( !blocked || groups.empty() ) {
        for( int i = 0; i < nConstraints; ++i ) {
            for( auto& entry : constraints[i]->entries ) {
                c( i, colMap[entry.first] ) = entry.second;
            }
        }
    } else {
        //int cnt = 0;
        int32_t col(0);
        int32_t row(0);
        std::set<size_t> hashes;
        for( auto& group : groups ) {
            auto ret = hashes.insert( group.entriesHash );
            if( ret.second == false ) continue;
            row += group.constraints.size();
            col += group.indices.size();
        }

        c.resize( row, col );
        c.zero();
        col = row = 0;
        hashes.clear();
        for( auto& group : groups ) {
            if( ! group.dense() ) {
                // TODO: Silent for now, but should throw or warn.
            }
            auto ret = hashes.insert( group.entriesHash );
            if( ret.second == false ) continue;
            int32_t nConstr = group.constraints.size();
            int32_t nPar = group.indices.size();
            int32_t firstIndex = *group.indices.begin();
            int32_t rowBegin = row;
            int32_t colBegin = col;
            row += nConstr;
            col += nPar;
            Array<int16_t> cc( c, rowBegin, row-1, colBegin, col-1 );  // sub-array of c
            for( int i = 0; i < nConstr; ++i ) {
                for( auto& entry : group.constraints[i]->entries ) {
                    if( entry.first - firstIndex >= nPar ) {
                        // TODO: Skip silently for now, but should throw or warn.
                        continue;
                    }
                    cc( i, entry.first - firstIndex ) = entry.second;
                }
            }
        }
    }

    return c;

}


Array<double> Constraints::getReducedNullMatrix( void ) const {

    if( nFreeParameters == 0 || nParameters == 0 ) {
        return Array<double>();
    }

    int32_t col(0);
    int32_t row(0);
    std::set<size_t> hashes;
    for( auto& group : groups ) {
        if( group.nullspace ) {
            auto ret = hashes.insert( group.nullspace->entriesHash );
            if( ret.second == false ) continue;
            row += group.nullspace->nParameters;
            col += group.nullspace->nParameters - group.nullspace->nConstraints;
        }
    }
    Array<double> ns( row, col );
    ns.zero();
    col = row = 0;
    hashes.clear();
    for( auto& group : groups ) {
        if( group.nullspace ) {
            auto ret = hashes.insert( group.nullspace->entriesHash );
            if( ret.second == false ) continue;
            int32_t nConstr = group.nullspace->nConstraints;
            int32_t nPar = group.nullspace->nParameters;
            int32_t rowBegin = row;
            int32_t colBegin = col;
            row += nPar;
            col += nPar-nConstr;
            Array<double> sns( ns, rowBegin, row-1, colBegin, col-1 );  // sub-array of c
            sns.copyFrom<double>( group.nullspace->ns.get() );
        }
    }
    
    return ns;

}


Array<int16_t> Constraints::getSubMatrix( int32_t groupID ) const {

    if( groupID >= (int32_t)groups.size() ) {
        return Array<int16_t>();
    }

    const Group& g = groups[groupID];
    if( ! g.dense() ) {
        return Array<int16_t>();
    }

    int32_t nConstr = g.constraints.size();
    int32_t nPar = g.indices.size();
    int32_t first = *g.indices.begin();
    Array<int16_t> c( nConstr, nPar );
    c.zero();
    for( int i = 0; i < nConstr; ++i ) {
        for( auto & entry : g.constraints[i]->entries ) {
            c( i, entry.first - first ) = entry.second;
        }
    }

    return c;

}


void Constraints::makeRowsCols (void) {

    for( auto& entry: ns_entries ) {
        ns_rows[entry.first.y].push_back( std::pair<size_t,double>(entry.first.x, entry.second) ); 
        ns_cols[entry.first.x].push_back( std::pair<size_t,double>(entry.first.y, entry.second) ); 
    }
//     int cnt(0);
// for ( auto& r: ns_cols ) {
//     if( cnt++ == 5 ) return;
//     cout << r.first << ":" << endl;
//     for(auto& e: r.second ) {
//         cout << "   { " << e.first << " -> " << e.second << " }" << endl;
//     }
// }
//cout << printArray( ns_cols,"cols" ) << endl;

}


void Constraints::read( void ) {

}


void Constraints::write( void ) {

}


void Constraints::dump( const string& tag ) const {

    LOG_DEBUG << "Dumping constraints." << ende;
    if( parameterOrder ) { // if we have an ordering-list (only when the constraints were generated, not after sending or loading)
        Ana::write( tag + "_order.f0", parameterOrder.get(), nParameters );
        Ana::write( tag + ".f0", getMatrix(false) );
        Ana::write( tag + "_ordered.f0", getMatrix(true) );
        Ana::write( tag + "_reduced.f0", getReducedMatrix(true) );
    }
    Ana::write( tag + "_ns.f0", getNullMatrix() );
    Ana::write( tag + "_reduced_ns.f0", getReducedNullMatrix() );

}


