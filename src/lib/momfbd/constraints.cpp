#include "redux/momfbd/constraints.hpp"

#include "redux/logger.hpp"
#include "redux/momfbd/momfbdjob.hpp"
#include "redux/file/fileana.hpp"
#include "redux/math/linalg.hpp"
#include "redux/util/datautil.hpp"

#include <algorithm>
#include <iostream>

#include <boost/functional/hash.hpp>

using namespace redux::momfbd;
using namespace redux::file;
using namespace redux::math;
using namespace redux::util;
using namespace redux;

using namespace std;

#define lg Logger::mlg
namespace {

    const string thisChannel = "constraints";

}


void Constraints::Constraint::replaceIndices( const map<uint32_t, uint32_t>& indexMap ) {
    map<uint32_t, int8_t> newEntries;
    for( auto & it : entries ) {
        newEntries[indexMap.at( it.first )] = it.second;
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


Constraints::NullSpace::NullSpace(const std::map<uint32_t, int8_t>& e, uint32_t np, uint32_t nc) :
    c_entries(e), entriesHash(0), nParameters(np), nConstraints(nc) {
        
}


#define NS_THRESHOLD 1E-10
void Constraints::NullSpace::mapNullspace(void) {
    ns_entries.clear();
    int nCols = nParameters - nConstraints;
    LOG_DEBUG << "Mapping (" << nParameters << "x" << nCols << ") nullspace.";
    int count = 0;
    for( auto &it: ns ) {
        if( abs(it) > NS_THRESHOLD ) {
            ns_entries.insert(make_pair(Point16(count/nCols,count%nCols),it));
        }
        count++;
    }
    //static int cnt(0);
    //Ana::write( "nullmatrix" + to_string( cnt++ ) + ".f0", ns );
}
#undef NS_THRESHOLD


#define NS_ANALYTIC 0
void Constraints::NullSpace::calculateNullspace(bool store) {
    
    if(NS_ANALYTIC) {
        LOG_DETAIL << "Finding nullspace basis for (" << nConstraints << "x" << nParameters << ") group.";
        ns.resize(nParameters,nParameters-nConstraints);
        if (nParameters == nConstraints+1) {    // 1D nullspace
            ns = 1.0/sqrt(nParameters);
        } else {
            // TODO implement
        }
        
    } else {
        LOG_DETAIL << "Calculating QR-decomposition for (" << nConstraints << "x" << nParameters << ") group.";
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
    isLoaded = true;
    if(store) cacheStore();
    
}
#undef NS_ANALYTIC


bool Constraints::NullSpace::verify(const std::map<uint32_t, int8_t>& e, uint32_t nP, uint32_t nC) {
    
    if(nP != nParameters || nC != nConstraints || c_entries != e) return false;
    
    if( entriesHash == 0 ) {
        entriesHash = boost::hash_range(c_entries.begin(), c_entries.end());
        setPath("Nullspace_"+to_string(nP)+"x"+to_string(nC)+"_"+to_string(entriesHash));
    }
    
    if( !isLoaded )  {
        if( cacheLoad() ) {
            if(nP == nParameters && nC == nConstraints && c_entries == e) {
                LOG_DEBUG << "Found matching constraint-group in file: " << itemPath;
                return true;
            }
        } else {
            LOG_DEBUG << "No matching constraint-group found.";
            ns.resize(0);
            ns_entries.clear();
        }
    }
        
    return true;
}


size_t Constraints::NullSpace::csize(void) const {
    size_t sz = sizeof(entriesHash) + sizeof(nParameters) + sizeof(nConstraints);
    sz += c_entries.size()*(sizeof(uint32_t)+sizeof(int8_t)) + sizeof(uint64_t);
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


bool Constraints::NullSpace::operator<(const NullSpace& rhs) {
    if( entriesHash == rhs.entriesHash) {
        return (c_entries < rhs.c_entries);
    }
    return (entriesHash < rhs.entriesHash);
}


Constraints::Group::Group( shared_ptr<Constraint>& con ) : nParameters(0), entriesHash(0) {
    add( con );
}


void Constraints::Group::add( const shared_ptr<Constraint>& con ) {
    if( find( constraints.begin(), constraints.end(), con ) == constraints.end() ) { // don't add the same row twice
        constraints.push_back( con );
        for( auto it : con->entries ) {
            indices.insert( it.first );
        }
    }
}


void Constraints::Group::addConnectedConstraints( vector<shared_ptr<Constraint>>& cons ) {
    for( auto it = cons.begin(); it < cons.end(); ) {
        bool connected( false );
        for( auto entry : ( *it )->entries ) {  // if one of the column-indices matches this group, add the constraint to group.
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
void Constraints::Group::blockify( uint32_t* columnOrdering, uint32_t& rOffset, uint32_t& cOffset ) {

    groupOffset.y = rOffset;    // store the offset of this group (in the nullspace matrix).
    groupOffset.x = cOffset;
    
    set<uint32_t> tmpIndices;
    map<uint32_t, uint32_t> indexMap;
    for( auto ind : indices ) {
        tmpIndices.insert( rOffset );
        indexMap[ind] = rOffset;
        std::swap( columnOrdering[rOffset++], ind );
    }
    indices = std::move( tmpIndices );
    nParameters = indices.size();
    
    std::map<uint32_t, int8_t> tmpEntries;
    uint32_t indexOffset = 0;
    uint32_t firstIndex = *indices.begin();
    for( auto & it : constraints ) {
        it->replaceIndices( indexMap );
        for( auto & ind : it->entries ) {
            //tmpEntries.insert(make_pair(ind.first-indexOffset,ind.second));
            tmpEntries.insert(make_pair(indexOffset-firstIndex+ind.first,ind.second));
        }
        indexOffset += nParameters;     // next row
    }
    entries = std::move( tmpEntries );
    entriesHash = boost::hash_range(entries.begin(), entries.end());
    
    cOffset += (nParameters - constraints.size());
    
    mapNullspace();
    
    

 
    // Now that this group is "dense", sort columns by column-sum.
//     for( auto ind1: indices ) {
//         for( auto ind2: indices ) {
//             if ( ind1 != ind2) {
//                 float colSum1(0);
//                 float colSum2(0);
//                 for( auto &it: constraints ) {
//                     colSum1 += it->entry(ind1);
//                     colSum2 += it->entry(ind2);
//                 }
//                 if( colSum1 > colSum2 ) {
//                     for( auto &it: constraints ) {
//                         it->swap(ind1,ind2);
//                         //std::swap(columnOrdering[ind1],columnOrdering[ind2]);
//                     }
//                 }
//             }
//         }
//     }

    // ...and rows by number of first entry -> nEntries -> value
//     std::sort( constraints.begin(), constraints.end(),
//             []( const shared_ptr<Constraint> &lhs, const shared_ptr<Constraint> &rhs ) {
//                 return ( (*lhs) > (*rhs) );
//             }
//         );

//    cout << "SortE: " << offset << printArray(columnOrdering,"  currentOrdering") << printArray(indices,"  indices") << endl;
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

    uint32_t first = *indices.begin();
    uint32_t last = *indices.rbegin();
    return ( ( last - first + 1 ) == indices.size() );

}


void Constraints::Group::mapNullspace(void) {
    
    shared_ptr<NullSpace> tmpNS( new NullSpace(entries,nParameters,constraints.size()) );
    nullspace = redux::util::Cache::get( entriesHash, tmpNS );
    
    bool storeNS = true;
    if( !nullspace->verify(entries,nParameters,constraints.size()) ) {
        LOG_WARN << "NullSpace verification failed: possible hash-collision. Using local instance.";
        LOG_DEBUG << printArray(entries,"\ngroupEntries") << printArray(nullspace->c_entries,"\nnsEntries");
        nullspace = tmpNS;
        storeNS = false;    // storing will overwrite the group we collided with.
        nullspace->ns.resize(0);
    }
    
    if( nullspace->ns.nDimensions() < 1 ) {
        nullspace->calculateNullspace(storeNS);
    }
        
    if( nullspace->ns_entries.empty() ) {
        nullspace->mapNullspace();
    }
    
    ns_entries.clear();
    for( auto& it: nullspace->ns_entries ) {
        ns_entries.insert(make_pair(it.first + groupOffset,it.second));
    }


 
}


Constraints::Constraints( const MomfbdJob& j ) : job( j ), blockified(false) {

    init();

}


Constraints::~Constraints() {

}


void Constraints::blockifyGroups( void ) {
    uint32_t cOffset = 0;
    uint32_t rOffset = 0;
    const uint32_t* colMap = parameterOrder.get();
    for( auto & g : groups ) {
        g.blockify( parameterOrder.get(), cOffset, rOffset );
        for( auto& it: g.ns_entries ) {     // map indices back to original parameter order, so ns_entries can be applied directly to alpha.
            auto index = it.first;
            index.y = colMap[index.y];
            ns_entries.insert(make_pair(index,it.second));
        }
    }
    blockified = true;
}


void Constraints::groupConnectedVariables( void ) {

    auto tmpC = constraints;
    while( tmpC.size() ) {
        Group g( tmpC[0] );
        g.addConnectedConstraints( tmpC );  // adds connected constraints to group and removes them from tmpC
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
    uint64_t sz = sizeof(nParameters) + sizeof(nConstrainedParameters);
    sz += sizeof(uint64_t) + ns_entries.size()*(Point16::size()+sizeof(double));
    return sz;
}


uint64_t Constraints::pack( char* ptr ) const {
    using redux::util::pack;
    uint64_t count = pack( ptr, nParameters );
    count += pack( ptr+count, nConstrainedParameters );
    count += pack( ptr+count, (uint64_t)ns_entries.size() );
    for(auto& it: ns_entries) {
        count += it.first.pack( ptr+count );
        count += pack( ptr+count, it.second );
    }
    return count;
}


uint64_t Constraints::unpack( const char* ptr, bool swap_endian ) {
    using redux::util::unpack;
    uint64_t count = unpack( ptr, nParameters, swap_endian );
    count += unpack( ptr+count, nConstrainedParameters, swap_endian );
    uint64_t nEntries;
    count += unpack( ptr+count, nEntries, swap_endian );
    while(nEntries--) {
        pair<Point16,double> tmp;
        count += tmp.first.unpack( ptr+count, swap_endian);
        count += unpack( ptr+count, tmp.second );
        ns_entries.insert(tmp);
    }
    return count;
}


std::string Constraints::name( const string& base ) const {
    string ret = base;

    cout << "constraint file: \"" << ret << "\"" << endl;
    return ret;
}


void Constraints::init( void ) {
    
    uint32_t nModes = job.modeNumbers.size();
    nParameters = job.nImages() * nModes;

    parameterOrder.reset( new uint32_t[nParameters] );
    uint32_t n = 0;
    generate( parameterOrder.get(), parameterOrder.get()+nParameters, [&] { return n++; } );

    if( nParameters ) {
        
        //const vector<shared_ptr<Object>>& objects = job.getObjects();
        for( uint modeIndex = 0; modeIndex < job.modeNumbers.size(); ++modeIndex ) {
            uint16_t modeNumber = job.modeNumbers[modeIndex];
            uint32_t imageCount = 0;
            if( modeNumber == 2 || modeNumber == 3 ) {      // tilts
                shared_ptr<Constraint> c( new Constraint( imageCount * nModes + modeIndex, 1 ) );
                for( auto wf UNUSED: job.getObjects()[0]->getChannels()[0]->imageNumbers ) {  // only for first object and channel
                    c->addEntry( imageCount * nModes + modeIndex, 1.0 );
                    imageCount++;
                }
                constraints.push_back( c );

                map<uint32_t, Constraint> wfCons;
                uint32_t kOffset = 0;
                for( auto obj : job.getObjects() ) {            // \alpha_{tkm} - \alpha_{t1m} - \alpha_{1km} + \alpha_{11m} = 0
                    for( auto ch : obj->getChannels() ) {
                        if( kOffset ) { // skip reference channel
                            uint32_t tOffset = 0;
                            for( auto wf UNUSED: ch->imageNumbers ) {
                                if( tOffset ) { // skip reference wavefront
                                    shared_ptr<Constraint> c( new Constraint( modeIndex, 1 ) ); // \alpha_{11m}
                                    c->addEntry( tOffset + modeIndex, -1 );                     // \alpha_{t1m}
                                    c->addEntry( kOffset + tOffset + modeIndex, 1 );            // \alpha_{tkm}
                                    c->addEntry( kOffset + modeIndex, -1 );                     // \alpha_{1km}
                                    constraints.push_back( c );
                                }
                                tOffset += nModes;
                            }
                        }
                        kOffset += ch->nImages() * nModes;
                    }
                }

            } else {                // for non-tilts, the mode coefficients are the same for co-temporal images (same wavefront)
                map<uint32_t, Constraint> wfCons;
                for( auto obj : job.getObjects() ) {
                    for( auto ch : obj->getChannels() ) {
                        for( auto wf : ch->imageNumbers ) {
                            uint32_t parameterOffset = imageCount * nModes;
                            auto ret = wfCons.insert( make_pair( wf, Constraint( parameterOffset + modeIndex, 1 ) ) );
                            if( !ret.second ) { // wavefront already existed in wfCons => constrain this image/mode coefficient
                                shared_ptr<Constraint> c( new Constraint( ret.first->second ) );
                                c->addEntry( parameterOffset + modeIndex, -1 );
                                constraints.push_back( c );
                            }
                            imageCount++;
                        }
                    }
                }
            }
        }

        //cout << "Constraints::Constraints()  " << constraints.size() << "x" << nParameters << endl;
        
        nConstrainedParameters = nParameters - constraints.size();

        //dump( "constraints" );

        groupConnectedVariables();
        
        blockifyGroups();

        //dump( "blockified" );

    }
}


Array<int16_t> Constraints::getMatrix( void ) const {

    uint32_t nConstraints = constraints.size();
    if( nConstraints == 0 || nParameters == 0 ) {
        return Array<int16_t>();
    }

    Array<int16_t> c( nConstraints, nParameters );
    c.zero();
    if( !blockified || groups.empty() ) {
        for( uint i = 0; i < nConstraints; ++i ) {
            for( auto & entry : constraints[i]->entries ) {
                c( i, entry.first ) = entry.second;
            }
        }
    } else {
        uint32_t rowEnd = 0;
        uint32_t colEnd = 0;
        //int cnt = 0;
        for( auto &group : groups ) {
            if( ! group.dense() ) {
                // TODO: Silent for now, but should throw or warn.
            }
            uint32_t nConstr = group.constraints.size();
            uint32_t nPar = group.indices.size();
            uint32_t firstIndex = *group.indices.begin();
            uint32_t rowBegin = rowEnd;
            uint32_t colBegin = colEnd;
            rowEnd = rowBegin + nConstr - 1;
            colEnd = colBegin + nPar - 1;
            Array<int16_t> cc( c, rowBegin, rowEnd, colBegin, colEnd );  // sub-array of c
            for( uint i = 0; i < nConstr; ++i ) {
                //uint32_t col = 0;
                for( auto & entry : group.constraints[i]->entries ) {
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

    //uint32_t nConstraints = constraints.size();
    if( nConstrainedParameters == 0 || nParameters == 0 ) {
        return Array<double>();
    }

    Array<double> ns( nParameters, nConstrainedParameters );
    ns.zero();
    
    for( auto &it: ns_entries ) {
        ns(it.first.y,it.first.x) = it.second;
    }
    
    return ns;

}


Array<int16_t> Constraints::getSubMatrix( uint32_t groupID ) const {

    if( groupID >= groups.size() ) {
        return Array<int16_t>();
    }

    const Group& g = groups[groupID];
    if( ! g.dense() ) {
        return Array<int16_t>();
    }

    uint32_t nConstr = g.constraints.size();
    uint32_t nPar = g.indices.size();
    uint32_t first = *g.indices.begin();
    //uint32_t last = *g.indices.rbegin();
    //cout << "Constraints::getSubMatrix(" << groupID << ")  " << nConstr << "x" << nPar << "   first = " << first << "   last = " << last << endl;
    Array<int16_t> c( nConstr, nPar );
    c.zero();
    for( uint i = 0; i < nConstr; ++i ) {
        //  cout << "Constraints::getSubMatrix(" << groupID << ")  i = " << i << endl;
       // uint32_t col = 0;
        for( auto & entry : g.constraints[i]->entries ) {
            //      cout << "Constraints::getSubMatrix(" << groupID << ")  e = " << entry.first << endl;
            c( i, entry.first - first ) = entry.second;
        }
    }

    return c;

}


void Constraints::read( void ) {

}


void Constraints::write( void ) {

}


void Constraints::dump( string tag ) {

    cout << "Dumping constraints." << endl;
    cout << printArray( parameterOrder.get(), nParameters, "ordering" ) << endl;
    Ana::write( tag + ".f0", getMatrix() );
    Ana::write( tag + "_ns.f0", getNullMatrix() );
    for( uint i = 0; i < groups.size(); ++i ) {
        Ana::write( tag + "_sub_" + to_string( i ) + ".f0", getSubMatrix( i ) );
    }
}


