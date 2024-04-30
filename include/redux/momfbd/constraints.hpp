#ifndef REDUX_MOMFBD_CONSTRAINTS_HPP
#define REDUX_MOMFBD_CONSTRAINTS_HPP

#include "redux/util/stringutil.hpp"
#include "redux/util/array.hpp"
#include "redux/util/cacheitem.hpp"
#include "redux/util/point.hpp"
#include "redux/util/trace.hpp"

#include <string>
#include <map>
#include <memory>
#include <set>
#include <vector>

#include <gsl/gsl_permutation.h>

namespace redux {

    namespace logging {
        class Logger;
    }
    
    namespace momfbd {

        /*! @ingroup momfbd
         *  @{
         */

        class MomfbdJob;

        /*! @brief Container for Linear Equality Constraints.
         *
         */
        struct Constraints
#ifdef RDX_TRACE_PROC
            : public redux::util::TraceObject<Constraints>
#endif
        {

            enum ConstraintType { CT_UNDEF = 0, CT_CALIB, CT_OLD, CT_NEW };
            const std::string ConstraintTag[4] = {"", "CALIB", "MOJPDSF", "MOMFBD"};

            /*!
             * Contains one constraint (= one row in the constraint matrix)
             */
            struct Constraint {
                Constraint (logging::Logger& l, const Constraint& rhs) : entries (rhs.entries),logger(l) {};
                Constraint (Constraint && rhs) : entries (std::move (rhs.entries)),logger(rhs.logger) {};
                Constraint (logging::Logger& l, int32_t i, int8_t v) : logger(l)  { entries[i] = v; }
                void addEntry (int32_t i, int8_t v) { entries[i] = v; }
                int32_t firstEntry (void) { return entries.begin()->first; }
                int8_t entry (int32_t ind) { return (entries.find (ind) != entries.end()) ? entries[ind] : 0.0; }
                bool has (int32_t ind) { return (entries.find (ind) != entries.end()); }
                void replaceIndices (const std::map<int32_t, int32_t>& indexMap);
                bool operator> (const Constraint& rhs) const;
                std::map<int32_t, int8_t> entries;      //!< Index and value of non-zero matrix elements
                logging::Logger& logger;
            };
            
            
            struct NullSpace : public redux::util::CacheItem {
                
                NullSpace(const std::map<int32_t, int8_t>& e, int32_t np, int32_t nc, bool old_ns);
                
                void mapNullspace( void );
                void calculateNullspace( logging::Logger&, bool store=true, bool nc=false );
                bool verify( logging::Logger&,const std::map<int32_t, int8_t>&, int32_t, int32_t );
                std::unique_lock<std::mutex> getLock( bool trylock=false ) {
                    if(trylock) return std::unique_lock<std::mutex>(nsMtx,std::try_to_lock);
                    return std::unique_lock<std::mutex>(nsMtx);
                }
        
                
                size_t csize(void) const override;
                uint64_t cpack(char*) const override;
                uint64_t cunpack(const char*, bool) override;
                void cclear(void) override;
                bool operator<( const NullSpace& ) const;

                std::map<int32_t, int8_t> c_entries;                       //!< Index and value of (non-zero) constraint-matrix elements
                size_t entriesHash;
                int32_t nParameters;
                int32_t nConstraints;
                redux::util::Array<double> ns;
                std::map<redux::util::PointI, double> ns_entries;  //!< (row,col) and value of (non-zero) nullspace-matrix elements
                std::mutex nsMtx;
            };

            /*!
             * Contains connected constraints
             */
            struct Group {
                
                explicit Group (logging::Logger& l, std::shared_ptr<Constraint>& con, bool old_ns=false);
                
                void add (const std::shared_ptr<Constraint>& con);
                void addConnectedConstraints (std::vector<std::shared_ptr<Constraint>>& cons);
                void blockify (int32_t* columnOrdering, int32_t& cOffset, int32_t& rOffset);
                void sortRows(void);
                bool dense (void) const;

                void mapNullspace(void);
                
                std::vector<std::shared_ptr<Constraint>> constraints;
                std::set<int32_t> indices;
                std::map<int32_t, int8_t> entries;                         //!< Index and value of (non-zero) constraint-matrix elements
                std::map<redux::util::PointI, double> ns_entries;          //!< (row,col) and value of (non-zero) nullspace-matrix elements
                int32_t nParameters;
                redux::util::PointI groupOffset;                           //!< Offset of this group in the nullspace matrix (y=row,x=col)
                size_t entriesHash;
                std::shared_ptr<NullSpace> nullspace;
                logging::Logger& logger;
                bool old_ns;
            };

            explicit Constraints(MomfbdJob&);
            ~Constraints();

            void blockifyGroups (void);
            void groupConnectedVariables (void);
            void sortConstraints (bool blockwise=false);

            uint64_t size (void) const;
            uint64_t pack (char*) const;
            uint64_t unpack (const char*, bool);
            
            template <typename T> void apply(const T* in, T* out) const {
                memset(out,0,nFreeParameters*sizeof(T));
                for( auto& entry: ns_entries ) {
                    out[entry.first.x] += entry.second * in[entry.first.y];
                }
            }
            template <typename T> void reverse(const T* in, T* out) const {
                memset(out,0,nParameters*sizeof(T));
                for( auto& entry: ns_entries ) {
                    out[entry.first.y] += entry.second * in[entry.first.x];
                }
            }

            template <typename T> void reverseAndAdd(const T* in, const T* offset, T* out) const {
                memcpy( out, offset, nParameters*sizeof(T) );
                for( auto& entry: ns_entries ) {
                    out[entry.first.y] += entry.second * in[entry.first.x];
                }
            }

            template <typename T> void apply(const T* in, T* out, double scale) const {
                memset(out,0,nFreeParameters*sizeof(T));
                for( auto& entry: ns_entries ) {
                    out[entry.first.x] += entry.second * scale * in[entry.first.y];
                }
            }
            template <typename T> void reverse(const T* in, T* out, double scale) const {
                memset(out,0,nParameters*sizeof(T));
                for( auto& entry: ns_entries ) {
                    out[entry.first.y] += entry.second * scale * in[entry.first.x];
                }
            }

            int32_t nConstraints (void) const { return constraints.size(); };

            void init (void);
            redux::util::Array<int16_t> getMatrix (bool blocked=false) const;
            redux::util::Array<double> getNullMatrix (void) const;
            redux::util::Array<int16_t> getReducedMatrix (bool blocked=false) const;
            redux::util::Array<double> getReducedNullMatrix (void) const;
            redux::util::Array<int16_t> getSubMatrix (int32_t groupID) const;
            void makeRowsCols (void);
            void read (void);
            void write (void);
            inline bool verify (void) const { return loaded; }

            void dump( const std::string& tag="constraints" ) const;

            std::vector<std::shared_ptr<Constraint>> constraints;
            std::vector<Group> groups;

            ConstraintType type;

            const MomfbdJob& job;
            logging::Logger& logger;
            bool blockified;
            bool loaded;

            int32_t nParameters;
            int32_t nFreeParameters;
            std::unique_ptr<int32_t[]> parameterOrder;
            std::map<redux::util::PointI, double> ns_entries;        //!< (row,col) and value of (non-zero) nullspace-matrix elements
            std::map<size_t,std::vector<std::pair<size_t,double>>> ns_rows, ns_cols;

        };


        /*! @} */


    }

}

#endif  // REDUX_MOMFBD_CONSTRAINTS_HPP


