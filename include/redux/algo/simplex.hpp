#ifndef REDUX_ALGO_SIMPLEX_HPP
#define REDUX_ALGO_SIMPLEX_HPP

//#include <Eigen/Dense>
#include <eigen3/Eigen/Dense>

namespace redux {

    namespace algo {
        
        /*!
         * 
         * @class Simplex
         *
         *  \brief Implementation of the Simplex Algorithm for solving linear programs
         * 
         *  Given an M-by-N matrix A, an M-length vector b, and an
         *  N-length vector c, solve the  LP { max cx : Ax <= b, x >= 0 }.
         *  Assumes that b >= 0 so that x = 0 is a basic feasible solution.
         *
         *  Creates an (M+1)-by-(N+M+1) simplex tableaux with the 
         *  RHS in column M+N, the objective function in row M, and
         *  slack variables in columns M through M+N-1.
         * 
         */
        class Simplex {
          
        public:
            Simplex(Eigen::MatrixXd& a, Eigen::VectorXd& b, Eigen::VectorXd& c);
            
            double value(void);
            
        private:
            
            void pivot(int p, int q);
            void solve(void);
            int bland(void);
            int dantzig(void);
            int minRatioRule(int q);

          Eigen::MatrixXd m_Tableu;
          Eigen::VectorXd m_b;
          Eigen::VectorXd m_c;
          
          int m_nConstraints;
          int m_nVariables;
            
        };
        
    }
    
}


#endif  // REDUX_ALGO_SIMPLEX_HPP
