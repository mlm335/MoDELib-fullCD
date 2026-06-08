#ifndef model_SpatialODESolver_H_
#define model_SpatialODESolver_H_

#include<Eigen/SparseCore>
#include<Eigen/SparseLU>
#include<Eigen/IterativeLinearSolvers>
#include<TrialBase.h>

namespace model
{
    // template<typename TrialFunctionType>
    class SpatialODESolver
    {

        typedef Eigen::SparseMatrix<double,Eigen::RowMajor> SparseMatrixType;

        SparseMatrixType A;
        const bool use_directSolver;
        const double tolerance;
        const size_t gSize;  // Store gSize directly

        double lastError;
        int lastIterations;

// #ifdef CHOLMOD_H // SuiteSparse Cholmod (LLT) module
//         typedef Eigen::CholmodSimplicialLDLT<SparseMatrixType> DirectSolverType; //If problem with speed try this
// #else
//         typedef Eigen::SimplicialLLT<SparseMatrixType> DirectSolverType;
// #endif
        // typedef Eigen::ConjugateGradient<SparseMatrixType> IterativeSolverType;

// #ifdef UMFPACK_H // SuiteSparse Umfpack (LU) module
//         typedef Eigen::UmfPackLU<SparseMatrixType> DirectSolverType;
// #else
//         typedef Eigen::SparseLU<SparseMatrixType> DirectSolverType;
// #endif
//         typedef Eigen::BiCGSTAB<SparseMatrixType> IterativeSolverType;

        typedef Eigen::SimplicialLLT<SparseMatrixType> DirectSolverType;
        typedef Eigen::ConjugateGradient<SparseMatrixType> IterativeSolverType;

        DirectSolverType directSolver;
        IterativeSolverType iterativeSolver;

    public:
        
        template<typename TrialFunctionType, typename IntegrationDomainType>
        SpatialODESolver(const TrialFunctionType& trial, const IntegrationDomainType& dV, const bool& use_directSolver_in, const double& tol):
        /* init */use_directSolver(use_directSolver_in)
        /* init */,tolerance(tol)
        /* init */,gSize(trial.gSize())
        /* init */,lastError(0.0)
        /* init */,lastIterations(0)
        {                    
            // const size_t gSize=trial.gSize();
            const auto bwfAi((test(trial),trial)*dV);
            std::vector<Eigen::Triplet<double> > globalTripletsAi(bwfAi.globalTriplets());
            A.resize(gSize,gSize);
            A.setFromTriplets(globalTripletsAi.begin(),globalTripletsAi.end());
            
            if(use_directSolver)
            {
                directSolver.compute(A);
                assert(directSolver.info()==Eigen::Success && "SOLVER  FAILED");
            }
            else
            {
                iterativeSolver.compute(A);
                assert(iterativeSolver.info()==Eigen::Success && "SOLVER  FAILED");
            }
        }
        
        // Eigen::VectorXd solve(const TrialFunctionType& trial, const Eigen::VectorXd& b)
        Eigen::VectorXd solve(const Eigen::VectorXd& b)
        {/*!@param[in] b the rhs of A*x=b
          * @param[in] y the guess for x
          */            

            if(use_directSolver)
            {
                assert(directSolver.info()==Eigen::Success && "SOLVER  FAILED");
                Eigen::VectorXd x = directSolver.solve(b);
                lastError = 0.0;
                lastIterations = 1;
                return x;                
            }
            else
            {
                iterativeSolver.setTolerance(tolerance);
                // const size_t gSize=trial.gSize();
                Eigen::VectorXd g(Eigen::VectorXd::Zero(gSize));
                assert(iterativeSolver.info()==Eigen::Success && "SOLVER  FAILED");                
                Eigen::VectorXd x = iterativeSolver.solveWithGuess(b,g);
                lastError = iterativeSolver.error();
                lastIterations = iterativeSolver.iterations();
                return x;
            }
        }
        
        double error() const
        {
            return lastError;
        }

        int iterations() const
        {
            return lastIterations;
        }
    };

}
#endif
