/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_SparsifiedMatrix_h_
#define model_SparsifiedMatrix_h_

#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <vector>
#include <algorithm>
#include <cmath>


namespace model
{
    
class SparsifiedMatrix : public Eigen::SparseMatrix<double>
{
public:
    using Scalar = double;

    template<typename MatrixType, unsigned int UpLo>
    SparsifiedMatrix(
        const Eigen::SparseSelfAdjointView<MatrixType, UpLo>& A_,
        double diagonalCompensation,
        int nonDiagValues)
    {// Construct from symmetric dense matrix
        const MatrixType& A = A_.derived();

        const int n = static_cast<int>(A.rows());

        this->resize(n, n);

        Eigen::VectorXd diag = A.diagonal();

        // rowsToInsert[j] contains row indices to insert in column j
        std::vector<std::vector<int>> rowsToInsert(n);

        // Pass 1: each row selects its top-k off-diagonal neighbors.
        // For every selected symmetric edge, store both directions.
        for (int i = 0; i < n; ++i)
        {
            std::vector<std::pair<double, int>> rowValues;
            rowValues.reserve(n - 1);

            for (int j = 0; j < n; ++j)
            {
                if (j != i)
                {
                    rowValues.emplace_back(std::abs(A.coeff(i, j)), j);
                }
            }

            const int nkeep =
                std::min<int>(nonDiagValues, rowValues.size());

            std::partial_sort(
                rowValues.begin(),
                rowValues.begin() + nkeep,
                rowValues.end(),
                [](const auto& a, const auto& b)
                {
                    return a.first > b.first;
                }
            );

            for (int k = 0; k < nkeep; ++k)
            {
                const int j = rowValues[k].second;

                if (j != i)
                {
                    // Store full symmetric matrix pattern.
                    rowsToInsert[j].push_back(i);
                    rowsToInsert[i].push_back(j);
                }
            }
        }

        // Sort and remove duplicate row indices in each column.
        for (auto& rows : rowsToInsert)
        {
            std::sort(rows.begin(), rows.end());
            rows.erase(std::unique(rows.begin(), rows.end()), rows.end());
        }

        // Pass 2: diagonal compensation for dropped symmetric pairs.
        // Only scan lower triangle once.
        for (int j = 0; j < n; ++j)
        {
            const auto& rows = rowsToInsert[j];

            for (int i = j + 1; i < n; ++i)
            {
                const bool kept =
                    std::binary_search(rows.begin(), rows.end(), i);

                if (!kept)
                {
                    const double absVal = std::abs(A.coeff(i, j));

                    diag(i) += 2.0 * diagonalCompensation * absVal;
                    diag(j) += 2.0 * diagonalCompensation * absVal;
                }
            }
        }

        // Reserve exact number of entries per column: off-diagonal rows + diagonal.
        Eigen::VectorXi reservePerCol(n);

        for (int j = 0; j < n; ++j)
        {
            reservePerCol(j) =
                static_cast<int>(rowsToInsert[j].size()) + 1;
        }

        this->reserve(reservePerCol);

        // Pass 3: insert column-by-column in increasing row order.
        for (int j = 0; j < n; ++j)
        {
            const auto& rows = rowsToInsert[j];

            bool diagInserted = false;

            for (const int i : rows)
            {
                if (i != j)
                {
                    if (!diagInserted && j < i)
                    {
                        this->insert(j, j) = diag(j);
                        diagInserted = true;
                    }

                    this->insert(i, j) = A.coeff(i, j);

                }
            }

            if (!diagInserted)
            {
                this->insert(j, j) = diag(j);
            }
        }

        this->makeCompressed();
    }
    
//     template<typename TripletContainer>
//     SparsifiedMatrix(
//         int n,
//         const TripletContainer& triplets,
//         double diagonalCompensation,
//         int nonDiagValues)
//     {/*! Symmetric triplet constructor.
//          Input triplets must contain diagonal entries and at most one of (i,j), (j,i)
//          for each off-diagonal symmetric pair. */
        
//         this->resize(n, n);

//         Eigen::VectorXd diag = Eigen::VectorXd::Zero(n);
//         std::vector<std::vector<std::pair<double,int>>> rowValues(n);

//         for (const auto& t : triplets)
//         {
//             const int i = t.row();
//             const int j = t.col();
//             const double v = t.value();

//             if (i == j)
//             {
//                 diag(i) += v;
//             }
//             else
//             {
//                 rowValues[i].emplace_back(std::abs(v), j);
//                 rowValues[j].emplace_back(std::abs(v), i);
//             }
//         }

//         std::vector<std::vector<int>> rowsToInsert(n);

//         // Pass 1: each row votes for its top-k neighbors
//         for (int i = 0; i < n; ++i)
//         {
//             auto& row = rowValues[i];

//             const int nkeep =
//                 std::min<int>(nonDiagValues, row.size());

//             // std::cout<< "KEEPING "<<nkeep<<" OUT OF "<<row.size()<<" ENTRIES IN ROW "<<i<<std::endl;

//             std::partial_sort(
//                 row.begin(),
//                 row.begin() + nkeep,
//                 row.end(),
//                 [](const auto& a, const auto& b)
//                 {
//                     return a.first > b.first;
//                 });

//             for (int k = 0; k < nkeep; ++k)
//             {
//                 const int j = row[k].second;

//                 rowsToInsert[j].push_back(i);
//                 rowsToInsert[i].push_back(j);
//             }
//         }

//         for (auto& rows : rowsToInsert)
//         {
//             std::sort(rows.begin(), rows.end());
//             rows.erase(std::unique(rows.begin(), rows.end()), rows.end());
//         }

//         // Need values accessible by edge
//         std::map<std::pair<int,int>, double> values;

//         for (const auto& t : triplets)
//         {
//             const int i = t.row();
//             const int j = t.col();

//             if (i != j)
//             {
//                 const int a = std::min(i,j);
//                 const int b = std::max(i,j);
//                 values[{a,b}] += t.value();
//             }
//         }

//         // Diagonal compensation
//         for (const auto& [edge, value] : values)
//         {
//             const int i = edge.first;
//             const int j = edge.second;

//             const bool kept =
//                 std::binary_search(
//                     rowsToInsert[j].begin(),
//                     rowsToInsert[j].end(),
//                     i);

//             if (!kept)
//             {
//                 const double absVal = std::abs(value);
//                 diag(i) += 2.0 * diagonalCompensation * absVal;
//                 diag(j) += 2.0 * diagonalCompensation * absVal;
// //                diag(i) += 2.0 * diagonalCompensation * value;
// //                diag(j) += 2.0 * diagonalCompensation * value;

//             }
//         }

//         Eigen::VectorXi reservePerCol(n);

//         for (int j = 0; j < n; ++j)
//             reservePerCol(j) = static_cast<int>(rowsToInsert[j].size()) + 1;

//         this->reserve(reservePerCol);

//         for (int j = 0; j < n; ++j)
//         {
//             bool diagInserted = false;

//             for (const int i : rowsToInsert[j])
//             {
//                 if (!diagInserted && j < i)
//                 {
//                     this->insert(j,j) = diag(j);
//                     diagInserted = true;
//                 }

//                 const int a = std::min(i,j);
//                 const int b = std::max(i,j);

//                 this->insert(i,j) = values[{a,b}];
//             }

//             if (!diagInserted)
//                 this->insert(j,j) = diag(j);
//         }

//         this->makeCompressed();
//     }
    
    template<typename Derived>
    SparsifiedMatrix(
        const Eigen::DenseBase<Derived>& A,
        double diagonalCompensation,
        int nonDiagValues)
    {// Construct from non-symmetric dense matrix
        const int rows = static_cast<int>(A.rows());
        const int cols = static_cast<int>(A.cols());

        this->resize(rows, cols);

        Eigen::VectorXd diag =
            Eigen::VectorXd::Zero(std::min(rows, cols));

        for (int i = 0; i < diag.size(); ++i)
            diag(i) = A(i,i);

        std::vector<std::vector<int>> rowsToInsert(cols);

        // Row-wise selection
        for (int i = 0; i < rows; ++i)
        {
            std::vector<std::pair<double,int>> rowValues;
//            rowValues.reserve(cols - 1);
            rowValues.reserve(cols > 0 ? cols - 1 : 0);

            for (int j = 0; j < cols; ++j)
            {
                if (j != i)
                {
                    rowValues.emplace_back(std::abs(A(i,j)), j);
                }
            }

            const int nkeep =
                std::min<int>(nonDiagValues, rowValues.size());

            std::partial_sort(
                rowValues.begin(),
                rowValues.begin() + nkeep,
                rowValues.end(),
                [](const auto& a, const auto& b)
                {
                    return a.first > b.first;
                });

            std::vector<bool> keep(cols,false);

            for (int k=0; k<nkeep; ++k)
            {
                const int j = rowValues[k].second;
                keep[j] = true;

                rowsToInsert[j].push_back(i);
            }

            for (int j=0; j<cols; ++j)
            {
                if (!(j==i || keep[j]))
                {
                    if (i < diag.size())
                    {
                        diag(i) += diagonalCompensation * std::abs(A(i,j));
                    }
                }
            }
        }

        for(auto& rows : rowsToInsert)
            std::sort(rows.begin(), rows.end());

        Eigen::VectorXi reserve(cols);

        for(int j=0; j<cols; ++j)
        {
            reserve(j) =
                static_cast<int>(rowsToInsert[j].size())
                + (j < diag.size());
        }

        this->reserve(reserve);

        for(int j=0; j<cols; ++j)
        {
            bool diagInserted = false;

            for(int i : rowsToInsert[j])
            {
                if(j < diag.size() && !diagInserted && j < i)
                {
                    this->insert(j,j) = diag(j);
                    diagInserted = true;
                }

                this->insert(i,j) = A(i,j);
            }

            if(j < diag.size() && !diagInserted)
            {
                this->insert(j,j) = diag(j);
            }
        }

        this->makeCompressed();
    }
    
    template<typename SparseDerived>
    SparsifiedMatrix(
        const Eigen::SparseMatrixBase<SparseDerived>& A_,
        double diagonalCompensation,
        int nonDiagValues)
    {// Construct from non-symmetric sparse matrix
        const SparseDerived& A = A_.derived();

        const int rows = static_cast<int>(A.rows());
        const int cols = static_cast<int>(A.cols());

        this->resize(rows, cols);

        Eigen::VectorXd diag =
            Eigen::VectorXd::Zero(std::min(rows, cols));

        for(int i=0; i<diag.size(); ++i)
            diag(i) = A.coeff(i,i);

        std::vector<std::vector<int>> rowsToInsert(cols);

        // Build row structure
        std::vector<std::vector<std::pair<double,int>>> rowValues(rows);

        for(int j=0; j<cols; ++j)
        {
            for(typename SparseDerived::InnerIterator it(A,j); it; ++it)
            {
                const int i = it.row();

                if(i != j)
                {
                    rowValues[i].emplace_back(
                        std::abs(it.value()),
                        j);
                }
            }
        }

        for(int i=0; i<rows; ++i)
        {
            auto& row = rowValues[i];

            const int nkeep =
                std::min<int>(nonDiagValues,row.size());

            std::partial_sort(
                row.begin(),
                row.begin()+nkeep,
                row.end(),
                [](const auto& a,const auto& b)
                {
                    return a.first > b.first;
                });

//            std::vector<bool> keep(row.size(),false);

            for(int k=0;k<nkeep;++k)
            {
                const int j = row[k].second;

                rowsToInsert[j].push_back(i);
  //              keep[k]=true;
            }

            for(int k=nkeep;k<int(row.size());++k)
            {
                if(i < diag.size())
                {
                    diag(i)+=
                    diagonalCompensation*row[k].first;
                }
            }
        }

        for(auto& rows : rowsToInsert)
            std::sort(rows.begin(),rows.end());

        Eigen::VectorXi reserve(cols);

        for(int j=0;j<cols;++j)
        {
            reserve(j)=
                static_cast<int>(rowsToInsert[j].size())
                +(j<diag.size());
        }

        this->reserve(reserve);

        for(int j=0;j<cols;++j)
        {
            bool diagInserted=false;

            for(int i : rowsToInsert[j])
            {
                if(j<diag.size() && !diagInserted && j<i)
                {
                    this->insert(j,j)=diag(j);
                    diagInserted=true;
                }

                this->insert(i,j)=A.coeff(i,j);
            }

            if(j<diag.size() && !diagInserted)
            {
                this->insert(j,j)=diag(j);
            }
        }

        this->makeCompressed();
    }
};
    
    
} // close namespace model
#endif
