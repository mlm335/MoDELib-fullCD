/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_MDStackingFaultNoise_cpp
#define model_MDStackingFaultNoise_cpp

#include <MDStackingFaultNoise.h>
#include <cmath>

namespace model
{
MDStackingFaultNoise::MDStackingFaultNoise(const PolycrystallineMaterialBase& mat,
                                           const std::string& tag,
                                           const std::string& correlationFile_in,
                                           const int& seed,
                                           const GridSizeType& gridSize,
                                           const GridSpacingType& gridSpacing,
                                           const Eigen::Matrix<double,2,2>& latticeBasis
                                           ) :
  /* init */ GlidePlaneNoiseBase<1>("MDStackingFaultNoise"+tag,seed,gridSize,gridSpacing,nonOrthogonalBasisReader(correlationFile_in))
  /* init */,correlationFile(correlationFile_in)
{
  // read the dimension of the original correlation
  const auto originalDimensions(readVTKfileDimension(correlationFile.c_str()));

  if(originalDimensions(2)!=1)
  {
      throw std::runtime_error("vtk stacking fault correlationFile 'DIMENSIONS' should have 3rd component == 1.");
  }

  const int originalNX = originalDimensions(0);
  const int originalNY = originalDimensions(1);
  const int originalNZ = originalDimensions(2);
  const int originalNR = originalNX*originalNY*originalNZ;

  // compute variable energy correction factor from zero padding
  eNormFactor = (this->NR!=originalNR ? static_cast<double>(this->NR)/static_cast<double>(originalNR) : 1.0);
  // fftScaleFactor is the scaling factor for FFT
  fftScaleFactor = std::pow(static_cast<double>(this->NR), 2.0);

  // original input correlation array
  REAL_SCALAR *Rr_original = (REAL_SCALAR*) fftw_malloc(sizeof(REAL_SCALAR)*originalNR);
  // explicitly initialize to zeros, avoid garbarge propagation
  std::fill(Rr_original, Rr_original+originalNR, REAL_SCALAR{0.0});

  // populate Rr_original with the correlation data
  StackingFaultCorrelationReader(correlationFile, Rr_original, originalNR);

  // unit conversion (from J^2/m^4 to unitless)
  for (int i = 0; i < originalNR; ++i)
  {
    const REAL_SCALAR unitconvert = mat.b_SI*mat.mu_SI;
    Rr_original[i] /= (unitconvert*unitconvert);
  }

  // allocate correlation array with zero padding in real space
  REAL_SCALAR *Rr = (REAL_SCALAR*) fftw_malloc(sizeof(REAL_SCALAR)*this->NR);
  std::fill(Rr, Rr+this->NR, REAL_SCALAR{0.0});

  // if the assigned grid size is the same as the actual data
  if(this->NR == originalNR)
  {
    // point to the same memory location
    Rr = Rr_original;
  }
  // pad with zeros if the assigned grid size is larger than the actual data grid size
  else
  {
    int start_y = (this->NY - originalNY) / 2;
    int start_x = (this->NX - originalNX) / 2;
    // 0-pading from centere
    for (int i = 0; i < originalNY; ++i)
    {
      for (int j = 0; j < originalNX; ++j)
      {
        Rr[(start_y + i) * this->NX + (start_x + j)] = Rr_original[i * originalNX + j];
      }
    }
    // Free allocated memory
    fftw_free(Rr_original);
  }

  // allocate correlation array in fourier space
  Rk = (COMPLEX*) fftw_malloc(sizeof(COMPLEX)*this->NK);
  std::fill(Rk, Rk+this->NK, COMPLEX{0.0, 0.0});

  // Execute FFTW plans to populate Rk
  fftw_plan plan_R_r2c = fftw_plan_dft_r2c_2d(this->NY, this->NX, Rr, reinterpret_cast<fftw_complex*>(Rk), FFTW_ESTIMATE);
  fftw_execute(plan_R_r2c);

  // Destroy FFTW plans
  fftw_destroy_plan(plan_R_r2c);

  // Free allocated memory
  fftw_free(Rr);
}

std::array<MDStackingFaultNoise::COMPLEX,1> MDStackingFaultNoise::kCorrelations(const Eigen::Matrix<double, 3, 1> &kv, const Eigen::Matrix<int, 3, 1> &index) const
{
    int idx=(this->NY/2+1)*NZ*index(0) + index(1)*NZ + index(2);
    // eNormFactor is variable energy correction factor from zero padding (1 if not padded)
    // fftScaleFactor is the Scaling factor for FFT (match the original input scale)
    std::array<MDStackingFaultNoise::COMPLEX,1> temp{eNormFactor*Rk[idx]/fftScaleFactor};
    return temp;
}

Eigen::Matrix<double,2,2> MDStackingFaultNoise::nonOrthogonalBasisReader(const std::string& fileName_vtk) const
{
    typedef Eigen::Matrix<double,2,1> VectorDimD;
    std::deque<VectorDimD> basis1;
    std::deque<VectorDimD> basis2;
    std::ifstream vtkFile(fileName_vtk); //access vtk file
    // error check
    if (!vtkFile.is_open()) {
        throw std::runtime_error("Error opening stacking fault VTK correlation file!");
    }
    // begin parsing structured_grid vtk file for lines
    std::string line;
    while (std::getline(vtkFile, line)) 
    {
        if(line.find("VECTORS")!=std::string::npos) 
        {
            std::getline(vtkFile, line);
            std::stringstream ss(line);
            double x, y, z;
            if(basis1.empty())
            {
                ss >> x >> y >> z;
                VectorDimD vec = (VectorDimD() << x, y).finished();
                vec.normalize(); // Normalize the vector
                basis1.push_back(vec);
            }
            else
            {
                ss >> x >> y >> z;
                VectorDimD vec = (VectorDimD() << x, y).finished();
                vec.normalize(); // Normalize the vector
                basis2.push_back(vec);
            }

        }
    }

    Eigen::Matrix<double,2,2> nonOrthoBasisMatrix;
    // fill the nonOrthoBasisMatrix with a block that has the size of <2,1> at the beginning of the first column
    nonOrthoBasisMatrix.block<2,1>(0,0) = basis1.front();
    // fill the nonOrthoBasisMatrix with a block that has the size of <2,1> at the beginning of the second column
    nonOrthoBasisMatrix.block<2,1>(0,1) = basis2.front();

    // transpose the matrix to create a transition matrix
    return nonOrthoBasisMatrix;
}

void MDStackingFaultNoise::StackingFaultCorrelationReader(const std::string& correlationFile, REAL_SCALAR* Rr, const int& NR)
{
  std::cout << "Reading stacking fault correlation" << std::endl;

  std::ifstream vtkFile(correlationFile); // access vtk file
  // error check
  if (!vtkFile.is_open())
  {
    throw std::runtime_error("Error opening stacking fault VTK correlation file!");
  }

  std::string line;
  while (std::getline(vtkFile, line))
  {
    // if the "POINT_DATA" string is read, read the following data
    if (line.find("POINT_DATA") != std::string::npos)
    {
      const size_t numOfHeaders = 2;
      // get the number of points in the file
      const size_t firstSpace(line.find(' '));
      const size_t numOfPoints = std::atoi(line.substr(firstSpace+1).c_str());
      // read the point coordinates
      for (size_t n = 0; n < numOfPoints + numOfHeaders; ++n)
      {
        std::getline(vtkFile, line);
        // ignore the headers that come right after point_data
        if (n < numOfHeaders)
          continue;
        const int ind = n - numOfHeaders;
        // correlationCoeffs.push_back(std::atoi(line.c_str()));
        if (ind >= NR) 
        {
          throw std::runtime_error("Index out of bounds while populating the original correlation array.");
        }

        try
        {
          double value = std::stod(line);
          //float value = std::atof(line.c_str());
          Rr[ind] = value;
        }
        catch(const std::invalid_argument& e) 
        {
          std::cerr << "Invalid correlation data in line: " << line << std::endl;
        }
        catch(const std::out_of_range& e)
        {
          std::cerr << "Out of range correlation value in line: " << line << std::endl;
        }
      }
    }
  }
  vtkFile.close(); // Close the file after reading
}

typename MDStackingFaultNoise::GridSizeType MDStackingFaultNoise::readVTKfileDimension(const char *fname)
{
    FILE *InFile=fopen(fname,"r");

    if (InFile == NULL)
    {
        fprintf(stderr, "Can't open stacking fault correlation VTK file %s\n",fname);
        exit(1);
    }
    // return the 5th line of the vtk file
    char line[200];
    for(int i=0;i<4;i++)
    {
        fgets(line, 200, InFile);
    }
    // scan the returned line
    int NXX, NYY, NZZ;
    fscanf(InFile, "%s %d %d %d\n", line, &(NXX), &(NYY), &(NZZ));
    return (GridSizeType()<<NXX,NYY,NZZ).finished();
}

}
#endif
