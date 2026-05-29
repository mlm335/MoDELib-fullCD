/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_MDSolidSolutionNoise_cpp
#define model_MDSolidSolutionNoise_cpp

#include <MDSolidSolutionNoise.h>

namespace model
{
MDSolidSolutionNoise::MDSolidSolutionNoise(const PolycrystallineMaterialBase& mat,
                                           const std::string& tag,
                                           const std::string& correlationFile_xz,
                                           const std::string& correlationFile_yz,
                                           const int& seed,
                                           const GridSizeType& gridSize,
                                           const GridSpacingType& gridSpacing,
                                           const Eigen::Matrix<double,2,2>& latticeBasis,
                                           const double& a_Cai_in
                                           ) :
  /*init*/ GlidePlaneNoiseBase<2>("MDSolidSolutionNoise"+tag,seed,gridSize,gridSpacing,latticeBasis)
  /*init*/,a_cai(a_Cai_in)
{
  // read the dimensions of the original correlation sampled from MD
  const auto originalDimensions_xz(readVTKfileDimension(correlationFile_xz.c_str()));
  const auto originalDimensions_yz(readVTKfileDimension(correlationFile_yz.c_str()));

  if((originalDimensions_xz-originalDimensions_yz).matrix().squaredNorm()>0)
  {
    throw std::runtime_error("correlationFile_xz and correlationFile_yz have different grid sizes.");
  }

  if(originalDimensions_xz(2)!=1 || originalDimensions_yz(2)!=1)
  {
    throw std::runtime_error("vtk solid solution correlationFiles 'DIMENSIONS' should have 3rd component == 1.");
  }

  const int originalNX = originalDimensions_xz(0);
  const int originalNY = originalDimensions_xz(1);
  const int originalNZ = originalDimensions_xz(2);
  const int originalNR = originalNX*originalNY*originalNZ;

  // compute variable energy correction factor from zero padding
  eNormFactor = (this->NR!=originalNR ? static_cast<double>(this->NR)/static_cast<double>(originalNR) : 1.0);
  // fftScaleFactor is the scaling factor for FFT
  fftScaleFactor = std::pow(static_cast<double>(this->NR),2.0);

  REAL_SCALAR *Rr_xz_original = (REAL_SCALAR*) fftw_malloc(sizeof(REAL_SCALAR)*originalNR);
  REAL_SCALAR *Rr_yz_original = (REAL_SCALAR*) fftw_malloc(sizeof(REAL_SCALAR)*originalNR);
  // explicitly initialize to zeros, avoid garbarge propagation
  std::fill(Rr_xz_original, Rr_xz_original+originalNR, REAL_SCALAR{0.0});
  std::fill(Rr_yz_original, Rr_yz_original+originalNR, REAL_SCALAR{0.0});

  // populate Rr_original with the correlation data from MD
  SolidSolutionCorrelationReader(correlationFile_xz, Rr_xz_original, originalNR);
  SolidSolutionCorrelationReader(correlationFile_yz, Rr_yz_original, originalNR);

  // Divide the values in Rr_xz_original and Rr_yz_original by mat.mu^2
  for (int i = 0; i < originalNR; ++i)
  {
    Rr_xz_original[i] /= (mat.mu_SI*mat.mu_SI); // divide by mu^2
    Rr_yz_original[i] /= (mat.mu_SI*mat.mu_SI); // divide by mu^2
  }

  // allocate real space correlation
  REAL_SCALAR *Rr_xz = (REAL_SCALAR*) fftw_malloc(sizeof(REAL_SCALAR)*this->NR);
  std::fill(Rr_xz, Rr_xz+this->NR, REAL_SCALAR{0.0});
  REAL_SCALAR *Rr_yz = (REAL_SCALAR*) fftw_malloc(sizeof(REAL_SCALAR)*this->NR);
  std::fill(Rr_yz, Rr_yz+this->NR, REAL_SCALAR{0.0});

  // if the assigned grid size is the same as the actual data
  if(this->NR == originalNR)
  {
    // point to the same memory location
    Rr_xz = Rr_xz_original;
    Rr_yz = Rr_yz_original;
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
        Rr_xz[(start_y + i) * this->NX + (start_x + j)] = Rr_xz_original[i * originalNX + j];
        Rr_yz[(start_y + i) * this->NX + (start_x + j)] = Rr_yz_original[i * originalNX + j];
      }
    }
    // Free allocated memory
    fftw_free(Rr_xz_original);
    fftw_free(Rr_yz_original);
  }

  // initialize the correlation in fourier space
  Rk_xz = (COMPLEX*) fftw_malloc(sizeof(COMPLEX)*this->NK);
  std::fill(Rk_xz, Rk_xz+this->NK, COMPLEX{0.0, 0.0});
  Rk_yz = (COMPLEX*) fftw_malloc(sizeof(COMPLEX)*this->NK);
  std::fill(Rk_yz, Rk_yz+this->NK, COMPLEX{0.0, 0.0});

  // Initialize FFTW plans as member variables
  fftw_plan plan_R_xz_r2c = fftw_plan_dft_r2c_2d(this->NY, this->NX, Rr_xz, reinterpret_cast<fftw_complex*>(Rk_xz), FFTW_ESTIMATE);
  fftw_plan plan_R_yz_r2c = fftw_plan_dft_r2c_2d(this->NY, this->NX, Rr_yz, reinterpret_cast<fftw_complex*>(Rk_yz), FFTW_ESTIMATE);

  // Execute FFTW plans to populate Rk_xz and Rk_yz
  fftw_execute(plan_R_xz_r2c);
  fftw_execute(plan_R_yz_r2c);

  // Destroy FFTW plans
  fftw_destroy_plan(plan_R_xz_r2c);
  fftw_destroy_plan(plan_R_yz_r2c);

  // Free allocated memory
  fftw_free(Rr_xz);
  fftw_free(Rr_yz);
}

std::array<MDSolidSolutionNoise::COMPLEX, 2> MDSolidSolutionNoise::kCorrelations(const Eigen::Matrix<double, 3, 1> &kv, const Eigen::Matrix<int, 3, 1> &index) const
{
  int idx=(this->NY/2+1)*NZ*index(0) + index(1)*NZ + index(2);

  // eNormFactor is variable energy correction factor from zero padding (1 if not padded)
  // fftScaleFactor is the Scaling factor for FFT (match the original input scale)
  std::array<MDSolidSolutionNoise::COMPLEX, 2> temp{
    eNormFactor*Rk_xz[idx]/fftScaleFactor,
    eNormFactor*Rk_yz[idx]/fftScaleFactor
  };

  if(a_cai>0.0)
  {
    const double wkc2(this->Wk_Cai_squared(kv(0),kv(1),kv(2), a_cai)); // using the square because this is before the square root
    temp[0]*=wkc2;
    temp[1]*=wkc2;
  }
  return temp;
}

void MDSolidSolutionNoise::SolidSolutionCorrelationReader(const std::string& correlationFile, REAL_SCALAR* Rr, const int& NR)
{
  std::cout << "Reading solid solution correlation" << std::endl;

  std::ifstream vtkFile(correlationFile); //access vtk file
  // error check
  if (!vtkFile.is_open())
  {
    throw std::runtime_error("Error opening solid solution VTK correlation file!");
  }

  std::string line;
  while (std::getline(vtkFile, line))
  {
    //if the "POINT_DATA" string is read, read the following data
    if(line.find("POINT_DATA")!=std::string::npos)
    {
      const size_t numOfHeaders = 2;
      // get the number of points in the file
      const size_t firstSpace(line.find(' '));
      const size_t numOfPoints = std::atoi(line.substr(firstSpace+1).c_str());
      // read the point coordinates
      for(size_t n=0; n<numOfPoints+numOfHeaders; ++n)
      {
        std::getline(vtkFile, line);
        // ignore the headers that come right after point_data
        if(n<numOfHeaders)
          continue;
        const int ind = n-numOfHeaders;

        if (ind >= NR)
        {
          throw std::runtime_error("Index out of bounds while populating the original correlation array.");
        }

        try
        {
          double value = std::stod(line);
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
  vtkFile.close();
}

typename MDSolidSolutionNoise::GridSizeType MDSolidSolutionNoise::readVTKfileDimension(const char *fname)
{
  FILE *InFile=fopen(fname,"r");

  if (InFile == NULL)
  {
    fprintf(stderr, "Can't open solid solution correlation VTK file %s\n",fname);
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
