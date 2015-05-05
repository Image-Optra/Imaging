/**
 *  @file  PatchExtractor.cpp
 *
 *  @brief  Implementation of the PatchExtractor.
 *
 *  Implementation of the PatchExtractor.  The PatchExtractor is a command-line program
 *  for extracting and correcting patches from a runfile into a desired directory. It
 *  operates in two modes. In the first mode, the expert classification file is not used
 *  and all the patches in each runfile are extracted to a directory created specifically
 *  for that runfile.  The second mode uses the expert classification file.  All the
 *  patches that belong to a particular class are extracted to a directory for that class
 *  regardless of the runfile from which the patch originated.  In the first mode, there
 *  is a directory for each runfile, so the labeling of patches only reflects their
 *  positions within their source runfile.  In the second mode, all the patches of a
 *  given class are stored together in one directory, so the labeling of the patches
 *  reflects the runfiles from which the patches come.
 *
 *  Copyright &copy; 2014  -  IRIS International, Inc.  -  All rights reserved
 */

  #include <ISL/APR/Runfile.h>
  #include <ISL/APR/Particle.h>
  #include <ISL/APR/Features.h>
  #include <ISL/APR/Calculators.h>

  #include <ISL/Image/BayerImage.h>
  #include <ISL/Image/Debayering.h>
  #include <ISL/Image/DirectImage.h>
  #include <ISL/Image/GrayscaleImage.h>
  #include <ISL/Image/RGB_Image.h>
  #include <ISL/Math/Matrix.h>

  #include <boost/filesystem.hpp>

  #include <iostream>
  #include <iomanip>

  #include "ClassificationList.h"


//-----------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------

    namespace APRT
      {

/**
 *  @brief  A class for extracting patches from runfiles and storing the patches in
 *          directories corresponding to either their runfiles or their patch
 *          classifications.
 */

        class PatchExtractor
          {
            public:
              PatchExtractor(const std::string destination,
                             const uint8_t     sample);
                /**< @brief  creates a PatchExtractor for a
                             runfilelist and subsample number */

            public:
              void  Sort(const std::string runfilelist);
                /**< @brief  a driver function used to iterate through a
                             runfile list to extract all the patches of a specific
                             classification to a single directory for that type of
                             patch, ideal for optimizing classifiers and feature
                             generators over particular classes/types of patches */

            private:
              void  WriteSort(const std::string runfilename);
                /**< @brief  a worker function that writes the contents of a
                             runfile to directories created for their patch types */

            private:
              std::string  outputdirectory;
                /**< @brief  the output directory containing images */
              std::string  inputdirectory;
                /**< @brief  the input directory containing runfiles */
              const uint8_t subsamplenumber;
                /**< @brief  the runfile subsample (stream) to write */
          };

/**
 *  @brief  An external function to create and run a PatchExtractor to write particles
 *          contained in all the runfiles listed on a runfilelist into directories
 *          associated created for their particle types.
 */

        void Sort(const std::string runfilelist,
                  const std::string destination,
                  const uint8_t     sample);
      }


//-----------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------

/**
 *  Creates an APRT::PatchExtractor object for extracting patches contained in runfiles
 *  listed on a runfilelist to either their respective expert classification folders or
 *  to a single folder for each runfile on the list.
 *
 *  @param [in]  destination  the output destination
 *  @param [in]  runfilelist  the subsample number
 */

  APRT::PatchExtractor::PatchExtractor(const std::string destination,
                                       const uint8_t     sample)
   : outputdirectory(destination),
     subsamplenumber(sample)
      {
        ;
      }


//-----------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------

/**
 *  A driver function used to iterate through a runfile list to extract all the patches
 *  of a specific classification to a single directory for that type of patch. This form
 *  of output is ideal for optimizing classifiers and feature generators over particular
 *  classes/types of patches.
 *
 *  @param [in]  runfilelist  the input list of runfiles
 */

  void APRT::PatchExtractor::Sort(const std::string runfilelist)
    {
//
//  Read the input list of runfiles ...
//
      std::ifstream runfileliststream(runfilelist.c_str());
//
//  Get the output runfile directory ...
//
      std::getline(runfileliststream,this->inputdirectory);
//
//  Process each listed runfile in turn ...
//
      while (!runfileliststream.eof())
        {
//
//  Get the runfilename to process ...
//
          std::string nextline;
          std::getline(runfileliststream,nextline);
          std::cout << "Processing -> "
                    << nextline.c_str()
                    << std::endl;
          this->WriteSort(nextline.c_str());
        }
    }


//-----------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------

/**
 *  A worker function that writes the contents of a runfile to directories common to their
 *  patch types. This is ideal for optimizing the features and classifiers on all the
 *  particles of a particular class contained in a group of runfiles.
 *
 *  @param [in]  runfilename  the input runfile name
 */

  void APRT::PatchExtractor::WriteSort(const std::string runfilename)
    {
//
//  Read the classification file ...
//
      std::ifstream
          pclfilestream((this->inputdirectory +
                         runfilename +
                         ".pcl").c_str(),
                        std::ios_base::in);
      APRT::ClassificationList
          pclpatchlist = ClassificationList(pclfilestream);

	  std::ifstream
          aclfilestream((this->inputdirectory +
                         runfilename +
                         ".acl").c_str(),
                        std::ios_base::in);
      APRT::ClassificationList
          aclpatchlist = ClassificationList(aclfilestream);
//
//  Schedule the particles in the runfile subsample in turn ...
//
	  ISL::Math::Matrix<int32_t,2> conmatrix(26,26);
      uint32_t count = 0;
      while ((count < pclpatchlist.Classifications()[this->subsamplenumber-1].size()) &&
		     (count < aclpatchlist.Classifications()[this->subsamplenumber-1].size()))
        {
          std::string pclclassification =
              std::string(pclpatchlist.Classifications()[this->subsamplenumber-1][count].classification);
		  std::string aclclassification =
              std::string(aclpatchlist.Classifications()[this->subsamplenumber-1][count].classification);
		  uint32_t pclindex = 25,aclindex = 25;
		  if (pclclassification.compare("RBC") == 0) pclindex = 0;
		  else if (pclclassification.compare("DRBC") == 0) pclindex = 1;
		  else if (pclclassification.compare("RBCC") == 0) pclindex = 2;
		  else if (pclclassification.compare("WBC") == 0)  pclindex = 3;
		  else if (pclclassification.compare("WBCC") == 0) pclindex = 4;
		  else if (pclclassification.compare("BACT") == 0) pclindex = 5;
		  else if (pclclassification.compare("SQEP") == 0) pclindex = 6;
		  else if (pclclassification.compare("NSE") == 0)  pclindex = 7;
		  else if (pclclassification.compare("TREP") == 0) pclindex = 8;
		  else if (pclclassification.compare("REEP") == 0) pclindex = 9;
		  else if (pclclassification.compare("CAOX") == 0) pclindex = 10;
		  else if (pclclassification.compare("URIC") == 0) pclindex = 11;
		  else if (pclclassification.compare("TPO4") == 0) pclindex = 12;
		  else if (pclclassification.compare("CAPH") == 0) pclindex = 13;
		  else if (pclclassification.compare("CYST") == 0) pclindex = 14;
		  else if (pclclassification.compare("LEUC") == 0) pclindex = 15;
		  else if (pclclassification.compare("AMOR") == 0) pclindex = 16;
		  else if (pclclassification.compare("CELL") == 0) pclindex = 17;
		  else if (pclclassification.compare("GRAN") == 0) pclindex = 18;
		  else if (pclclassification.compare("MUCS") == 0) pclindex = 19;
		  else if (pclclassification.compare("SPRM") == 0) pclindex = 20;
		  else if (pclclassification.compare("BYST") == 0) pclindex = 21;
		  else if (pclclassification.compare("HYST") == 0) pclindex = 22;
		  else if (pclclassification.compare("TRCH") == 0) pclindex = 23;
		  else if (pclclassification.compare("BUBB") == 0) pclindex = 24;
		  else if (pclclassification.compare("NONE") == 0) pclindex = 25;

		  if (aclclassification.compare("RBC") == 0) aclindex = 0;
		  else if (aclclassification.compare("DRBC") == 0) aclindex = 1;
		  else if (aclclassification.compare("RBCC") == 0) aclindex = 2;
		  else if (aclclassification.compare("WBC") == 0)  aclindex = 3;
		  else if (aclclassification.compare("WBCC") == 0) aclindex = 4;
		  else if (aclclassification.compare("BACT") == 0) aclindex = 5;
		  else if (aclclassification.compare("SQEP") == 0) aclindex = 6;
		  else if (aclclassification.compare("NSE") == 0)  aclindex = 7;
		  else if (aclclassification.compare("TREP") == 0) aclindex = 8;
		  else if (aclclassification.compare("REEP") == 0) aclindex = 9;
		  else if (aclclassification.compare("CAOX") == 0) aclindex = 10;
		  else if (aclclassification.compare("URIC") == 0) aclindex = 11;
		  else if (aclclassification.compare("TPO4") == 0) aclindex = 12;
		  else if (aclclassification.compare("CAPH") == 0) aclindex = 13;
		  else if (aclclassification.compare("CYST") == 0) aclindex = 14;
		  else if (aclclassification.compare("LEUC") == 0) aclindex = 15;
		  else if (aclclassification.compare("AMOR") == 0) aclindex = 16;
		  else if (aclclassification.compare("CELL") == 0) aclindex = 17;
		  else if (aclclassification.compare("GRAN") == 0) aclindex = 18;
		  else if (aclclassification.compare("MUCS") == 0) aclindex = 19;
		  else if (aclclassification.compare("SPRM") == 0) aclindex = 20;
		  else if (aclclassification.compare("BYST") == 0) aclindex = 21;
		  else if (aclclassification.compare("HYST") == 0) aclindex = 22;
		  else if (aclclassification.compare("TRCH") == 0) aclindex = 23;
		  else if (aclclassification.compare("BUBB") == 0) aclindex = 24;
		  else if (aclclassification.compare("NONE") == 0) aclindex = 25;
		  ++conmatrix(pclindex,aclindex); 
          ++count;
        }
	   std::string basefolder =
          std::string(this->outputdirectory + "/");
      boost::filesystem::path base(basefolder);
	  std::ostringstream confile;
          confile
              << basefolder
              << "ConfusionMatrix.txt";
             
	  FILE *fp = fopen(confile.str().c_str(),"a+");
	  for (uint32_t i = 0; i < conmatrix.Dim1(); ++i)
	  {
		  for (uint32_t j = 0; j < conmatrix.Dim2(); ++j)
		  {
			  fprintf(fp,"%d\t",conmatrix(i,j)); 
		  }
		  fprintf(fp,"\n");
	  }
	  
	  fclose(fp);
    }


//-----------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------

/**
 *  An external function to create and run a PatchExtractor to write particles
 *  contained in all the runfiles listed on a runfilelist into directories associated
 *  created for their particle types. This is ideal for opitmizing a feature selector
 *  or classifier over images of particular types/classes obtained from a collection
 *  of runfiles in a runfilelist.
 *
 *  @param [in]  runfilelist  the list of runfiles to extract
 *  @param [in]  destination  the output image directory
 *  @param [in]  sample       the runfile sample number of interest
 */

  void APRT::Sort(const std::string runfilelist,
                  const std::string destination,
                  const uint8_t     sample)
    {
//
//  Extract the patches contained in the runfile listed in the runfilelist
//  into the output image directories ...
//
      PatchExtractor extractor = PatchExtractor(destination,sample);
      extractor.Sort(runfilelist);
//
//  Characterize the contents of the output directories ...
//
      /* Skipped for now.  This would produce the "Runfile List
         Statistics Report" file. */
    }


//-----------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------

/**
 *  The main entry point to the program.
 *
 *  @param [in]  argc  the number of input arguments
 *  @param [in]  argv  the strings of input arguments
 *
 *  @return  EXIT_FAILURE (upon an exception or forced termination)
 */

  int main(int argc, char* argv[])
    {
      try
        {
          //if (argc == 5)
            {
              const std::string runfilelist = "E:/rundata/runfile.txt";//argv[1];
              const std::string destination = "D:/Dataout/";//argv[2];
              const int         subsample   = 1;//std::atoi(argv[3]);

              std::cout << "Readying "
                        << runfilelist
                        << " for processing."
                        << std::endl;
              APRT::Sort(runfilelist,destination,subsample);
            }
          //else
            {
              std::cout << "Invalid argument list. Try again." << std::endl;
            }
        }

      catch (const std::runtime_error& e)
        {
          std::cout << e.what() << std::endl;
        }

      catch (...)
        {
          std::cout << "Oops! Not good. Not good at all!" << std::endl;
          return (EXIT_FAILURE);
        }

      return (EXIT_FAILURE);
    }




              



