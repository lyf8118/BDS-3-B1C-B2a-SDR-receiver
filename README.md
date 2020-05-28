An Open-Source BDS3 B1C/B2a SDR 
===============================================================================



Authors
-------------------------------------------------------------------------------
Yafeng Li
E-Mail: <lyf8118@126.com> <yali2822@colorado.edu>

Nagaraj Channarayapatna Shivaramaiah
E-Mail: <nagarajcs@colorado.edu>

Dennis Akos  
E-Mail: <dma@colorado.edu>
HP: <http://www.colorado.edu/aerospace/dennis-akos>




Features
-------------------------------------------------------------------------------
* GNSS signal processing functions written in MATLAB
    * Code generations
    * Signal acquisition / tracking (data + pilot)
    * Decoding navigation messages 
    * Pseudo-range mesurement generation
    * Calculation of position
* Support following signals
	* Beidou Phase III B2a
	* Beidou Phase III B1C
* Support RF binary file for post processing
    * All the SDRs have been tested with IF signals collected by the NUT4NT 
	sampler of the Amungo Navigation company 

	

Directory and Files
-------------------------------------------------------------------------------
./Doc                     Summary PowerPoint documents for each SDR receiver
./Common                  Common functions between differnt SDR receivers
./IF_Data_Set             Folder containing the IF data sets to be processed and 
                          corresponding Metadata files
./BDS_B1C                 Beidou B1/B2 SDR receiver
    ./include             Functions related to NAV data decoding and PVT computation  
    ./init.m              Strating function of the receiver  
    ./initSettings.m      Parameter configurations of the receiver    
    ./postProcessing.m    Top level processing function of the receiver  
    ./acquisition.m       Signal acquisition function  
    ./tracking.m          Signal tracking function 
    ./postNavigation.m    PVT computation
./BDS_B2a                 Beidou B1/B2 SDR receiver
    ./include             Functions related to NAV data decoding and PVT computation  
    ./init.m              Strating function of the receiver  
    ./initSettings.m      Parameter configurations of the receiver    
    ./postProcessing.m    Top level processing function of the receiver  
    ./acquisition.m       Signal acquisition function  
    ./tracking.m          Signal tracking function 
    ./postNavigation.m    PVT computation
 


Software Dependencies
-------------------------------------------------------------------------------
* All the SDRs have been tested on MATLAB R2016b and R2013b. Some functions 
  of the MATLAB Communications System Toolbox are needed:
  -- Create Galois field array: gf()
  -- BCH decoder: bchdec()
  -- Detect errors in input data using CRC: comm.CRCDetector()
  -- Convert convolutional code polynomials to trellis description: poly2trellis()
  -- Convolutionally decode binary data using Viterbi algorithm: vitdec()
  -- Detect errors in input data using CRC: step() 
  
  
 
How to use
-------------------------------------------------------------------------------
* Step 1: Copy the IF data file into folder "IF_Data_Set";
* Step 2: Configue parameters related to IF data file in function "initSettings.m";
* Step 2: Start processing by runing function "init.m".



Implementation details
-------------------------------------------------------------------------------
See the reference:
Li, Y., Shivaramaiah, N.C. & Akos, D.M. Design and implementation of an open-source BDS-3 B1C/B2a SDR receiver. GPS Solut (2019) 23: 60. https://doi.org/10.1007/s10291-019-0853-z
