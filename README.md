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




Signals to be downloaded is in the following link
-------------------------------------------------------------------------------
----------- Beidou_B3I_IF_signal.bin --------------------------
link: https://pan.baidu.com/s/14X937mHKtB4pxZKPPlDx6A password: p9v9

-- dataType: schar
-- fileType: 8 bit real samples
-- IF: 8.52e6 Hz
-- sampling Frequency: 40e6 Hz
----------- B1C_fs_99.375_if14.58.bin --------------------------
(link:): https://pan.baidu.com/s/12X4Ov-eSDb7w11CydI6ryg password: srgh

-- dataType: schar
-- fileType: 8 bit real samples
-- IF: 14.58e6 Hz
-- sampling Frequency: 99.375e6 Hz
----------- B1I_IF3.902_Fs23.8.bin --------------------------
(link:): https://pan.baidu.com/s/1sjvP9ONP7fh9vR6a3_bGzg password: z95r

-- dataType: schar
-- fileType: 8 bit real samples
-- IF: 14.58e6 Hz
-- sampling Frequency: 99.375e6 Hz
----------- GPS_L2C_IF_signal.bin --------------------------
(link:): https://pan.baidu.com/s/1TB4NBtFG0Bgx3iOboyHObg password: 7kqe

-- dataType: schar
-- fileType: 8 bit real samples
-- IF: 7.4e6 Hz
-- sampling Frequency: 53e6 Hz
----------- GPS_L1_CA_IF_signal.bin --------------------------
(link:) https://pan.baidu.com/s/1Y78Qlxo2rgIcGXDHiPQ_aw password: 2dr5

-- dataType: schar
-- fileType: 8 bit real samples
-- IF: 14.58e6 Hz
-- sampling Frequency: 53e6 Hz
----------- IF_Data_Set_L5C.bin --------------------------
(link:): https://pan.baidu.com/s/1L7yyjzLAHsUKkRdLqEIWOQ password: 7r94

-- dataType: schar
-- fileType: 8 bit real samples
-- IF: 13.55e6 Hz
-- sampling Frequency: 99.375e6 Hz
----------- Galileo_E5a_IF_signal.bin --------------------------
link: https://pan.baidu.com/s/1e0RJsq1c8sfoxQvs6ItO2Q password: 1j33

-- dataType: schar
-- fileType: 8 bit real samples
-- IF: 13.55e6 Hz
-- sampling Frequency: 99.375e6 Hz
----------- BDS_B2a_IF_signal.bin --------------------------
(link:): https://pan.baidu.com/s/1jqrCt9BcakgqDFUFiKl50Q password: kf4c

-- also for Galileo E5a
-- dataType: schar
-- fileType: 8 bit real samples
-- IF: 13.55e6 Hz
-- sampling Frequency: 99.375e6 Hz
----------- BDS_B2a_IF_signal.bin --------------------------
(link:): https://pan.baidu.com/s/1z-7jWyRvdNVqHM51-Ac16Q password: 5a7d

-- also for Galileo E5b
-- dataType: schar
-- fileType: 8 bit real samples
-- IF: 17.14e6 Hz
-- sampling Frequency: 99.375e6 Hz
