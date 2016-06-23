README FILE: Writing pkls from raw .ams or .amz spectra files

To write .pkl files (compatible with FATES) run redopkls if the raw data is in the .ams file type or redopkls_AMZ if the raw individual spectra are saved as .amz files.  These (.ams and .amz) file types are commonly output from ATOFMS and TSI-ATOFMS instruments.  The workflow called by the scripts is as follows

user calls redopkls
redopkls iterates through  top directory looking for .ams files.
When .ams files are found AMStoPKL_freshStart is called.
AMStoPKL_freshStart aggregates data within all .ams files (time space) within a folder into a single .pkl (m/z space) file.
	AMStoPKL_freshStart calls get_spectrumAMS to import raw .ams data into MATLAB
	AMStoPKL_freshStart then calls PeakList_gen_* (either 	LVNFinal or SLYFinal) to convert
	 the raw spectra (often consisting of 15,000 points) to a "peak picked" m/z calibrated spectra only 
	containg data points relevant to peaks, thus minimizing storage demands.  PeakList_gen_* performs 
	baselining, noise removal, and peak picking using thresholds and algorithims empirically 
	determined to be effective for raw spectra generated with ATOFMS within Kim Prather's lab at UCSD.
	PEAKList_gen_LVNFinal is used for linear reflectron TOFs 
	PEAKList_gen_SLYFinal is used for TofWerks Z-tof
	
IMPORTANT NOTES: 
The user must specify within AMStoPKL_freshStart whether to call PEAKList_gen_LVNFinal or PEAKList_gen_SLYFinal! 
All baselining and peak picking within these files has been empirically optimized for use with ATOFMS within Kim
Prather's lab.  The thresholds and algorithms used are not necessarilly appropriate for all SPMS raw spectra!
	
