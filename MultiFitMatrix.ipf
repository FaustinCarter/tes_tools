#pragma rtGlobals=3		// Use modern global access method and strict wave access.
// Function designed to take as input the name of a Multi-peak Fit 2
// set folder containing three fit waves: Baseline, Peak 0 and Peak 1, along with a
// data matrix. This function loops through the matrix using
// the coefficient waves found in the data folder as the initial guesses
// for each of the data waves found in the matrix.


Function MultiPeakBatchFit(fitFolder, dataMat)
	String fitFolder
	Wave dataMat
	
	Variable nfits = DimSize(dataMat,1)
	
	
	// The structure required by MPF2_DoMPFit()
	STRUCT MPFitInfoStruct MPStruct

	Variable i, err

	String saveDF = GetDataFolder(1)

		//Types of coefficient waves; modify as needed
		MPStruct.ListOfFunctions = "Constant;ExpConvExp;ExpConvExp;"

		// Create the list of coefficient waves and hold strings
		MPStruct.ListOfCWaveNames = ""
		
		//Fix everything but the peak heights
		//The ordering is different for different coef wave types
		//Be careful!
		MPStruct.ListOfHoldStrings = "0;1011;1011"

		//Names of the coefficient waves
		MPStruct.ListOfCWaveNames += "'Baseline Coefs';"+"'Peak 0 Coefs';"+"'Peak 1 Coefs';"
		

		MPStruct.NPeaks = 2
		MPStruct.XPointRangeBegin = 0
		// this function always fits the entire data set
		MPStruct.FitCurvePoints = 1000
		MPStruct.fitOptions = 4	// you might want to change this to 4
	
		// loop through the matrix, calling MPF2_DoMPFit
		// on each column. After each call, the results are checked and acted
		// on appropriately.
	
	//Create a temporary holder to copy each column into
	Make/O/N=(DimSize(dataMat,0)) temp
	CopyScales dataMat, temp
	
	SetDataFolder fitFolder
	
	//fast and slow are arbitrary names for peak 0 and peak 1
	Wave fastCoefs, slowCoefs
	
	if(waveexists(fastCoefs))
		KillWaves fastCoefs
	endif
	
	if(waveexists(slowCoefs))
		KillWaves slowCoefs
	endif
	
	for (i = 0; i < nfits; i += 1)

		// Most of the contents of the structure don't need to change.
		// Naturally, we have to put in a wave reference to the
		// data set being fit this time through the loop
		
		//pull out the appropriate column of the matrix
		temp=dataMat[p][i]
		Wave MPStruct.yWave = temp


		// just in case the new wave has a different number of points
		MPStruct.XPointRangeEnd = numpnts(MPStruct.yWave)-1
		
		err = MPF2_DoMPFit(MPStruct, fitFolder)
		if (err)
			// error return from MPF2_DoMPFit generally indicates
			// a programmer error
			DoAlert 0, "Error calling MPF2_DoMPFit: "+num2str(err)
			return err
		endif
		if (MPStruct.fitQuitReason == 2)
			// if the user aborts a fit, we assume that the whole process
			// should be aborted
			DoAlert 0, "User aborted batch fit"
			return -1
		endif
		if (MPStruct.fitError)
			// Something went wrong with the current fit, and it
			// failed to converge. We note this fact to the user via
			// an alert, then move on to the next one. You may wish
			// to do something more sophisticated than this.
			
			// Avoid a long line by concatenating the message
			// in small pieces
			String alertMsg = "Error doing fit to "
			alertMsg += NameofWave(dataMat)+"[]["+num2str(i)+"]: "
			alertMsg += num2str(MPStruct.fitError)
			alertMsg += "; continuing with next fit."
			DoAlert 0, alertMsg
		endif
		
		//make a matrix containing all the fit coefficients for each peak type
		String s=fitFolder+"'Peak 0 Coefs'"
		Concatenate {$s}, fastCoefs
		
		s=fitFolder+"'Peak 1 Coefs'"
		Concatenate {$s}, slowCoefs
		
		//Super low-tech progress bar
		if(mod(i, 1000)==0)
			print num2str(i)
		endif
		
	endfor
	
	//The fitting routine creates a lot of crap I don't care about. Reference all of it, then kill it!
	wave blepswave, 'Peak 0 Coefseps', 'Peak 1 Coefseps', fit_temp, Res_temp, M_Covar, W_sigma, W_sigma_0, W_sigma_1, W_sigma_2
	KillWaves temp, blepswave, 'Peak 0 Coefseps', 'Peak 1 Coefseps', fit_temp, Res_temp, M_Covar, W_sigma, W_sigma_0, W_sigma_1, W_sigma_2
	
	//step back out to whatever datafolder was current when this whole thing started.
	SetDataFolder saveDF
end