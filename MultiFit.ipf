#pragma rtGlobals=3		// Use modern global access method and strict wave access.
// Function designed to take as input the name of a Multi-peak Fit 2
// set folder found in root:Packages:MultiPeakFit2:, along with a
// list of input data waves. This function loops through the list using
// the coefficient waves found in the data folder as the initial guesses
// for each of the data waves found in ywavelist, xwavelist.

//To create ywavelist, use the function below: MakeList

Function MultiPeakBatchFit(datafolderName,destinationFolder, ywavelist[, xwavelist])
	String datafolderName, destinationFolder, ywavelist, xwavelist
	
	Variable nfits = ItemsInList(ywavelist)
	
	Make/O/N=(nfits) pulseHeights
	//Make/O/N=(nfits) pulseAreas
	Make/O/T/N=(nfits) pulseNames
	
	// The structure required by MPF2_DoMPFit()
	STRUCT MPFitInfoStruct MPStruct

	Variable i, err

	// Switch to the Multi-peak fit folder and look up some stuff
	String saveDF = GetDataFolder(1)
	//SetDataFolder datafolderName
		MPStruct.ListOfFunctions = "Constant;ExpModGauss;ExpModGauss;"

		// Create the list of coefficient waves and hold strings
		MPStruct.ListOfCWaveNames = ""
		MPStruct.ListOfHoldStrings = "0;1101;1111"

		MPStruct.ListOfCWaveNames += "'Baseline Coefs';"+"'Peak 0 Coefs';"+"'Peak 1 Coefs';"
		
		
	//SetDataFolder saveDF
	
	// fill in more members of the structure
	MPStruct.NPeaks = 2
	MPStruct.XPointRangeBegin = 0
	// this function always fits the entire data set
	MPStruct.FitCurvePoints = 1000
	MPStruct.fitOptions = 4	// you might want to change this to 4
	
	// loop through each of the data sets, calling MPF2_DoMPFit
	// on each one. After each call, the results are checked and acted
	// on appropriately.
	for (i = 0; i < nfits; i += 1)
		// We make a copy of the Multi-peak fit folder for
		// each data set. That also duplicates all the coefficient waves
		// so that we preserve the results for each fit. The duplicated data folder
		// is created inside the current data folder.
		
		String BatchDFName = destinationFolder + "BatchPeakFit_"+num2str(i)
		DuplicateDataFolder $datafolderName, $BatchDFName
		
		setDataFolder BatchDFName
		
		// Most of the contents of the structure don't need to change.
		// Naturally, we have to put in a wave reference to the
		// data set being fit this time through the loop
		Wave MPStruct.yWave = $StringFromList(i, ywavelist)

		// the X wave for the current data set won't exist if the data
		// is waveform data.
		//Wave/Z MPStruct.xWave = $StringFromList(i, xwavelist)

		// just in case the new wave has a different number of points
		MPStruct.XPointRangeEnd = numpnts(MPStruct.yWave)-1
		
		err = MPF2_DoMPFit(MPStruct, BatchDFName+":")
		//err = MPF2_DoMPFit(MPStruct, datafolderName)
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
			alertMsg += StringFromList(i, ywavelist)+": "
			alertMsg += num2str(MPStruct.fitError)
			alertMsg += "; continuing with next fit."
			DoAlert 0, alertMsg
		endif
		
		String s=BatchDFName + ":'Peak 0 Coefs'"
		//String s=datafolderName + ":'Peak 0 Coefs'"
		Wave p0 = $s
		
		pulseHeights[i]=p0[2]
		//pulseAreas[i]=area(p0)
		pulseNames[i]=StringFromList(i,ywavelist)
		SetDataFolder saveDF
	endfor
end


//Makes a list of data traces to put into ywavelist above.
//Assumes data is called pulse# where # is consecutive
//and stored in a numbered folder inside of a folder called $prefix.
//Example: $prefix = "data"
//the data folder contains folders named "1", "2", "3", etc...
//Each numbered folder contains X traces with names "pulse1", "pulse2", etc.
//This is stupid. Put the data in a matrix and use MultiFitMatrix.ipf
Function/S MakeList(prefix, numFolders, numPulses)

	String prefix
	Variable numFolders, numPulses
	
	Variable i,j
	
	String outputList = ""
	
	
	for (i=1; i<=numFolders; i+=1)
		for (j=1;j<=numPulses;j+=1)
			outputList = outputList + "root:"+prefix +":'"+ num2str(i) +"':pulse"+num2str(j)+ ";"
		endfor
	endfor
	
	return outputList
End