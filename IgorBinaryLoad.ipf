#pragma rtGlobals=3		// Use modern global access method and strict wave access.


//Pulls data out of one specific folder, for one specific trace
//folder variable is interna folder to store data in eg. "root:data:"
//numTrace variable is the waveform number (ie scope channel)
function LoadFolderBin(folder, numTrace)
	
	String folder
	Variable numTrace
	String oldFolder
	
	oldFolder = GetDataFolder(1)
	
		
	NewDataFolder/O $folder
	SetDataFolder $folder
	
	variable i=0
	string wname, fname, temp
	

	getfilefolderinfo/D
	

	
	string baseFolder=S_path

		
		NewPath/O/Q dataFolder, baseFolder
		
		i=0
		
		string filelist=IndexedFile(dataFolder, -1, ".bin")
		
		do
			fname=baseFolder+StringFromList(i, filelist)
				
				LoadAgilentBin(fname, numTrace)
				Wave w= voltageVector
				Wave wTime=timeWave
				
				//Uncomment this line and comment the following to import a subrange
				//Duplicate/O/FREE/R=[24000,28000] w, shortVector
				Duplicate/O/FREE w, shortVector
				
				Concatenate {shortVector}, allData
				Concatenate/NP {wTime}, cDates
				
			//endif
			
			i+=1
			
			if(mod(i, 1000)==0)
				print i
			endif
			
		//while(i<ItemsInList(filelist))
		while(i<1000)
	KillWaves w, wTime
	
	SetDataFolder oldFolder
end

//Same as load folder, except grabs data from all four channels
function LoadFolderBinAll()
	
	String folder="data"
	//Variable numTrace
	String oldFolder
	
	oldFolder = GetDataFolder(1)
	
		
	NewDataFolder/O $folder
	SetDataFolder $folder
	
	variable i=0
	variable counter=0
	string wname, fname, temp
	

	getfilefolderinfo/D
	

	
	string baseFolder=S_path

		
		NewPath/O/Q dataFolder, baseFolder
		
		i=0
		
		string filelist=IndexedFile(dataFolder, -1, ".bin")
		
		do
			fname=baseFolder+StringFromList(i, filelist)
			GetFileFolderInfo/Q fname
			Wave wTime=timeWave
			
				
				LoadAgilentBin(fname, 1)
				Wave w= voltageVector
				
				Concatenate {w}, ptData
				
				LoadAgilentBin(fname, 2)
				Wave w= voltageVector
				
				Concatenate {w}, tiData
				
				LoadAgilentBin(fname, 3)
				Wave w= voltageVector
				
				Concatenate {w}, auData
				
				LoadAgilentBin(fname, 4)
				Wave w= voltageVector
				
				Concatenate {w}, pmtData
				Concatenate/NP {wTime}, cDates

			
			i+=1
			
			if(mod(i, 1000)==0)
				print i
			endif
			
		while(i<ItemsInList(filelist))
		
	KillWaves w, wTIme
	
	SetDataFolder oldFolder
end

//This will pull all the channel X data from many files that are stored in consecutively numbered folders.
//numTrace is the scope channel number. isDir is 1 if you've already called GetFileFolderInfo
function LoadTracesBin(folder, numFolders, numTrace, isDir)
	
	String folder
	Variable numFolders, numTrace, isDir
	String oldFolder
	
	oldFolder = GetDataFolder(1)
	
		
	NewDataFolder/O $folder
	SetDataFolder $folder
	
	variable i=0
	variable j=1
	string wname, fname, temp
	
	if (isDir==0)
		getfilefolderinfo/D
	endif
	
	string baseFolder=S_path
	string currentFolder
	
	for(j=1;j<=numFolders;j+=1)
		currentFolder=baseFolder+num2str(j)+":"
		
		NewPath/O/Q dataFolder, currentFolder
		
		i=0
		
		string filelist=IndexedFile(dataFolder, -1, ".bin")
		
		do
			fname=currentFolder+StringFromList(i, filelist)
			
			LoadAgilentBin(fname, numTrace)
			Wave w=voltageVector
			Concatenate {w}, allData
			
			i+=1
		while(i<ItemsInList(filelist))
		print j
	endfor
	
	KillWaves w
	
	scaleVolts(allData, numTrace)
	
	print "All Done"
	
	SetDataFolder oldFolder
end


//Basic procedure to load a file and pull one channel's worth of data out of it. Called by the above procedures.
function LoadAgilentBin(fileName, waveformSelect)
	
	String fileName
	Variable waveformSelect
	Variable fileId
	
	
	//Set up all the variable names along with a summy string variable
	String fileCookie, fileVersion, dateString, timeSTring, frameSTring, waveformString, str="0"
	Variable fileSize, nWaveforms, headerSize, waveformType, nWaveformBuffers, nPoints
	Variable count, xDisplayRange, xDisplayOrigin, xIncrement, xOrigin, xUnits, yUnits, timeTag, segmentIndex
	Variable bytesLeft, bufferType, bufferSize, bytesPerPoint

	
	//Igor reads until the string runs out of memory when it reads from a binary file, so set up all the strings to the correct size
	fileCookie=PadString(str, 2, 0)
	fileVersion=PadString(str, 2, 0)
	dateString=PadString(str, 16, 0)
	timeString=PadString(str, 16, 0)
	frameString=PadString(str, 24, 0)
	waveformString=PadString(str, 16, 0)

	
	//Open the file for reading
	Open/R fileID as fileName
	
	//Read some basics about the data in the file. Check to make sure it is good
	FBinRead fileID, fileCookie
	FBinRead fileID, fileVersion
	FBinRead/F=3 fileID, fileSize
	FBinRead/F=3 fileID, nWaveforms
	
	//Spit out how many waveforms live in the file
	//print nWaveforms
	
	Variable waveformIndex, bufferIndex, i
	
	//Check to make sure it's the right kind of file
	if(cmpstr(fileCookie, "AG")==0)
		//Read in all the header data
		for(waveformIndex=1;waveformIndex<=nWaveforms; waveformIndex+=1)
			FBinRead/F=3 fileID, headerSize
			bytesLeft=headerSize-4
			FBinRead/F=3 fileID, waveformType
			bytesLeft-=4
			FBinRead/F=3 fileID, nWaveformBuffers
			bytesLeft-=4
			FBinRead/F=3 fileID, nPoints
			bytesLeft-=4
			FBinRead/F=3 fileID, count
			bytesLeft-=4
			FBinRead/F=4 fileID, xDisplayRange
			bytesLeft-=4
			FBinRead/F=5 fileID, xDisplayOrigin
			bytesLeft-=8
			FBinRead/F=5 fileID, xIncrement
			bytesLeft-=8
			FBinRead/F=5 fileID, xOrigin
			bytesLeft-=8
			FBinRead/F=3 fileID, xUnits
			bytesLeft-=4
			FBinRead/F=3 fileID, yUnits
			bytesLeft-=4
			FBinRead fileID, dateString
			bytesLeft-=16
			FBinRead fileID, timeString
			bytesLeft-=16
			FBinRead fileID, frameString
			bytesLeft-=24
			FBinRead fileID, waveformString
			bytesLeft-=16
			FBinRead/F=5 fileID, timeTag
			bytesLeft-=8
			FBinRead/U/F=3 fileID, segmentIndex
			bytesLeft-=4
			
			//skip all remaining header data (if any)
			FStatus fileID
			FSetPos fileID, V_filePos+bytesLeft


			//generate voltageVector and set scaling from xIncrement and xOrigin values
    			if (waveformIndex == waveformSelect)
        			//Don't know why there would be more than one waveform buffer, but whatever...
        			if(nWaveformBuffers>1)
        				Make/O/N=(nPoints,nWaveformBuffers) voltageVector
        			else
        				Make/O/N=(nPoints) voltageVector
        			endif
        			SetScale/P x, xOrigin, xIncrement, voltageVector
   			endif

    			for(bufferIndex = 1;bufferindex<=nWaveformBuffers;bufferIndex+=1)
       		// read waveform buffer header
        			FBinRead/F=3 fileID, headerSize
        			bytesLeft = headerSize - 4
        			FBinRead/F=2 fileID, bufferType
        			bytesLeft = bytesLeft - 2
        			FBinRead/F=2 fileID, bytesPerPoint
        			bytesLeft = bytesLeft - 2
        			FBinRead/F=3 fileID, bufferSize
        			bytesLeft = bytesLeft - 4

        			// skip over any remaining data in the header
       			FStatus fileID
				FSetPos fileID, V_filePos+bytesLeft
				
				Variable temp
				
				//This whole mess just outputs the file creation time to the variable igorDate
        			if (waveformIndex == waveformSelect)
        				Variable day, month, year, hrs, mins, secs
					String monthS
					sscanf dateString, "%d %s %d", day, monthS, year
					sscanf timeString, "%d:%d:%d", hrs, mins, secs
					
					strswitch(monthS)
						case "JAN":
							month=1
							break
						case "FEB":
							month=2
							break
						case "MAR":
							month=3
							break
						case "APR":
							month=4
							break
						case "MAY":
							month=5
							break
						case "JUN":
							month=6
							break
						case "JUL":
							month=7
							break
						case "AUG":
							month=8
							break
						case "SEP":
							month=9
							break
						case "OCT":
							month=10
							break
						case "NOV":
							month=11
							break
						case "DEC":
							month=12
							break
					endswitch
					
					Variable igorDate
					igorDate=date2secs(year, month, day)+secs+mins*60+hrs*60*60
					Make/O/N=1 timeWave=igorDate
					//End of code that gets the time

            				if ((bufferType == 1) | (bufferType == 2) | (bufferType == 3))
                				// bufferType is PB_DATA_NORMAL, PB_DATA_MIN, or PB_DATA_MAX (float)
                				for(i=0;i<nPoints;i+=1)
                					FBinRead/F=4 fileID, temp
                					if(nWaveformBuffers>1)
                						voltageVector[i][bufferIndex]=temp
                					else
                						voltageVector[i]=temp
                					endif
                				endfor
            				elseif (bufferType == 4)
                				// bufferType is PB_DATA_COUNTS (int32)
                				for(i=0;i<nPoints;i+=1)
                					FBinRead/F=3 fileID, temp
                					if(nWaveformBuffers>1)
                						voltageVector[i][bufferIndex]=temp
                					else
                						voltageVector[i]=temp
                					endif
                				endfor
            				elseif (bufferType == 5)
                			// bufferType is PB_DATA_LOGIC (int8)
                				for(i=0;i<nPoints;i+=1)
                					FBinRead/U/F=1 fileID, temp
                					if(nWaveformBuffers>1)
                						voltageVector[i][bufferIndex]=temp
                					else
                						voltageVector[i]=temp
                					endif
                				endfor
					else
                			// unrecognized bufferType read as unformated bytes
                				for(i=0;i<bufferSize;i+=1)
                					FBinRead/U/F=1 fileID,temp
                					if(nWaveformBuffers>1)
                						voltageVector[i][bufferIndex]=temp
                					else
                						voltageVector[i]=temp
                					endif
                				endfor
            				endif
        			else
        				//this is not the data you are looking for, move to the end of the buffer
        				FStatus fileID
					FSetPos fileID, V_filePos+bufferSize
        			endif
   			 endfor
		endfor
		//Close the file up
		FStatus fileID
		FSetPos fileID, V_logEOF
		Close fileID

	else
		print fileCookie	
		FStatus fileID
		FSetPos fileID, V_logEOF
		Close fileID
	endif
	
	
end	

//Scale channel from scope voltage to TES current
//Can pass a zero for no scaling

function scaleVolts(dataMat, chan)
wave dataMat
variable chan

//Added some gain recently, need to take care of that
dataMat/=5

if (chan==1)
	dataMat*=-2/245455
elseif (chan==2)
	dataMat*=-2/242000
elseif (chan==3)
	dataMat*=-2/257800
else
	dataMat*=1
endif

end

function LoadIbIs(chan, temps, currs)
	Variable chan
	Wave temps, currs
	Variable temp, biasCurrent, i

	GetFileFolderInfo/D
	
	string baseFolder=S_path
	string fileName
	
	For(i=0;i<DimSize(temps, 0);i+=1)
	
		temp=temps[i]
		biasCurrent=currs[i]
	
		fileName=S_path+num2str(temps[i])+"mK:CH"+num2str(chan)+"_"+num2str(currs[i])+"uA_50OhmInput_Gain1_10kHzLPF_256Avg.bin"
	
		LoadAgilentBin(fileName, 1)
		
		Wave voltageVector
		String wName
	
		Duplicate/FREE voltageVector, Is
	
		scaleVolts(Is, chan)
	
		Is*=5 //Remove the gain by 5, since this is IbIs curve
	
		LoadAgilentBin(fileName, 2)
	
	
		Duplicate/FREE voltageVector, Ib
	
		WaveStats/Q Ib
	
		if(V_minRowLoc<V_maxRowLoc)
		WaveStats/Q/R=[0, V_minRowLoc] Ib
		endif

		wName="Ib"+num2str(temp)
	
		Redimension/N=(V_minRowLoc-V_maxRowLoc+1) Ib
		
		Ib=(1e-6)*biasCurrent*((DimSize(Ib,0)-p)/DimSize(Ib,0))
		
		Duplicate/O Ib, $wName
		
		wName="Is"+num2str(temp)
	
		Duplicate/O/R=[V_maxRowLoc, V_minRowLoc] Is, $wName
	EndFor
	
	KillWaves voltageVector
	
end
	
