#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#include <XY Pair To Waveform Panel>
#include <Wave Arithmetic Panel>
#include <Multi-peak fitting 2.0>
#include <Waves Average>
#include <Function Grapher>
#include <Autosize Images>

//INDEX (not necessarily in order)

//***For TES stuff
//scaleArea(dataVec, chan)
//alignIbIs(Is, Ib, headCoef,tailCoef)

//***Viewing histograms in different ways
//peakAverage(dataMat,refHist,low, high)
//peakForest(dataMat, refHist, startX, endX, dX)

//***Making histograms in different ways
//pulseInt(dataMat, endtime, chan)
//dotAvg(dataMat, avgWave)
//pulseHeights(dataMat)
//pulseHeightsInt(dataMat2, endVal)

//***Specific functions for two-pixel TES
//threshCheck(dataMat, level, startX, endX)
//makeBanana(dataMat1, dataMat2, mask1, mask2, times1, times2)
//pulseHeightsIntBanana(dataMat2, dataMat3, endVal)

//***Some experimental optimal filter crap
//covMat(dataMat)
//wdotAvg(dataMat, wMat, avgWave)
//pulseCorrelate(dataMat, testWave)
//matCorrelate(dataMat1, dataMat2)

//***Some routines for dealing with basic matrices full of time traces
//filterPulses(dataMat)
//subtractWave(dataMat, subWave)
//removeBaseline(dataMat, endpoint)
//avgWaves(dataMat)
//Disp(dataMat, n)
//removeConditional(dataMat1, dataMat2)

//Poisson Fit Function


//This function loads data, processes it, and makes some plots
function LoadAndProcess(folder, numTrace)

	String folder
	Variable numTrace

	//Get channel data and scale it
	print "Loading data"
	LoadFolderBin(folder, numTrace)
	
	//Step into data folder
	string oldFolder = GetDataFolder(1)
	SetDataFolder $folder
	
	//Scale to current from voltage
	print "Scaling Data"
	Wave w = allData
	ScaleVolts(w, numTrace)
	
	w*=-1
	
	//Align the pre-trigger noise to zero average
	print "Removing Baseline"
	removeBaseline(w, 200)
	
	//Make an average of the data
	print "Averaging Waves"
	avgWaves(w)
	
	//Filter data with 3MHz LPF
	String filt_folder = oldFolder+folder + ":filtered"
	
	NewDataFolder/O $filt_folder
	SetDataFolder $filt_folder
	
	Duplicate/O w, allDataF
	
	print "Filtering Data"
	filterPulses(allDataF, 50e6, 3e6)
	
	//Create some useful quantities
	print "Making analysis waves"
	pulseHeightsInt(allDataF, 0.5, numTrace)
	riseTimesHist(allDataF)
	
	SetDataFolder oldFolder
end

//This takes scaled IbIs curves and computes voltage and resistance and power
function RandV(Is, Ib, Rp, Rs, chan)
	Wave Is, Ib
	Variable Rp, Rs, chan
	
	if (chan==1)
		Duplicate/O Is, ch1Rtes
		Duplicate/O Is, ch1Vtes
		Duplicate/O Is, ch1Ptes
		
		ch1Rtes=Rs*(Ib/Is-1)-Rp
		ch1Vtes=ch1Rtes*Is
		ch1Ptes=ch1Vtes*Is
		
	elseif (chan==2)
		Duplicate/O Is, ch2Rtes
		Duplicate/O Is, ch2Vtes
		Duplicate/O Is, ch2Ptes
			
		ch2Rtes=Rs*(Ib/Is-1)-Rp
		ch2Vtes=ch2Rtes*Is
		ch2Ptes=ch2Vtes*Is
	elseif (chan==3)
		Duplicate/O Is, ch3Rtes
		Duplicate/O Is, ch3Vtes
		Duplicate/O Is, ch3Ptes
		
		ch3Rtes=Rs*(Ib/Is-1)-Rp
		ch3Vtes=ch3Rtes*Is
		ch3Ptes=ch3Vtes*Is
	endif
end

//Scale Areas from coulombs to eV
function scaleArea(dataVec, chan)
wave dataVec
variable chan



if (chan==1)
	//reflecto-device
	//dataVec*=5.8e-6/1.6e-19
	
	//Shielded device (March) at 117 mK
	dataVec*=3e-6/1.6e-19
elseif (chan==2)
	
	//reflecto-device
	//dataVec*=4.58e-6/1.6e-19
	
	//Shielded device (March) at 117 mK
	dataVec*=3.25e-6/1.6e-19
	
elseif (chan==3)
	//reflecto-device
	//dataVec*=4.37e-6/1.6e-19
	
	//Shielded device (March) at 117 mK
	dataVec*=3.55e-6/1.6e-19
else
	print "Must choose channel 1-3"
endif

end

//Convert a scatter plot to a 2-D density histogram
function scatter2density(xWav, yWav, logscalex, nBins, xMin, xMax, yMin, yMax)
wave xWav, yWav
variable logscalex, nBins

Variable xMin, xMax, yMin, yMax
variable xBin, yBin

WaveStats/Q xWav

//xMax=V_max
//xMin=V_min
//
//WaveStats/Q yWav
//
//yMax=V_max
//yMin=V_min

Make/O/N=(nBins, nBins) densityMatrix=0
Make/O/N=(nBins+1) densityXscale
Make/O/N=(nBins+1) densityYscale

densityXscale=xMin+((xMax-xMin)/(1+nBins))*p

variable logMax, logMin

if(logscalex==1)
	logMax=log(yMax)
	logMin=log(yMin)
	densityYscale=10^(logMin+((logMax-logMin)/(1+nBins))*p)
else
	densityYscale=yMin+((yMax-yMin)/(1+nBins))*p
endif

variable i,xPos, yPos

for(i=0; i<DimSize(xWav, 0);i+=1)
	
	xPos=trunc((xWav[i]-xMin)*(nBins)/(xMax-xMin))
	
	if(logscalex==1)
		yPos=trunc((log(yWav[i])-logMin)*(nBins)/(logMax-logMin))
	else
		yPos=trunc((yWav[i]-yMin)*(nBins)/(yMax-yMin))
	endif
	
	if(numtype(xPos)!=2 && numtype(yPos)!=2)
		if(xPos==nBins)
			xPos=nBins-1
		endif
		if(yPos==nBins)
			yPos=nBins-1
		endif
		
		if((xPos>nBins)||(yPos>nBins)||(xPos<0)||(yPos<0))
		else
			densityMatrix[xPos][yPos]+=1
		endif
	endif
endfor

end




//Slice up a histogram and make pulse averages of each slice
function peakForest(dataMat, refHist, startX, endX, dX)
wave dataMat, refHist
variable startX, endX, dX

variable i, j, k

string newName




Make/FREE/O/N=(DimSize(dataMat,0)) temp

Display

for(k=startX; k<=endX; k+=dX)
	
	j=0
	Make/O/N=(DimSize(dataMat,0)) peakAvg
	CopyScales dataMat, peakAvg

	for(i=0;i<DimSize(dataMat,1);i+=1)
		if(refHist[i]>=k && refHist[i]<(k+dX))
			temp=dataMat[p][i]
			CopyScales dataMat, temp
			peakAvg=peakAvg+temp
			j+=1;
		endif
	endfor

	peakAvg/=j
	newName = "peakAvg" + num2str(k)
	Rename peakAvg, $newName
	
	AppendToGraph $newName
	
endfor

end
	
//Remove some waves based on some condtion
function removeConditional(dataMat1)
wave dataMat1
variable i, flag

Make/FREE/N=(DimSize(dataMat1,0)) temp

i=0
do
	temp=dataMat1[p][i]
	CopyScales dataMat1, temp
	
	WaveStats/Q/R=(-0.5e-6, 0.5e-6) temp
	if (V_min <= -25e-9)
		flag=1
	endif
	
	if(flag==1)
		DeletePoints/M=1 i, 1, dataMat1
		flag=0
	else
		i+=1
	endif
	
	
while(i<DimSize(dataMat1, 1))
end

//Looks for waves over a threshold. Creates two waves: a mask wave that is 0 or 1 for each wave depending
//on whether or not the threshold was crossed. The other wave is the time value the threshold was crossed.
function threshCheck(dataMat, level, startX, endX)
wave dataMat
variable level, startX, endX
variable i
string name

name=NameOfWave(dataMat)+"_mask"

Make/O/N=(DimSize(dataMat, 1)) $name
wave mask = $name

name=NameOfWave(dataMat)+"_times"
Make/O/N=(DimSize(dataMat, 1)) $name
wave times = $name

Make/FREE/N=(DimSize(dataMat,0)) temp

for(i=0;i<DimSize(dataMat,1);i+=1)
	
	temp=dataMat[p][i]
	CopyScales dataMat, temp
	
	FindLevel/B=3/EDGE=1/R=(startX, endX)/Q temp,level
	
	mask[i]=V_flag
	times[i]=V_LevelX
endfor

end

//Take the data and masks from threshCheck and make a bananaPlot
function makeBanana(dataMat1, dataMat2, mask1, mask2, times1, times2)
wave dataMat1, dataMat2
wave mask1, mask2, times1, times2
variable size, i

size=0

for(i=0;i<DimSize(dataMat1,1);i+=1)
	if(mask1[i]==0 || mask2[i]==0)
		size+=1
	endif
endfor

Make/O/N=(size) bLeft
Make/O/N=(size) bRight
//Make/O/N=(size) tLeft
//Make/O/N=(size) tRight
Make/O/N=(size) deltaT


Make/FREE/N=(DimSize(dataMat1,0)) temp1
Make/FREE/N=(DimSize(dataMat2,0)) temp2


for(i=0;i<size;i+=1)
	if(mask1[i]==0 || mask2[i]==0)
		temp1=dataMat1[p][i]
		CopyScales dataMat1, temp1
		
		temp2=dataMat2[p][i]
		CopyScales dataMat2, temp2
		
		//bLeft[i]=area(temp1, 0, +Inf)
		//bRight[i]=area(temp2, 0, +inf)
		
		Smooth/B=5 temp1
		Smooth/B=5 temp2
		WaveStats/Q/R=(0,2e-6) temp1
		deltaT[i]=V_maxloc
		bLeft[i]=V_max
		WaveStats/Q/R=(0,2e-6) temp2
		deltaT[i]-=V_maxloc
		bRight[i]=V_max
		
	endif
endfor

end

//Align an IbIs curve
//First do a fit to each linear region, then pass the coeficient matricies to this function
function alignIbIs(Is, Ib, headCoef,tailCoef)
wave Is, Ib, headCoef, tailCoef
variable newX, newY

newX=(headCoef[0]-tailCoef[0])/(tailCoef[1]-headCoef[1])

newY=(headCoef[0]*tailCoef[1]-tailCoef[0]*headCoef[1])/(tailCoef[1]-headCoef[1])

Is-=newY
Ib-=newX

printf "Parasitic resistance is: %g Ohms\r", ((1/headCoef[1]-1)*0.2)




end

//Make a covariance matrix
function covMat(dataMat)
wave dataMat


MatrixOp/O allData_Cov =  (subtractMean(dataMat,2) x subtractMean(dataMat,2)^t)
allData_Cov/=(DimSize(dataMat,0)-1)
//MatrixOp/O allData_W = Inv(allData_Cov)


end

//Average some subset of waves together
function peakAverage(dataMat,refHist,low, high)
wave dataMat, refHist
variable low, high
variable i, j

j=0


Make/O/N=(DimSize(dataMat,0)) temp

Make/O/N=(DimSize(dataMat,0)) peakAvg
CopyScales dataMat, peakAvg

for(i=0;i<DimSize(dataMat,1);i+=1)

	if(refHist[i]>=low && refHist[i]<=high)
		temp=dataMat[p][i]
		CopyScales dataMat, temp
		peakAvg=peakAvg+temp
		j+=1;
	endif
endfor

peakAvg/=j

KillWaves temp

end

//Make a sub-matrix using a mask to pull out certain columns
function makeSubMatrix(dataMat, mask)
wave dataMat, mask
variable i, length

length=DimSize(dataMat, 0)

Wave subMatrix
If(WaveExists(subMatrix))
	Rename subMatrix, subMatrix_old
Endif

For(i=0;i<(area(mask));i+=1)
	Make/O/FREE/N=(length) temp



	If(mask[i]==1)
		temp=dataMat[p][i]
		CopyScales dataMat, temp
		Concatenate {temp}, subMatrix
	EndIf

EndFor

End

//Find the first 10 PMT pulses for each trace
Function findPMT(dataMat)
Wave dataMat

Duplicate/O/FREE dataMat, dataMat_smth_diff
Smooth/S=4 55, dataMat_smth_diff
Differentiate/DIM=0 dataMat_smth_diff

Make/FREE/N=(DimSize(dataMat_smth_diff,0)) temp
CopyScales dataMat_smth_diff, temp

Variable i,j

Make/O/N=(10, (DimSize(dataMat, 1))) pmtLocs = nan

For(i=0; i<(DimSize(dataMat_smth_diff,1)); i+=1)
	temp=dataMat_smth_diff[p][i]
	
	Wave levels
	FindLevels/D=levels/N=10/EDGE=2/Q temp, -25e3
	
	For(j=0; j<DimSize(levels, 0);j+=1)
		pmtLocs[j][i]=levels[j]
	EndFor
	
	//Concatenate {levels}, pmtLocs
	
EndFor

End

//Make a mask based on the PMT histogram
//Can also include an optional mask to further refinet the PMT mask
Function makePMTMask(dataMat, x1, x2, [mask])
Wave dataMat
Variable x1, x2
wave mask

Variable i, j

if(ParamIsDefault(mask))


	
	Make/O/N=(DimSize(dataMat, 0)) pmtmask=0
	
	For(i=0; i<DimSize(dataMat, 0); i+=1)
		For(j=0; j<DimSize(dataMat, 1); j+=1)
			If((dataMat[i][j]<x2)&&(dataMat[i][j]>x1))
				pmtmask[i]=1
			EndIf
		EndFor
	EndFor
Else
	
	Make/O/N=(DimSize(dataMat, 0)) pmtmask=0
	
	For(i=0; i<DimSize(dataMat, 0); i+=1)
		For(j=0; j<DimSize(dataMat, 1); j+=1)
			If((dataMat[i][j]<x2)&&(dataMat[i][j]>x1)&&(mask[i]==1))
				pmtmask[i]=1
			EndIf
		EndFor
	EndFor
EndIf
End
		


//Make a PSD from a dataset and a mask wave
//The mask contains the locations of the noise traces for the psd
//Might want to try putting together an average psd for the model too...
function makePSD(dataMat, mask0)
wave dataMat, mask0
variable i, length

If(mod(DimSize(dataMat, 0),2)==1)
	length=DimSize(dataMat,0)+1
Else
	length=DimSize(dataMat, 0)
EndIf

Wave tempPSD
If(WaveExists(tempPSD))
	Rename tempPSD temp_PSD_Old
Endif

For(i=0; i<(area(mask0)); i+=1)
Make/O/FREE/N=(DimSize(dataMat, 0)) temp
	
	If(mask0[i]==1)
		temp=dataMat[p][i]
		CopyScales dataMat, temp
		
		//Not going to bother converting to per Hz. Leaving in W.
		Wave/C W_FFT
		FFT/PAD=(length)/DEST=W_FFT temp
		W_FFT=W_FFT*conj(W_FFT)
		//W_FFT/=(length/2)
		//W_FFT[0]/=2
		//W_FFT[length/2]/=2
		Concatenate/C {W_FFT}, tempPSD
		KillWaves W_FFT
	EndIf
EndFor

MatrixOp/C/O/FREE noiseFFT2=sumRows(tempPSD)/numCols(tempPSD)
Redimension/C/N=(length/2+1) noiseFFT2
CopyScales tempPSD, noiseFFT2
KillWaves tempPSD


Make/N=(length/2+1)/O noisePSD
CopyScales noiseFFT2, noisePSD
noisePSD = real(noiseFFT2)

Killwaves noiseFFT2
End

//Takes a noise PSD, a pulse model, and a dataset, and applies the filter
function applyFilter(dataMat, filter)
Wave dataMat, filter

Variable length=DimSize(dataMat, 0)

if(mod(length, 2)==1)
	length+=1
endif

Wave/C W_FFT
FFT/PAD=(length)/DEST=W_FFT dataMat

W_FFT*=filter[p]

Wave filteredData
IFFT/DEST=filteredData W_FFT

KillWaves W_FFT


End

//Average some subset of waves together
function peakAverageM(dataMat, mask)
wave dataMat, mask
variable i, j

j=0


Make/O/N=(DimSize(dataMat,0)) temp

Make/O/N=(DimSize(dataMat,0)) peakAvg
CopyScales dataMat, peakAvg

for(i=0;i<DimSize(dataMat,1);i+=1)

	if(mask[i]==1)
		temp=dataMat[p][i]
		CopyScales dataMat, temp
		peakAvg=peakAvg+temp
		j+=1;
	endif
endfor

peakAvg/=j

KillWaves temp

end

//Integrate each wave and add to histogram
//chan sets the scaling to eV
function pulseInt(dataMat, starttime, endtime, chan)
wave dataMat
variable starttime, endtime, chan
variable i

Make/O/FREE/N=(DimSize(dataMat,0)) temp

Make/O/N=(DimSize(dataMat,1)) areas

for(i=0;i<DimSize(dataMat,1);i+=1)
	temp=dataMat[p][i]
	CopyScales dataMat, temp
	areas[i]=area(temp, starttime,endtime)
endfor

//ScaleArea(areas, chan)

end

//Dot each wave with average and add to histogram
function dotAvg(dataMat, avgWave)
wave dataMat, avgWave

Duplicate/FREE/R=() avgWave, avW
Duplicate/FREE/R=()() dataMat, dtM

MatrixOp/FREE n=(avW . avW)
MatrixOp/O dots = (avW^t x dtM)^t
dots/=n[0]

end

//Dot each wave with average subject to a weigt matrix and add to histogram
function wdotAvg(dataMat, wMat, avgWave)
wave dataMat, wMat, avgWave

MatrixOp/FREE n=(avgWave . (wMat x avgWave))
MatrixOp/O dots_W = (((wMat x avgWave)^t x dataMat)^t)

dots_W/=n[0]

end

//Take max height in some range and add to histogram
function pulseHeights(dataMat, first, last)
wave dataMat
variable first, last
variable i

Make/O/N=(DimSize(dataMat,0)) temp

Make/O/N=(DimSize(dataMat,1)) heights
Make/O/N=(DimSize(dataMat,1)) times

for(i=0;i<DimSize(dataMat,1);i+=1)
	temp = dataMat[p][i]
	CopyScales dataMat, temp
	
	WaveStats/Q/R=(first,last) temp
	
	heights[i]=V_max
	times[i]= V_maxloc
endfor

KillWaves temp

end

//Try and get the steepness of a trace
function riseTimesHist(dataMat)
wave dataMat
variable i, length, numTraces
string name

length=DimSize(dataMat, 0)
numTraces=DimSize(dataMat, 1)

name=NameOfWave(dataMat)+"_riseTimes"
Make/O/N=(numTraces) $name
Wave wRiseTimes = $name

name=NameOfWave(dataMat)+"_riseLocs"
Make/O/N=(numTraces) $name
Wave wRiseLocs = $name

Make/O/FREE/N=(length) temp

for(i=0; i<numTraces;i+=1)
	temp = dataMat[p][i]
	CopyScales dataMat, temp
	
	Smooth/B 5, temp
	Differentiate temp
	
	WaveStats/Q temp
	wRiseTimes[i]=V_max
	wRiseLocs[i]=V_maxloc
	
endfor

end

//Calibrate model pulse by applying filter and then running the width/area algorithm
function pulseCalibrate(model, filter, intCutoff, chan)
Wave model, filter
Variable intCutoff, chan

Variable length=DimSize(model, 0)
If(mod(length, 2)==1)
	length+=1
EndIf

	Wave/C model_FFT
	FFT/DEST=model_FFT/PAD=(length) model
	
	model_FFT*=filter
	
	Wave model_filtered
	IFFT/DEST=model_filtered model_FFT

	
	Duplicate/O model_filtered, temp, tempF1, tempF3, tempMax, tempThresh

	Variable fs, ff	
	fs=50e6
	ff=3e6
	
	FilterFIR/DIM=0/LO={ff/fs, ff*1.2/fs,101}/WINF=Hamming tempF3

	WaveStats/Q tempF3
	
	Variable height, hloc, t1, t2
	height=V_max
	hloc= V_maxloc
	//wMins[i]=V_min
	//wMinlocs[i]=V_minloc
	
	tempMax=V_max
	tempThresh=V_max*intCutoff
	

	
	FindLevel/Q/R=(V_maxloc, -inf) tempF3, (V_max*intCutOff)
	
	if(V_flag==0)
		t1=V_LevelX
	else
		t1=pnt2x(temp, V_startRow)
	endif
	
	ff=0.5e6
	FilterFIR/DIM=0/LO={ff/fs, ff*1.2/fs,101}/WINF=Hamming tempF1

	
	FindLevel/Q/R=(V_maxloc, inf) tempF1, (V_max*intCutOff)
	
	if(V_flag==0)
		t2=V_LevelX
	else
		t2=pnt2x(temp, V_endRow)
	endif
	
		Duplicate/O/R=(t1, t2) temp, temp2
		Duplicate/O temp, tempt
		tempt=0
		variable index
		WaveStats/Q temp2
		for(index=V_startRow;index<V_endRow;index+=1)
			tempt[x2pnt(tempt, pnt2x(temp2, index))]=temp2[index]
		endfor
		Display/N=diagnostics temp, temp2, tempF3, tempF1, tempMax, tempThresh
		ModifyGraph lSize=1
		ModifyGraph lStyle(tempMax)=3
		ModifyGraph lStyle(tempThresh)=3
		ModifyGraph rgb(tempMax)=(0,0,0)
		ModifyGraph rgb(tempThresh)=(0,0,0)
		ModifyGraph mode(temp2)=7
		ModifyGraph hbfill(temp2)=5
		ModifyGraph rgb(tempF1)=(10000,50000,10000)
		ModifyGraph rgb(tempF3)=(10000,10000,50000)
		ModifyGraph expand=3
		DoUpdate/W=diagnostics
		
		Variable ar=area(temp, t1,t2)
		Variable width=t2-t1
		
		if (chan==1)
			//reflecto-device
			//dataVec*=5.8e-6/1.6e-19
			
			//Shielded device (March) at 117 mK
			ar*=3e-6/1.6e-19
		elseif (chan==2)
			
			//reflecto-device
			//dataVec*=4.58e-6/1.6e-19
			
			//Shielded device (March) at 117 mK
			ar*=3.25e-6/1.6e-19
			
		elseif (chan==3)
			//reflecto-device
			//dataVec*=4.37e-6/1.6e-19
			
			//Shielded device (March) at 117 mK
			ar*=3.55e-6/1.6e-19
		else
			print "Must choose channel 1-3"
		endif
		
		printf "Area =  %g\r", ar
		printf "Width = %g\r", width
	
End


//Finds max height in some range, then walks down from there to compute energy integral
//(startVal, endVal) is the time range to look for a maximum peak
//intCutOff is the Y-value to stop integrating at
//chan is the channel. Use this for scaling.
function pulseHeightsInt(dataMat, intCutOff, chan)
wave dataMat
variable intCutOff, chan
variable i, t1, t2, length, numTraces
string name
variable flag=1


length=DimSize(dataMat, 0)
numTraces=DimSize(dataMat, 1)



Make/O/N=(length) temp

name=NameOfWave(dataMat)+"_heights"
Make/O/N=(numTraces) $name
Wave wHeights=$name

name=NameOfWave(dataMat)+"_times"
Make/O/N=(numTraces) $name
Wave wTimes=$name

name=NameOfWave(dataMat)+"_areas"
Make/O/N=(numTraces) $name
Wave wAreas=$name

name=NameOfWave(dataMat)+"_widths"
Make/O/N=(numTraces) $name
Wave wWidths=$name

//name=NameOfWave(dataMat)+"_spikiness"
//Make/O/N=(numTraces) $name
//Wave wSpikiness=$name
//
name=NameOfWave(dataMat)+"_mins"
Make/O/N=(numTraces) $name
Wave wMins=$name
//
name=NameOfWave(dataMat)+"_minLocs"
Make/O/N=(numTraces) $name
Wave wMinlocs=$name

for(i=0;i<DimSize(dataMat,1);i+=1)

//for(i=0;i<20;i+=1)
	temp = dataMat[p][i]
	CopyScales dataMat, temp
	
	Duplicate/O temp, tempF1, tempF3, tempMax, tempThresh

	Variable fs, ff	
	fs=50e6
	ff=3e6
	
	FilterFIR/DIM=0/LO={ff/fs, ff*1.2/fs,101}/WINF=Hamming tempF3

	WaveStats/Q tempF3
	
	wHeights[i]=V_max
	wTimes[i]= V_maxloc
	wMins[i]=V_min
	wMinlocs[i]=V_minloc
	
	tempMax=V_max
	tempThresh=V_max*intCutoff
	

	
	FindLevel/Q/R=(V_maxloc, -inf) tempF3, (V_max*intCutOff)
	
	if(V_flag==0)
		t1=V_LevelX
	else
		t1=pnt2x(temp, V_startRow)
	endif
	
	ff=0.5e6
	FilterFIR/DIM=0/LO={ff/fs, ff*1.2/fs,101}/WINF=Hamming tempF1

	
	FindLevel/Q/R=(V_maxloc, inf) tempF1, (V_max*intCutOff)
	
	if(V_flag==0)
		t2=V_LevelX
	else
		t2=pnt2x(temp, V_endRow)
	endif
	
	//if(t2-t1>40e-6)
		//flag=1
	//endif
	
	if(flag==1)
		
		Duplicate/O/R=(t1, t2) temp, temp2
		Duplicate/O temp, tempt
		tempt=0
		variable index
		WaveStats/Q temp2
		for(index=V_startRow;index<V_endRow;index+=1)
			tempt[x2pnt(tempt, pnt2x(temp2, index))]=temp2[index]
		endfor
		Display/N=diagnostics temp, temp2, tempF3, tempF1, tempMax, tempThresh
		ModifyGraph lSize=1
		ModifyGraph lStyle(tempMax)=3
		ModifyGraph lStyle(tempThresh)=3
		ModifyGraph rgb(tempMax)=(0,0,0)
		ModifyGraph rgb(tempThresh)=(0,0,0)
		ModifyGraph mode(temp2)=7
		ModifyGraph hbfill(temp2)=5
		ModifyGraph rgb(tempF1)=(10000,50000,10000)
		ModifyGraph rgb(tempF3)=(10000,10000,50000)
		ModifyGraph expand=3
		DoUpdate/W=diagnostics
		DoAlert 2, "t1="+num2str(t1)+"\r t2="+num2str(t2)
		flag=V_Flag
		KillWindow diagnostics
		if(flag==3)
			break
		endif
		//Concatenate {tempt}, allDataInts
	endif
	
	Variable ar=area(temp, t1,t2)
	
	if(ar>0)
		wAreas[i]=ar
		wWidths[i]=t2-t1
	else
		wAreas[i]=nan
		wWidths[i]=nan
	endif

endfor

scaleArea(wAreas, chan)

//wSpikiness= wHeights/wAreas

KillWaves temp, temp2, tempt, tempF1, tempF3, tempMax, tempThresh //, w, fitConstraints

end

//Find max height, then walk down from there and compute average
//This version does two pixels and computes some extra stuff
function pulseHeightsIntBanana(dataMat2, dataMat3, endVal)
wave dataMat2, dataMat3
variable endVal
variable i, t1, t2

Make/O/N=(DimSize(dataMat2,0)) temp
Make/O/N=(DimSize(dataMat2,1)) heights2
Make/O/N=(DimSize(dataMat2,1)) times2
Make/O/N=(DimSize(dataMat2,1)) hareas2
Make/O/N=(DimSize(dataMat2,1)) widths2

for(i=0;i<DimSize(dataMat2,1);i+=1)
	temp = dataMat2[p][i]
	CopyScales dataMat2, temp
	
	WaveStats/Q/R=(0,endVal) temp
	
	heights2[i]=V_max
	times2[i]= V_maxloc
	
	FindLevel/Q/R=(V_maxloc, 0) temp, (0.33*V_max)
	t1=V_LevelX
	
	FindLevel/Q/R=(V_maxloc, (V_maxloc+2*endVal)) temp, (0.33*V_max)
	t2=V_LevelX
	
	hareas2[i]=area(temp, t1, t2)

	widths2[i]=t2-t1
endfor

scaleArea(hareas2, 2)

//Now do the other channel

Make/O/N=(DimSize(dataMat3,0)) temp
Make/O/N=(DimSize(dataMat3,1)) heights3
Make/O/N=(DimSize(dataMat3,1)) times3
Make/O/N=(DimSize(dataMat3,1)) hareas3
Make/O/N=(DimSize(dataMat3,1)) widths3

for(i=0;i<DimSize(dataMat3,1);i+=1)
	temp = dataMat3[p][i]
	CopyScales dataMat3, temp
	
	WaveStats/Q/R=(0,endVal) temp
	
	heights3[i]=V_max
	times3[i]= V_maxloc
	
	FindLevel/Q/R=(V_maxloc, 0) temp, (0.33*V_max)
	t1=V_LevelX
	
	FindLevel/Q/R=(V_maxloc, (V_maxloc+2*endVal)) temp, (0.33*V_max)
	t2=V_LevelX
	
	hareas3[i]=area(temp, t1, t2)

	widths3[i]=t2-t1
endfor

scaleArea(hareas3, 3)

Duplicate/O times2, deltaT, deltaH, deltaW
deltaT=times2-times3
deltaH=heights2-heights3
deltaW=widths2-widths3

KillWaves temp

end

//Take correlation of some test wave and add it to histogram. Also find time difference.
function pulseCorrelate(dataMat, testWave)
wave dataMat, testWave
variable i

Make/O/N=(DimSize(dataMat,0)) temp

Make/O/N=(DimSize(dataMat,1)) Cheights
Make/O/N=(DimSize(dataMat,1)) Ctimes

for(i=0;i<DimSize(dataMat,1);i+=1)
	temp = dataMat[p][i]
	CopyScales dataMat, temp
	
	Duplicate/O/FREE temp, tempCorr
	
	Correlate testWave, tempCorr
	
	WaveStats/Q tempCorr
	
	Cheights[i]=V_max
	Ctimes[i]= V_maxloc
	
	if(mod(i, 1000)==0)
		print i
	endif
endfor

KillWaves temp

end

//Take correlation of each pair of pulses. Also find time difference. Assumes both matrices have same dimensions
function matCorrelate(dataMat1, dataMat2)
wave dataMat1, dataMat2
variable i
variable startTime

Make/O/FREE/N=(DimSize(dataMat1,0)) temp
CopyScales dataMat1, temp

Make/O/N=(DimSize(dataMat1,1)) Cheights
Make/O/N=(DimSize(dataMat1,1)) Ctimes

for(i=0;i<DimSize(dataMat1,1);i+=1)
	temp = dataMat1[p][i]
	startTime=pnt2x(temp,0)
	
	Duplicate/O/FREE temp, tempCorr
	tempCorr=dataMat2[p][i]
	
	Correlate temp, tempCorr
	
	WaveStats/Q tempCorr
	
	Cheights[i]=V_max
	Ctimes[i]= V_maxloc-startTime
	
	if(mod(i, 1000)==0)
		print i
	endif
endfor

end

//Apply a filter to the pulses (should probably remake the filter each time)
//fs is the sampling frequency
//ff is the LPF filter frequency
function filterPulses(dataMat, fs, ff)
wave dataMat
variable fs, ff
variable i, passBand, rejectBand

passBand=ff/fs
rejectBand=(ff+ff*0.2)/fs

Make/O/N=(DimSize(dataMat,0)) temp

for(i=0;i<DimSize(dataMat,1);i+=1)
	
	temp=dataMat[p][i]
	FilterFIR/DIM=0/LO={passBand, rejectBand,101}/WINF=Hamming temp
	dataMat[][i]=temp[p]

endfor

KillWaves temp

end

//subtracts some wave from the data Matrix
function subtractWave(dataMat, subWave)
wave dataMat, subWave
variable i

Make/O/N=(DimSize(dataMat,0)) temp

for(i=0;i<DimSize(dataMat,1);i+=1)
		
		temp=dataMat[p][i]
		CopyScales dataMat, temp
		temp= temp-subWave
		dataMat[][i]=temp[p]

endfor

KillWaves temp

end

//Removes the baseline from the data
function removeBaseline(dataMat, endpoint)
wave dataMat
variable endpoint
variable i

Make/O/N=(DimSize(dataMat,0)) baseline

for(i=0;i<DimSize(dataMat,1);i+=1)
		
		baseline=dataMat[p][i]
		WaveStats/Q/R=[0,endpoint] baseline
		dataMat[][i]-=V_avg
		
		if(mod(i,10000)==0)
			print i //helpful for big data sets so you know where you are
		endif

endfor

KillWaves baseline

end

//Averages all the columns of a Matrix together
//Also creates a wave containing the standard deviation
function avgWaves(dataMat)

wave dataMat
variable i
string name
variable length=DimSize(dataMat,1)

name=nameofwave(dataMat)+"_Avg"
Make/O/N=(DimSize(dataMat,0)) $name
wave wAvg = $name
CopyScales dataMat, wAvg

name=nameofwave(dataMat)+"_SDev"
Make/O/N=(DimSize(dataMat,0)) $name
wave wSDev = $name
CopyScales dataMat, wSDev


Make/O/N=(DimSize(dataMat,0)) temp
CopyScales dataMat, temp
temp=0

for(i=0;i<length;i+=1)

	temp+=dataMat[p][i]
	
endfor

Print "Avg Done"

wAvg=temp/length
temp=0

for(i=0;i<length;i+=1)

	temp+=(dataMat[p][i]-wAvg)^2
	
endfor

wSDev=sqrt(temp/length)

KillWaves temp



end

//Displays every nth wave in a matrix
function Disp(dataMat, n)
wave dataMat
variable n
variable i

Display

Make/O/N=(DimSize(dataMat,0)) temp

for(i=1;i<DimSize(dataMat,1);i+=1)

if(mod(i,n)==0)
temp=dataMat[p][i]

WaveStats/Q temp

//if(V_max > 0.006)

AppendToGraph dataMat[][i]

//Concatenate {temp}, bigPulses

//endif
endif

endfor

KillWaves temp
end

Function Poisson_x_scaled(w,x) : FitFunc
	Wave w
	Variable x

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ variable S0N=StatsPoissonPDF(0,mu)*StatsNormalPDF(x-x0,0*xs,sigma*xs)
	//CurveFitDialog/ variable S1N=StatsPoissonPDF(1,mu)*StatsNormalPDF(x-x0,1*xs,sigma*xs)
	//CurveFitDialog/ variable S2N=StatsPoissonPDF(2,mu)*StatsNormalPDF(x-x0,2*xs,sigma*xs)
	//CurveFitDialog/ variable S3N=StatsPoissonPDF(3,mu)*StatsNormalPDF(x-x0,3*xs,sigma*xs)
	//CurveFitDialog/ variable S4N=StatsPoissonPDF(4,mu)*StatsNormalPDF(x-x0,4*xs,sigma*xs)
	//CurveFitDialog/ variable S5N=StatsPoissonPDF(5,mu)*StatsNormalPDF(x-x0,5*xs,sigma*xs)
	//CurveFitDialog/ variable S6N=StatsPoissonPDF(6,mu)*StatsNormalPDF(x-x0,6*xs,sigma*xs)
	//CurveFitDialog/ variable S7N=StatsPoissonPDF(7,mu)*StatsNormalPDF(x-x0,7*xs,sigma*xs)
	//CurveFitDialog/ variable S8N=StatsPoissonPDF(8,mu)*StatsNormalPDF(x-x0,8*xs,sigma*xs)
	//CurveFitDialog/ variable S9N=StatsPoissonPDF(9,mu)*StatsNormalPDF(x-x0,9*xs,sigma*xs)
	//CurveFitDialog/ variable S10N=StatsPoissonPDF(10,mu)*StatsNormalPDF(x-x0,10*xs,sigma*xs)
	//CurveFitDialog/ variable S11N=StatsPoissonPDF(11,mu)*StatsNormalPDF(x-x0,11*xs,sigma*xs)
	//CurveFitDialog/ variable S12N=StatsPoissonPDF(12,mu)*StatsNormalPDF(x-x0,12*xs,sigma*xs)
	//CurveFitDialog/ variable S13N=StatsPoissonPDF(13,mu)*StatsNormalPDF(x-x0,13*xs,sigma*xs)
	//CurveFitDialog/ variable S14N=StatsPoissonPDF(14,mu)*StatsNormalPDF(x-x0,14*xs,sigma*xs)
	//CurveFitDialog/ variable S15N=StatsPoissonPDF(15,mu)*StatsNormalPDF(x-x0,15*xs,sigma*xs)
	//CurveFitDialog/ variable S16N=StatsPoissonPDF(16,mu)*StatsNormalPDF(x-x0,16*xs,sigma*xs)
	//CurveFitDialog/ variable S17N=StatsPoissonPDF(17,mu)*StatsNormalPDF(x-x0,17*xs,sigma*xs)
	//CurveFitDialog/ variable S18N=StatsPoissonPDF(18,mu)*StatsNormalPDF(x-x0,18*xs,sigma*xs)
	//CurveFitDialog/ variable S19N=StatsPoissonPDF(19,mu)*StatsNormalPDF(x-x0,19*xs,sigma*xs)
	//CurveFitDialog/ variable S20N=StatsPoissonPDF(20,mu)*StatsNormalPDF(x-x0,20*xs,sigma*xs)
	//CurveFitDialog/ 
	//CurveFitDialog/ 
	//CurveFitDialog/ f(x) = A*(S0N+S1N+S2N+S3N+S4N+S5N+S6N+S7N+S8N+S9N+S10N+S11N+S12N+S13N+S14N+S15N+S16N+S17N+S18N+S19N+S20N)
	//CurveFitDialog/ 
	//CurveFitDialog/ 
	//CurveFitDialog/ 
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ x
	//CurveFitDialog/ Coefficients 5
	//CurveFitDialog/ w[0] = mu
	//CurveFitDialog/ w[1] = xs
	//CurveFitDialog/ w[2] = A
	//CurveFitDialog/ w[3] = sigma
	//CurveFitDialog/ w[4] = x0

	variable S0N=StatsPoissonPDF(0,w[0])*StatsNormalPDF(x-w[4],0*w[1],w[3]*w[1])
	variable S1N=StatsPoissonPDF(1,w[0])*StatsNormalPDF(x-w[4],1*w[1],w[3]*w[1])
	variable S2N=StatsPoissonPDF(2,w[0])*StatsNormalPDF(x-w[4],2*w[1],w[3]*w[1])
	variable S3N=StatsPoissonPDF(3,w[0])*StatsNormalPDF(x-w[4],3*w[1],w[3]*w[1])
	variable S4N=StatsPoissonPDF(4,w[0])*StatsNormalPDF(x-w[4],4*w[1],w[3]*w[1])
	variable S5N=StatsPoissonPDF(5,w[0])*StatsNormalPDF(x-w[4],5*w[1],w[3]*w[1])
	variable S6N=StatsPoissonPDF(6,w[0])*StatsNormalPDF(x-w[4],6*w[1],w[3]*w[1])
	variable S7N=StatsPoissonPDF(7,w[0])*StatsNormalPDF(x-w[4],7*w[1],w[3]*w[1])
	variable S8N=StatsPoissonPDF(8,w[0])*StatsNormalPDF(x-w[4],8*w[1],w[3]*w[1])
	variable S9N=StatsPoissonPDF(9,w[0])*StatsNormalPDF(x-w[4],9*w[1],w[3]*w[1])
	variable S10N=StatsPoissonPDF(10,w[0])*StatsNormalPDF(x-w[4],10*w[1],w[3]*w[1])
	variable S11N=StatsPoissonPDF(11,w[0])*StatsNormalPDF(x-w[4],11*w[1],w[3]*w[1])
	variable S12N=StatsPoissonPDF(12,w[0])*StatsNormalPDF(x-w[4],12*w[1],w[3]*w[1])
	variable S13N=StatsPoissonPDF(13,w[0])*StatsNormalPDF(x-w[4],13*w[1],w[3]*w[1])
	variable S14N=StatsPoissonPDF(14,w[0])*StatsNormalPDF(x-w[4],14*w[1],w[3]*w[1])
	variable S15N=StatsPoissonPDF(15,w[0])*StatsNormalPDF(x-w[4],15*w[1],w[3]*w[1])
	variable S16N=StatsPoissonPDF(16,w[0])*StatsNormalPDF(x-w[4],16*w[1],w[3]*w[1])
	variable S17N=StatsPoissonPDF(17,w[0])*StatsNormalPDF(x-w[4],17*w[1],w[3]*w[1])
	variable S18N=StatsPoissonPDF(18,w[0])*StatsNormalPDF(x-w[4],18*w[1],w[3]*w[1])
	variable S19N=StatsPoissonPDF(19,w[0])*StatsNormalPDF(x-w[4],19*w[1],w[3]*w[1])
	variable S20N=StatsPoissonPDF(20,w[0])*StatsNormalPDF(x-w[4],20*w[1],w[3]*w[1])
	
	
	return w[2]*(S0N+S1N+S2N+S3N+S4N+S5N+S6N+S7N+S8N+S9N+S10N+S11N+S12N+S13N+S14N+S15N+S16N+S17N+S18N+S19N+S20N)
	
	
	
End

Function ModifiedPoisson(w,x) : FitFunc
	Wave w
	Variable x

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ variable S0N=StatsPoissonPDF(0,mu)*StatsNormalPDF(x-x0,0*xs,sigma*xs)
	//CurveFitDialog/ variable S1=StatsPoissonPDF(1,mu)*Exp((lambda/2)*(2*1*xs+lambda*(xs*sigma)^2-2*(2*1*xs-(x-x0))))*(lambda/2)*Erfc((1*xs+lambda*(xs*sigma)^2-(2*1*xs-(x-x0)))/(sqrt(2)*xs*sigma))
	//CurveFitDialog/ variable S2=StatsPoissonPDF(2,mu)*Exp((lambda/2)*(2*2*xs+lambda*(xs*sigma)^2-2*(2*2*xs-(x-x0))))*(lambda/2)*Erfc((2*xs+lambda*(xs*sigma)^2-(2*2*xs-(x-x0)))/(sqrt(2)*xs*sigma))
	//CurveFitDialog/ variable S3=StatsPoissonPDF(3,mu)*Exp((lambda/2)*(2*3*xs+lambda*(xs*sigma)^2-2*(2*3*xs-(x-x0))))*(lambda/2)*Erfc((3*xs+lambda*(xs*sigma)^2-(2*3*xs-(x-x0)))/(sqrt(2)*xs*sigma))
	//CurveFitDialog/ variable S4=StatsPoissonPDF(4,mu)*Exp((lambda/2)*(2*4*xs+lambda*(xs*sigma)^2-2*(2*4*xs-(x-x0))))*(lambda/2)*Erfc((4*xs+lambda*(xs*sigma)^2-(2*4*xs-(x-x0)))/(sqrt(2)*xs*sigma))
	//CurveFitDialog/ variable S5=StatsPoissonPDF(5,mu)*Exp((lambda/2)*(2*5*xs+lambda*(xs*sigma)^2-2*(2*5*xs-(x-x0))))*(lambda/2)*Erfc((5*xs+lambda*(xs*sigma)^2-(2*5*xs-(x-x0)))/(sqrt(2)*xs*sigma))
	//CurveFitDialog/ 
	//CurveFitDialog/ 
	//CurveFitDialog/ 
	//CurveFitDialog/ 
	//CurveFitDialog/ 
	//CurveFitDialog/ f(x) = A*(S0N+S1+S2+S3+S4+S5)
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ x
	//CurveFitDialog/ Coefficients 6
	//CurveFitDialog/ w[0] = A
	//CurveFitDialog/ w[1] = mu
	//CurveFitDialog/ w[2] = xs
	//CurveFitDialog/ w[3] = x0
	//CurveFitDialog/ w[4] = sigma
	//CurveFitDialog/ w[5] = lambda

	variable S0N=StatsPoissonPDF(0,w[1])*StatsNormalPDF(x-w[3],0*w[2],w[4]*w[2])
	variable S1=StatsPoissonPDF(1,w[1])*Exp((w[5]/2)*(2*1*w[2]+w[5]*(w[2]*w[4])^2-2*(2*1*w[2]-(x-w[3]))))*(w[5]/2)*Erfc((1*w[2]+w[5]*(w[2]*w[4])^2-(2*1*w[2]-(x-w[3])))/(sqrt(2)*w[2]*w[4]))
	variable S2=StatsPoissonPDF(2,w[1])*Exp((w[5]/2)*(2*2*w[2]+w[5]*(w[2]*w[4])^2-2*(2*2*w[2]-(x-w[3]))))*(w[5]/2)*Erfc((2*w[2]+w[5]*(w[2]*w[4])^2-(2*2*w[2]-(x-w[3])))/(sqrt(2)*w[2]*w[4]))
	variable S3=StatsPoissonPDF(3,w[1])*Exp((w[5]/2)*(2*3*w[2]+w[5]*(w[2]*w[4])^2-2*(2*3*w[2]-(x-w[3]))))*(w[5]/2)*Erfc((3*w[2]+w[5]*(w[2]*w[4])^2-(2*3*w[2]-(x-w[3])))/(sqrt(2)*w[2]*w[4]))
	variable S4=StatsPoissonPDF(4,w[1])*Exp((w[5]/2)*(2*4*w[2]+w[5]*(w[2]*w[4])^2-2*(2*4*w[2]-(x-w[3]))))*(w[5]/2)*Erfc((4*w[2]+w[5]*(w[2]*w[4])^2-(2*4*w[2]-(x-w[3])))/(sqrt(2)*w[2]*w[4]))
	variable S5=StatsPoissonPDF(5,w[1])*Exp((w[5]/2)*(2*5*w[2]+w[5]*(w[2]*w[4])^2-2*(2*5*w[2]-(x-w[3]))))*(w[5]/2)*Erfc((5*w[2]+w[5]*(w[2]*w[4])^2-(2*5*w[2]-(x-w[3])))/(sqrt(2)*w[2]*w[4]))
	
	
	
	
	
	return w[0]*(S0N+S1+S2+S3+S4+S5)
End

Function Poisson_x_scaled_b(w,x) : FitFunc
	Wave w
	Variable x

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ variable S0N=StatsPoissonPDF(0,mu)*StatsNormalPDF(x-x0,0*xs,sigma*xs)
	//CurveFitDialog/ variable S1N=StatsPoissonPDF(1,mu)*StatsNormalPDF(x-x0,1*xs,b*sigma*xs)
	//CurveFitDialog/ variable S2N=StatsPoissonPDF(2,mu)*StatsNormalPDF(x-x0,2*xs,b*sigma*xs)
	//CurveFitDialog/ variable S3N=StatsPoissonPDF(3,mu)*StatsNormalPDF(x-x0,3*xs,b*sigma*xs)
	//CurveFitDialog/ variable S4N=StatsPoissonPDF(4,mu)*StatsNormalPDF(x-x0,4*xs,b*sigma*xs)
	//CurveFitDialog/ variable S5N=StatsPoissonPDF(5,mu)*StatsNormalPDF(x-x0,5*xs,b*sigma*xs)
	//CurveFitDialog/ variable S6N=StatsPoissonPDF(6,mu)*StatsNormalPDF(x-x0,6*xs,sigma*xs)
	//CurveFitDialog/ variable S7N=StatsPoissonPDF(7,mu)*StatsNormalPDF(x-x0,7*xs,sigma*xs)
	//CurveFitDialog/ variable S8N=StatsPoissonPDF(8,mu)*StatsNormalPDF(x-x0,8*xs,sigma*xs)
	//CurveFitDialog/ variable S9N=StatsPoissonPDF(9,mu)*StatsNormalPDF(x-x0,9*xs,sigma*xs)
	//CurveFitDialog/ variable S10N=StatsPoissonPDF(10,mu)*StatsNormalPDF(x-x0,10*xs,sigma*xs)
	//CurveFitDialog/ variable S11N=StatsPoissonPDF(11,mu)*StatsNormalPDF(x-x0,11*xs,sigma*xs)
	//CurveFitDialog/ variable S12N=StatsPoissonPDF(12,mu)*StatsNormalPDF(x-x0,12*xs,sigma*xs)
	//CurveFitDialog/ variable S13N=StatsPoissonPDF(13,mu)*StatsNormalPDF(x-x0,13*xs,sigma*xs)
	//CurveFitDialog/ variable S14N=StatsPoissonPDF(14,mu)*StatsNormalPDF(x-x0,14*xs,sigma*xs)
	//CurveFitDialog/ variable S15N=StatsPoissonPDF(15,mu)*StatsNormalPDF(x-x0,15*xs,sigma*xs)
	//CurveFitDialog/ variable S16N=StatsPoissonPDF(16,mu)*StatsNormalPDF(x-x0,16*xs,sigma*xs)
	//CurveFitDialog/ variable S17N=StatsPoissonPDF(17,mu)*StatsNormalPDF(x-x0,17*xs,sigma*xs)
	//CurveFitDialog/ variable S18N=StatsPoissonPDF(18,mu)*StatsNormalPDF(x-x0,18*xs,sigma*xs)
	//CurveFitDialog/ variable S19N=StatsPoissonPDF(19,mu)*StatsNormalPDF(x-x0,19*xs,sigma*xs)
	//CurveFitDialog/ variable S20N=StatsPoissonPDF(20,mu)*StatsNormalPDF(x-x0,20*xs,sigma*xs)
	//CurveFitDialog/ 
	//CurveFitDialog/ 
	//CurveFitDialog/ f(x) = A*(S0N+S1N+S2N+S3N+S4N+S5N+S6N+S7N+S8N+S9N+S10N+S11N+S12N+S13N+S14N+S15N+S16N+S17N+S18N+S19N+S20N)
	//CurveFitDialog/ 
	//CurveFitDialog/ 
	//CurveFitDialog/ 
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ x
	//CurveFitDialog/ Coefficients 6
	//CurveFitDialog/ w[0] = mu
	//CurveFitDialog/ w[1] = xs
	//CurveFitDialog/ w[2] = A
	//CurveFitDialog/ w[3] = sigma
	//CurveFitDialog/ w[4] = x0
	//CurveFitDialog/ w[5] = b

	variable S0N=StatsPoissonPDF(0,w[0])*StatsNormalPDF(x-w[4],0*w[1],w[3]*w[1])
	variable S1N=StatsPoissonPDF(1,w[0])*StatsNormalPDF(x-w[4],1*w[1],w[5]*w[3]*w[1])
	variable S2N=StatsPoissonPDF(2,w[0])*StatsNormalPDF(x-w[4],2*w[1],w[5]*w[3]*w[1])
	variable S3N=StatsPoissonPDF(3,w[0])*StatsNormalPDF(x-w[4],3*w[1],w[5]*w[3]*w[1])
	variable S4N=StatsPoissonPDF(4,w[0])*StatsNormalPDF(x-w[4],4*w[1],w[5]*w[3]*w[1])
	variable S5N=StatsPoissonPDF(5,w[0])*StatsNormalPDF(x-w[4],5*w[1],w[5]*w[3]*w[1])
	variable S6N=StatsPoissonPDF(6,w[0])*StatsNormalPDF(x-w[4],6*w[1],w[3]*w[1])
	variable S7N=StatsPoissonPDF(7,w[0])*StatsNormalPDF(x-w[4],7*w[1],w[3]*w[1])
	variable S8N=StatsPoissonPDF(8,w[0])*StatsNormalPDF(x-w[4],8*w[1],w[3]*w[1])
	variable S9N=StatsPoissonPDF(9,w[0])*StatsNormalPDF(x-w[4],9*w[1],w[3]*w[1])
	variable S10N=StatsPoissonPDF(10,w[0])*StatsNormalPDF(x-w[4],10*w[1],w[3]*w[1])
	variable S11N=StatsPoissonPDF(11,w[0])*StatsNormalPDF(x-w[4],11*w[1],w[3]*w[1])
	variable S12N=StatsPoissonPDF(12,w[0])*StatsNormalPDF(x-w[4],12*w[1],w[3]*w[1])
	variable S13N=StatsPoissonPDF(13,w[0])*StatsNormalPDF(x-w[4],13*w[1],w[3]*w[1])
	variable S14N=StatsPoissonPDF(14,w[0])*StatsNormalPDF(x-w[4],14*w[1],w[3]*w[1])
	variable S15N=StatsPoissonPDF(15,w[0])*StatsNormalPDF(x-w[4],15*w[1],w[3]*w[1])
	variable S16N=StatsPoissonPDF(16,w[0])*StatsNormalPDF(x-w[4],16*w[1],w[3]*w[1])
	variable S17N=StatsPoissonPDF(17,w[0])*StatsNormalPDF(x-w[4],17*w[1],w[3]*w[1])
	variable S18N=StatsPoissonPDF(18,w[0])*StatsNormalPDF(x-w[4],18*w[1],w[3]*w[1])
	variable S19N=StatsPoissonPDF(19,w[0])*StatsNormalPDF(x-w[4],19*w[1],w[3]*w[1])
	variable S20N=StatsPoissonPDF(20,w[0])*StatsNormalPDF(x-w[4],20*w[1],w[3]*w[1])
	
	
	return w[2]*(S0N+S1N+S2N+S3N+S4N+S5N+S6N+S7N+S8N+S9N+S10N+S11N+S12N+S13N+S14N+S15N+S16N+S17N+S18N+S19N+S20N)
	
	
	
End

Function NModifiedPoisson(w,x) : FitFunc
	Wave w
	Variable x

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ variable S0N=StatsPoissonPDF(0,mu)*Exp((lambda/2)*(2*0*xs+lambda*(xs*sigma)^2-2*((x-x0))))*(lambda/2)*Erfc((0*xs+lambda*(xs*sigma)^2-((x-x0)))/(sqrt(2)*xs*sigma))
	//CurveFitDialog/ variable S1=StatsPoissonPDF(1,mu)*Exp((lambda/2)*(2*1*xs+lambda*(xs*sigma)^2-2*((x-x0))))*(lambda/2)*Erfc((1*xs+lambda*(xs*sigma)^2-((x-x0)))/(sqrt(2)*xs*sigma))
	//CurveFitDialog/ variable S2=StatsPoissonPDF(2,mu)*Exp((lambda/2)*(2*2*xs+lambda*(xs*sigma)^2-2*((x-x0))))*(lambda/2)*Erfc((2*xs+lambda*(xs*sigma)^2-((x-x0)))/(sqrt(2)*xs*sigma))
	//CurveFitDialog/ variable S3=StatsPoissonPDF(3,mu)*Exp((lambda/2)*(2*3*xs+lambda*(xs*sigma)^2-2*((x-x0))))*(lambda/2)*Erfc((3*xs+lambda*(xs*sigma)^2-((x-x0)))/(sqrt(2)*xs*sigma))
	//CurveFitDialog/ variable S4=StatsPoissonPDF(4,mu)*Exp((lambda/2)*(2*4*xs+lambda*(xs*sigma)^2-2*((x-x0))))*(lambda/2)*Erfc((4*xs+lambda*(xs*sigma)^2-((x-x0)))/(sqrt(2)*xs*sigma))
	//CurveFitDialog/ variable S5=StatsPoissonPDF(5,mu)*Exp((lambda/2)*(2*5*xs+lambda*(xs*sigma)^2-2*((x-x0))))*(lambda/2)*Erfc((5*xs+lambda*(xs*sigma)^2-((x-x0)))/(sqrt(2)*xs*sigma))
	//CurveFitDialog/ 
	//CurveFitDialog/ 
	//CurveFitDialog/ 
	//CurveFitDialog/ 
	//CurveFitDialog/ 
	//CurveFitDialog/ f(x) = A*(S0N+S1+S2+S3+S4+S5)
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ x
	//CurveFitDialog/ Coefficients 6
	//CurveFitDialog/ w[0] = A
	//CurveFitDialog/ w[1] = mu
	//CurveFitDialog/ w[2] = xs
	//CurveFitDialog/ w[3] = x0
	//CurveFitDialog/ w[4] = sigma
	//CurveFitDialog/ w[5] = lambda

	variable S0N=StatsPoissonPDF(0,w[1])*Exp((w[5]/2)*(2*0*w[2]+w[5]*(w[2]*w[4])^2-2*((x-w[3]))))*(w[5]/2)*Erfc((0*w[2]+w[5]*(w[2]*w[4])^2-((x-w[3])))/(sqrt(2)*w[2]*w[4]))
	variable S1=StatsPoissonPDF(1,w[1])*Exp((w[5]/2)*(2*1*w[2]+w[5]*(w[2]*w[4])^2-2*((x-w[3]))))*(w[5]/2)*Erfc((1*w[2]+w[5]*(w[2]*w[4])^2-((x-w[3])))/(sqrt(2)*w[2]*w[4]))
	variable S2=StatsPoissonPDF(2,w[1])*Exp((w[5]/2)*(2*2*w[2]+w[5]*(w[2]*w[4])^2-2*((x-w[3]))))*(w[5]/2)*Erfc((2*w[2]+w[5]*(w[2]*w[4])^2-((x-w[3])))/(sqrt(2)*w[2]*w[4]))
	variable S3=StatsPoissonPDF(3,w[1])*Exp((w[5]/2)*(2*3*w[2]+w[5]*(w[2]*w[4])^2-2*((x-w[3]))))*(w[5]/2)*Erfc((3*w[2]+w[5]*(w[2]*w[4])^2-((x-w[3])))/(sqrt(2)*w[2]*w[4]))
	variable S4=StatsPoissonPDF(4,w[1])*Exp((w[5]/2)*(2*4*w[2]+w[5]*(w[2]*w[4])^2-2*((x-w[3]))))*(w[5]/2)*Erfc((4*w[2]+w[5]*(w[2]*w[4])^2-((x-w[3])))/(sqrt(2)*w[2]*w[4]))
	variable S5=StatsPoissonPDF(5,w[1])*Exp((w[5]/2)*(2*5*w[2]+w[5]*(w[2]*w[4])^2-2*((x-w[3]))))*(w[5]/2)*Erfc((5*w[2]+w[5]*(w[2]*w[4])^2-((x-w[3])))/(sqrt(2)*w[2]*w[4]))
	
	
	
	
	
	return w[0]*(S0N+S1+S2+S3+S4+S5)
End