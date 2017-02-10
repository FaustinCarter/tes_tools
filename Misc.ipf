#pragma rtGlobals=3		// Use modern global access method and strict wave access.

function analyzerTemplate(chan)
	variable chan
	wave allDataF, allDataF_areas, allDataF_widths
	pulseHeightsInt(allDataF, 0.3, chan)
	Duplicate/O allDataF_areas allDataF_areasM
	Duplicate/O allDataF_areas, mask
	mask=0
	mask=(allDataF_widths>700e-9)&&(allDataF_widths<900e-9)?1:0
	allDataF_areasM=mask==1?allDataF_areas[p]:nan
	allDataF_areasM/=0.34
	
	Make/N=30/O allDataF_areasM_Hist;DelayUpdate
	Histogram/B=1 allDataF_areasM,allDataF_areasM_Hist
end

function runAnalyzerTemplate() //Macro, edit as needed
	string folders="t1700;t1800;t1850;t1900;t1950;t2000"
	variable i
	string fName
	
	//Display
	
	for(i=0;i<ItemsInList(folders);i+=1)
		fName="root:"+StringFromList(i, folders)+":filtered"
		SetDataFolder $fName
		
		analyzerTemplate(2)
		
		//wave allDataF_areasM_Hist
		//AppendToGraph allDataF_areasM_Hist
		
	endfor
	
	SetDataFolder root:
end