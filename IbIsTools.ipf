#pragma rtGlobals=3		// Use modern global access method and strict wave access.


Function ScaleStuff(Suffix)

	string Suffix
	string Rname, Pname, Vname
	
	Pname="Pwr"+Suffix
	Rname="Rt"+Suffix
	Vname="Vt"+Suffix

	wave wS = $("Is"+Suffix)
	wave wB = $("Ib"+Suffix)
	
	Duplicate wS wR
	Duplicate wS wP
	Duplicate wS, wV
	
	Rename wR, $Rname
	Rename wP, $Pname
	Rename wV, $Vname

	
	SetFormula $Rname,  "(Ib"+Suffix+"/Is"+Suffix+"-1)*0.2-0.00372"

	SetFormula $Pname, "Rt"+Suffix+"*Is"+Suffix+"^2"
	
	SetFormula $Vname, "Is"+Suffix+"*Rt"+Suffix
	
	AppendToGraph/W=Graph0 wS vs wB
	AppendToGraph/W=Graph1 wR vs wB
	AppendToGraph/W=Graph2 wP vs wB
	AppendToGraph/W=Graph3 wV vs wB
	
End

Function makeCurves(temps)

	Wave temps
	Variable i=0
	
	Display/N=Graph0
	Display/N=Graph1
	Display/N=Graph2
	Display/N=Graph3
	
	For(i=0; i<DimSize(temps, 0);i+=1)
		ScaleStuff(num2str(temps[i]))
	EndFor
End


Function OffsetAll(graphName, offset, slope)
	
	String graphName		// Must contain the name of a graph
	Variable offset, slope
	
	
	Variable index = 0
	
	do
		// Get wave reference for next Y wave in graph
		WAVE/Z wy = WaveRefIndexed(graphName,index,1)
		WAVE/Z wx = WaveRefIndexed(graphName,index,2)
		if (WaveExists(wy) == 0)
			break							// No more waves
		endif
		
		wy=wy+offset
		wx=wx+offset/slope	
		index += 1
	while(1)								// Loop till break above

End

Function OffsetOne(suffix, slope, offset)
	string suffix
	variable slope, offset
	
	string wName
	
	wName= "Is"+suffix
	Wave wy=$wName
	
	wName="Ib"+suffix
	Wave wx=$wName
	
	wy=wy+offset
	wx=wx+offset/slope
	
End


Function FindPowerAt(graphName, current)
	
	String graphName		// Must contain the name of a graph
	Variable current
	
	
	Variable index = 0
	Variable startPoint
	
	Make/O $("Pwrs"+num2istr(current))
	
	wave w=$("Pwrs"+num2istr(current))
	
	do
		// Get wave reference for next Y wave in graph
		WAVE/Z wy = WaveRefIndexed(graphName,index,1)
		WAVE/Z wx = WaveRefIndexed(graphName,index,2)
		if (WaveExists(wy) == 0)
			break							// No more waves
		endif
		
		startPoint=BinarySearchInterp(wx,current)
		WaveStats/Q/R=[(startPoint-5), (startPoint+5)] wy
		
		w[index]=V_avg
			
		index += 1
	while(1)								// Loop till break above
	
	Redimension/N=(index) w

End

Function FindResAt(graphName, current)
	
	String graphName		// Must contain the name of a graph
	Variable current
	
	
	Variable index = 0
	
	do
		WAVE/Z wy = WaveRefIndexed(graphName,index,1)
		if (WaveExists(wy) == 0)
			break							// No more waves
		endif
		index+=1
	while(1)

	index=0
	
	Make/O/N=(index-1) $("Rs"+num2istr(current))
	
	wave w=$("Rs"+num2istr(current))
	
	do
		// Get wave reference for next Y wave in graph
		WAVE/Z wy = WaveRefIndexed(graphName,index,1)
		WAVE/Z wx = WaveRefIndexed(graphName,index,2)
		if (WaveExists(wy) == 0)
			break							// No more waves
		endif
		
		w[index]=wy[BinarySearchInterp(wx,current)]
			
		index += 1
	while(1)								// Loop till break above

End

Function CalibrateIbIs(graphName, current)
	
	String graphName		// Must contain the name of a graph
	Variable current
	
	
	Variable index = 0
	
	Variable a1,b1,a2,b2,V_minloc, V_maxloc, V_startRow, V_LevelX, avgSlope
	wave W_coef

	avgSlope=0
	
	do
		// Get wave reference for next Y wave in graph
		WAVE/Z wy = WaveRefIndexed(graphName,index,1)
		WAVE/Z wx = WaveRefIndexed(graphName,index,2)
		if (WaveExists(wy) == 0)
			break							// No more waves
		endif
		
		WaveStats/Q wy
		
		CurveFit line wy[V_minloc-10,V_maxloc+10] /X=wx
		a1=W_coef[0]
		b1=W_coef[1]
		
		avgSlope+=b1
		
		FindLevel wx, current
		
		print V_LevelX
		
		CurveFit line wy[0,V_LevelX] /X=wx
		a2=W_coef[0]
		b2=W_coef[1]
		
		wy-=(-b1*(a1-a2)/(b1-b2)+a1)
		wx-=(-(a1-a2)/(b1-b2))
			
		index += 1
	while(1)								// Loop till break above

	avgSlope=avgSlope/(index)
	print avgSlope

End	

Function PowerDissipation(w,T) : FitFunc
	Wave w
	Variable T

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(T) = K*(Tc^n-T^n)
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ T
	//CurveFitDialog/ Coefficients 3
	//CurveFitDialog/ w[0] = K
	//CurveFitDialog/ w[1] = Tc
	//CurveFitDialog/ w[2] = n

	return w[0]*(w[1]^w[2]-T^w[2])
End

