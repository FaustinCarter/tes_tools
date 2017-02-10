#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#include <Multi-peak fitting 2.0>


function showFits(Coefs,scaleWave,type,ap)
//type 1 = Fast, type 2 = Slow
//set ap to zero for a new graph

wave Coefs, scaleWave
variable type, ap
variable i, size
string name

size=DimSize(Coefs, 1)

Make/O/N=(DimSize(Coefs,0)) tmpCoef
Make/O/N=2000 tmpFit
CopyScales scaleWave, tmpFit

if(type==1)
	name="fitWavesFast"
elseif(type==2)
	name="fitWavesSlow"
endif

wave w=$name

if(waveexists($name))
	KillWaves $name
endif

if(ap==0)
	Display
endif

for (i=0;i<size;i+=1)
	tmpCoef=Coefs[p][i]
	
	if(type==1)
		tmpFit=EMGFunction(tmpCoef,x)
	elseif(type==2)
		tmpFit=ExpConvExpFunc(tmpCoef,x)
	endif
	
	Concatenate {tmpFit}, $name
	AppendToGraph $name[][i]
endfor

KillWaves tmpCoef, tmpFit

end

function IntegrateFits(Coefs,scaleWave,type)
//type 1 = fast, ExpModGauss, type 2 = slow, ExpConvExp
wave Coefs, scaleWave
variable type
variable i, size
string name

size=DimSize(Coefs, 1)

Make/O/N=(DimSize(Coefs,0)) tmpCoef
Make/O/N=2000 tmpFit
CopyScales scaleWave, tmpFit

if(type==1)
	name="fitWavesFast"
elseif(type==2)
	name="fitWavesSlow"
endif

Make/O/N=(size) $name+"_Ints"
wave w=$name+"_Ints"


for (i=0;i<size;i+=1)
	tmpCoef=Coefs[p][i]
	if(type==1)
		tmpFit=EMGFunction(tmpCoef,x)
	elseif(type==2)
		tmpFit=ExpConvExpFunc(tmpCoef,x)
	endif
	
	Concatenate {tmpFit}, $name

	w[i]=area(tmpFit)
	
	if(mod(i,1000)==0)
		print i
	endif
endfor

KillWaves tmpCoef, tmpFit

end

function IntegrateFits2(Coefs,scaleWave,type)
//type 1 = fast, ExpModExp, type 2 = slow, ExpConvExp
wave Coefs, scaleWave
variable type
variable i, size
string name

size=DimSize(Coefs, 1)

Make/O/N=(DimSize(Coefs,0)) tmpCoef
Make/O/N=2000 tmpFit
CopyScales scaleWave, tmpFit

if(type==1)
	name="fitWavesFast"
elseif(type==2)
	name="fitWavesSlow"
endif

Make/O/N=(size) $name+"_Ints"
wave w=$name+"_Ints"


for (i=0;i<size;i+=1)
	tmpCoef=Coefs[p][i]

		tmpFit=ExpConvExpFunc(tmpCoef,x)
	
	Concatenate {tmpFit}, $name

	w[i]=area(tmpFit)
	
	if(mod(i,1000)==0)
		print i
	endif
endfor

KillWaves tmpCoef, tmpFit

end
