#pragma rtGlobals=1		// Use modern global access method.
function Vsense2Rtes(Rb, Vb, Rs, Rh, wRt,[gain, titletxt, Rl])
//If the RvT data includes the offset from the leads, don't specify Rl. If the RvT data goes to zero at T<Tc, specify Rl.

Variable Rb, Vb, Rs, Rh, Rl, gain
String titletxt
Wave wRt

if(ParamIsDefault(gain))
	gain=1
endif

if(ParamIsDefault(titletxt))
	titletxt="To set the title, call again\rwith the data folder name.\rEx: titletxt=\"page97\""
else
	titletxt="Data from experiment " + titletxt
endif

if(ParamIsDefault(Rl))
	Rl=0
endif

Variable Ib, V_npnts
WaveStats/Q wRt
Make/O/N=(V_npnts) wVs
Make/O/N=(V_npnts) wIs
CopyScales wRt wVs
CopyScales wRt wIs


Ib=Vb/Rb //bias current
wVs=1e6*Ib*Rh*Rs/(Rh+Rs+Rl+wRt(x)) //wave containing the voltage accross the sense resistor vs tes temperature in uV
wIs=1e6*Ib*Rh/(Rh+Rs+Rl+wRt(x)) //wave containing the current through the sense resistor vs tes temperature in uA



Display wVs vs Rtes //plot the voltage accross the sense resistor vs temperature
ModifyGraph tick=2,mirror=1
Label left "Voltage accross R\Bsense\M (uV)"
Label bottom "TES + lead resistance (Ohms)"
SetAxis/A/E=1 bottom
TextBox/C/N=descipt0/A=MC titletxt+"\r\rThis is the range of possible voltage values\rthat the scope or DMM could see over the\rfull range of the TES.\r\rMax Delta="+num2str(wVs[0]-wVs[V_npnts-1])+" uV"



end //end function
