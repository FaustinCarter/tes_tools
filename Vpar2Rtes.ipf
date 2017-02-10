#pragma rtGlobals=1		// Use modern global access method.
function Vpar2Rtes(Rb, Vb, Rh, wRt,[gain, titletxt, Rl])
//If the RvT data includes the offset from the leads, don't specify Rl. If the RvT data goes to zero at T<Tc, specify Rl.

Variable Rb, Vb, Rh, Rl, gain
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
Make/O/N=(V_npnts) wVh
CopyScales wRt wVh


Ib=Vb/Rb //bias current
wVh=1e6*Ib*Rh*(Rl+wRt(x))/(Rh+Rl+wRt(x)) //wave containing the voltage accross the sense resistor vs tes temperature in uV

Display wVh vs wRt //plot the voltage accross the sense resistor vs temperature
ModifyGraph tick=2,mirror=1
Label left "Voltage accross R\Bshunt\M (uV)"
Label bottom "TES + lead resistance (Ohms)"
SetAxis/A/E=1 bottom
TextBox/C/N=descipt0/A=MC titletxt+"\r\rThis is the range of possible voltage values\rthat the scope or DMM could see over the\rfull range of the TES.\r\rMax Delta="+num2str(wVh[0]-wVh[V_npnts-1])+" uV"



end //end function
