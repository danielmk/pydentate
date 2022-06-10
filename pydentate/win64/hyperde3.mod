TITLE hyperde3.mod  
 
COMMENT
Chen K, Aradi I, Thon N, Eghbal-Ahmadi M, Baram TZ, Soltesz I: Persistently
modified
h-channels after complex febrile seizures convert the seizure-induced
enhancement of
inhibition to hyperexcitability. Nature Medicine, 7(3) pp. 331-337, 2001.
(modeling by Ildiko Aradi, iaradi@uci.edu)
distal dendritic Ih channel kinetics for both HT and Control anlimals
ENDCOMMENT
 
UNITS {
        (mA) =(milliamp)
        (mV) =(millivolt)
        (uF) = (microfarad)
	(molar) = (1/liter)
	(nA) = (nanoamp)
	(mM) = (millimolar)
	(um) = (micron)
	FARADAY = 96520 (coul)
	R = 8.3134	(joule/degC)
}
 
? interface 
NEURON { 
SUFFIX hyperde3 
USEION hyf READ ehyf WRITE ihyf VALENCE 1
USEION hys READ ehys WRITE ihys VALENCE 1
USEION hyhtf READ ehyhtf WRITE ihyhtf VALENCE 1
USEION hyhts READ ehyhts WRITE ihyhts VALENCE 1
RANGE  ghyf, ghys, ghyhtf, ghyhts
RANGE ghyfbar, ghysbar, ghyhtfbar, ghyhtsbar
RANGE hyfinf, hysinf, hyftau, hystau
RANGE hyhtfinf, hyhtsinf, hyhtftau, hyhtstau, ihyf, ihys
}
 
INDEPENDENT {t FROM 0 TO 100 WITH 100 (ms)}
 
PARAMETER {
      v (mV) 
      celsius = 6.3 (degC)
      dt (ms) 

	ghyfbar (mho/cm2)
	ghysbar (mho/cm2)
	ehyf (mV)
	ehys (mV)
	ghyhtfbar (mho/cm2)
	ghyhtsbar (mho/cm2)
	ehyhtf (mV)
	ehyhts (mV)
}
 
STATE {
	hyf hys hyhtf hyhts
}
 
ASSIGNED {
         
  
	ghyf (mho/cm2)
 	ghys (mho/cm2)

	ghyhtf (mho/cm2)
	ghyhts (mho/cm2)

  
	ihyf (mA/cm2)
	ihys (mA/cm2)
	ihyhtf (mA/cm2)
	ihyhts (mA/cm2)

	hyfinf hysinf hyhtfinf hyhtsinf
 	hyftau (ms) hystau (ms) hyhtftau (ms) hyhtstau (ms)
	hyfexp hysexp hyhtfexp hyhtsexp     
} 

? currents
BREAKPOINT {

	SOLVE states

	ghyf = ghyfbar * hyf*hyf
	ihyf = ghyf * (v-ehyf)
	ghys = ghysbar * hys*hys
	ihys = ghys * (v-ehys)

	ghyhtf = ghyhtfbar * hyhtf* hyhtf
	ihyhtf = ghyhtf * (v-ehyhtf)
	ghyhts = ghyhtsbar * hyhts* hyhts
	ihyhts = ghyhts * (v-ehyhts)
		
		}
 
UNITSOFF
 
INITIAL {
	trates(v)
	
	hyf = hyfinf
      hys = hysinf
	hyhtf = hyhtfinf
	hyhts = hyhtsinf
	VERBATIM
	return 0;
	ENDVERBATIM
}

? states
PROCEDURE states() {	:Computes state variables m, h, and n 
        trates(v)	:      at the current v and dt.
        
        hyf = hyf + hyfexp*(hyfinf-hyf)
        hys = hys + hysexp*(hysinf-hys)
	  hyhtf = hyhtf + hyhtfexp*(hyhtfinf-hyhtf)
	  hyhts = hyhts + hyhtsexp*(hyhtsinf-hyhts)

        VERBATIM
        return 0;
        ENDVERBATIM
}
 
LOCAL q10

? rates
PROCEDURE rates(v) {  :Computes rate and other constants at current v.
                      :Call once from HOC to initialize inf at resting v.
        LOCAL  alpha, beta, sum
       q10 = 3^((celsius - 6.3)/10)
       
	:"hyf" FAST CONTROL Hype activation system
	hyfinf =  1 / (1 + exp( (v+91)/10 ))
	hyftau = 14.9 + 14.1 / (1+exp(-(v+95.2)/0.5))

	:"hys" SLOW CONTROL Hype activation system
	hysinf =  1 / (1 + exp( (v+91)/10 ))
	hystau = 80 + 172.7 / (1+exp(-(v+59.3)/-0.83))

		:"hyhtf" FAST HT Hypeht activation system
	hyhtfinf =  1 / (1 + exp( (v+87)/10 ))
	hyhtftau = 23.2 + 16.1 / (1+exp(-(v+91.2)/0.83))

		:"hyhts" SLOW HT Hypeht activation system
	hyhtsinf =  1 / (1 + exp( (v+87)/10 ))
	hyhtstau = 227.3 + 170.7*exp(-0.5*((v+80.4)/11)^2)
}
 
PROCEDURE trates(v) {  :Computes rate and other constants at current v.
                      :Call once from HOC to initialize inf at resting v.
	LOCAL tinc
      TABLE hyfinf, hyhtfinf, hyfexp, hyhtfexp, hyftau, hyhtftau, 
		hysinf, hyhtsinf, hysexp, hyhtsexp, hystau, hyhtstau	
	DEPEND dt, celsius FROM -120 TO 100 WITH 220
                           
	rates(v)	: not consistently executed from here if usetable_hh == 1
		: so don't expect the tau values to be tracking along with
		: the inf values in hoc

	       tinc = -dt * q10
        
        hyfexp = 1 - exp(tinc/hyftau)
	  hysexp = 1 - exp(tinc/hystau)
	  hyhtfexp = 1 - exp(tinc/hyhtftau)
	  hyhtsexp = 1 - exp(tinc/hyhtstau)
}
 
FUNCTION vtrap(x,y) {  :Traps for 0 in denominator of rate eqns.
        if (fabs(x/y) < 1e-6) {
                vtrap = y*(1 - x/y/2)
        }else{  
                vtrap = x/(exp(x/y) - 1)
        }
}
 
UNITSON

