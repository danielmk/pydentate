COMMENT
Revised 12/15/2000 in light of a personal communication 
from Misha Tsodyks that u is incremented _before_ x is 
converted to y--a point that was not clear in the paper.
If u is incremented _after_ x is converted to y, then 
the first synaptic activation after a long interval of 
silence will produce smaller and smaller postsynaptic 
effect as the length of the silent interval increases, 
eventually becoming vanishingly small.

Implementation of a model of short-term facilitation and depression 
based on the kinetics described in
  Tsodyks et al.
  Synchrony generation in recurrent networks 
  with frequency-dependent synapses
  Journal of Neuroscience 20:RC50:1-5, 2000.
Their mechanism represented synapses as current sources.
The mechanism implemented here uses a conductance change instead.

The basic scheme is

x -------> y    Instantaneous, spike triggered.
                Increment is u*x (see discussion of u below).
                x == fraction of "synaptic resources" that have 
                     "recovered" (fraction of xmtr pool that is 
                     ready for release, or fraction of postsynaptic 
                     channels that are ready to be opened, or some 
                     joint function of these two factors)
                y == fraction of "synaptic resources" that are in the 
                     "active state."  This is proportional to the 
                     number of channels that are open, or the 
                     fraction of max synaptic current that is 
                     being delivered. 
  tau_1
y -------> z    z == fraction of "synaptic resources" that are 
                     in the "inactive state"

  tau_rec
z -------> x

where x + y + z = 1

The active state y is multiplied by a synaptic weight to compute
the actual synaptic conductance (or current, in the original form 
of the model).

In addition, there is a "facilition" term u that 
governs the fraction of x that is converted to y 
on each synaptic activation.

  -------> u    Instantaneous, spike triggered, 
                happens _BEFORE_ x is converted to y.
                Increment is U*(1-u) where U and u both 
                lie in the range 0 - 1.
  tau_facil
u ------->      decay of facilitation

This implementation for NEURON offers the user a parameter 
u0 that has a default value of 0 but can be used to specify 
a nonzero initial value for u.

When tau_facil = 0, u is supposed to equal U.

Note that the synaptic conductance in this mechanism 
has the same kinetics as y, i.e. decays with time 
constant tau_1.

This mechanism can receive multiple streams of 
synaptic input via NetCon objects.  
Each stream keeps track of its own 
weight and activation history.

The printf() statements are for testing purposes only.
ENDCOMMENT

: THE NEURON BLOCK DEFINES HOW THE MECHANISM LOOKS LIKE FROM THE OUTSIDE
: IT TELLS US WHAT PARTS OF THE MECHANISM ARE EXPOSED TO OTHER PARTS OF THE NEURON LANGUAGE
NEURON {
	POINT_PROCESS tmgexp2syn
	: THE RANGE STATEMENT ASSERTS THAT THESE VARIABLES CAN BE ACCESSED IN HOC THROUGH RANGE VARIABLE SYNTAX
	RANGE e, i, g
	RANGE tau_1, tau_2, tau_rec, tau_facil, U, u0
	: THE NONSPECIFIC CURRENT DEFINITION HAS TWO CONSEQUENCES
	: FIRST, i HAS TO BE ACCOUNTED FOR IN CHARGE BALANCE EQUATIONS
	: SECOND, i DOES NOT CONTRIBUTE DIRECTLY TO MASS BALANCS OF ANY IONS
	NONSPECIFIC_CURRENT i
	: RANGE g
}
: UNITS BLOCK JUST DEFINES NEW NAMES FOR UNITS IN TERMS OF EXISTING UNITS IN THE UNIX UNITS DATABASE
UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(umho) = (micromho)
}

: PARAMETERS ARE NORMALLY SPECIFIED BY THE USER
: PARAMETERS GENERALLY REMAIN CONSTANT BUT THEY CAN BE CHANGED IN MID-RUN IF NECESSARY
: THE PARAMETERS BLOCK SETS DEFAULT VALUES
: IF A PARAMETER DOES NOT APPEAR IN A NEURON BLOCKS RANGE STATEMENT, IT WILL HAVE GLOBAL SCOPE
: GLOBAL SCOPE MEANS CHANGING ITS VALUE WILL AFFECT EVERY INSTANCE OF THAT MECHANISM
PARAMETER {
	: e = -90 mV for inhibitory synapses,
	:     0 mV for excitatory
	e = -90	(mV)
	: tau_1 was the same for inhibitory and excitatory synapses
	: in the models used by T et al.
	tau_1 = 2 (ms) < 1e-9, 1e9 >
	tau_2 = 20
	: tau_rec = 100 ms for inhibitory synapses,
	:           800 ms for excitatory
	tau_rec = 100 (ms) < 1e-9, 1e9 >
	: tau_facil = 1000 ms for inhibitory synapses,
	:             0 ms for excitatory
	tau_facil = 1000 (ms) < 0, 1e9 >
	: U = 0.04 for inhibitory synapses, 
	:     0.5 for excitatory
	: the (1) is needed for the < 0, 1 > to be effective
	:   in limiting the values of U and u0
	U = 0.04 (1) < 0, 1 >
	: initial value for the "facilitation variable"
	u0 = 0 (1) < 0, 1 >
}

: THE ASSIGNED BLOCK IS USEDTO DECLARE TWO KINDS OF VARIABLES
: ONE TYPE THAT IS GIVEN VALUES OUTSIDE THE MOD FILE
: THE OTHER TYPE THAT APPEARS ON THE LEFT HAND SIDE OF AN ASSIGNMENT IN THE MOD FILE
: HERE v IS ASSIGNED OUTSIDE THE MOD FILE. IT COMES FROM THE SEGMENT THE MECHANISM IS ATTACHED TO
: i IS ASSIGNED IN THE MOD FILE. IT IS THE CURRENT THAT FLOWS AS A CONSEQUENCE OF THE CONDUCTANCE CHANGE
: AND THE DIFFERENCE BETWEEN v AND THE REVERSAL POTENTIAL
: x IS ALSO ASSIGNED IN THE MOD FILE
ASSIGNED {
	v (mV)
	i (nA)
	g (umho)
	x
	factor
}

: IF A MODEL INVOLVES DIFFERENTIAL EQUATIONS THEIR DEPENDENT VARIABLES OR UNKNOWNS ARE TO BE LISTED
: IN THE STATE BLOCK
STATE {
	A (uS)
	B (uS)

}

: THE INITIAL BLOCK IS CALLED BE finitialize() TO SET GOOD INITIAL CONDITIONS OF THE MECHANISM
INITIAL {
	LOCAL tp
	g=0
	if (tau_1/tau_2 > 0.9999) {
		tau_1 = 0.9999*tau_2
	}
	if (tau_1/tau_2 < 1e-9) {
		tau_1 = tau_2*1e-9
	}
	A = 0
	B = 0
	tp = (tau_1*tau_2)/(tau_2 - tau_1) * log(tau_2/tau_1)
	factor = -exp(-tp/tau_1) + exp(-tp/tau_2)
	factor = 1/factor
}

: BREAKPOINT IS THE MAIN COMPUTATION BLOCK OF THE MECHANISM
: BY THE END OF THE BLOCK ALLE VARIABLES ARE CONSISTENT WITH THE NEW TIME
: SOLVE SPECIFIES THE MAGIC BY WHICH THE STATE IS SOLVED
BREAKPOINT {
	SOLVE state METHOD cnexp
	g = B - A
	i = g*(v - e)
}

: THE DERIVATIVE STATE BLOCK ASSIGNES VALUES TO THE DERIVATIVES
DERIVATIVE state {
	A' = -A/tau_1
	B' = -B/tau_2
	: g' = -g/tau_2
}

: THE NET_RECEIVE BLOCK SPECIFIES WHAT HAPPENS TO A NET_RECIEVE EVENT FROM A NETCON OBJECT
NET_RECEIVE(weight (umho), y, z, u, tsyn (ms)) {
INITIAL {
: these are in NET_RECEIVE to be per-stream

	y = 0
	z = 0
:	u = 0
	u = u0
	tsyn = t
: this header will appear once per stream
: printf("t\t t-tsyn\t y\t z\t u\t newu\t g\t dg\t newg\t newy\n")
}

	: first calculate z at event-
	:   based on prior y and z
	z = z*exp(-(t - tsyn)/tau_rec)
	z = z + ( y*(exp(-(t - tsyn)/tau_1) - exp(-(t - tsyn)/tau_rec)) / ((tau_1/tau_rec)-1) )
	: now calc y at event-
	y = y*exp(-(t - tsyn)/tau_1)

	x = 1-y-z

	: calc u at event--
	if (tau_facil > 0) {
		u = u*exp(-(t - tsyn)/tau_facil)
	} else {
		u = U
	}

: printf("%g\t%g\t%g\t%g\t%g", t, t-tsyn, y, z, u)

	if (tau_facil > 0) {
		state_discontinuity(u, u + U*(1-u))
	}

: printf("\t%g\t%g\t%g", u, g, weight*factor*x*u)

	state_discontinuity(g, g + weight*factor*x*u)
	state_discontinuity(y, y + x*u)

	tsyn = t
	A = A + g
	B = B + g

: printf("\t%g\t%g\n", g, y)
}
