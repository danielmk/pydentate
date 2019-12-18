: $Id: netstim.mod 2212 2008-09-08 14:32:26Z hines $
: adapted by A. Hanuschkin 2011 -> output of a 'nspk' spikes in a given interval, if activated!
: netstim is activated by a NET_RECEIVE event. 
: comments at end

NEURON	{ 
  ARTIFICIAL_CELL NetStimBox
  RANGE start, forcestop, status, nspk
  THREADSAFE : only true if every instance has its own distinct Random
  POINTER donotuse
}

PARAMETER {
	start		= 50 (ms)	: start of first spike
	forcestop 	= 200 (ms)	: stop of firing spikes
	status 		= 0		: if status=0, no spike is sent
	nspk		= 1		: number of spikes per PP input	
}

ASSIGNED {
	event (ms)
	on
	ispike
	donotuse
}


INITIAL {			: deactivated by default
	on = 0  : off
}	

FUNCTION invl(mean (ms)) (ms) {				      
}	
VERBATIM
double nrn_random_pick(void* r);
void* nrn_random_arg(int argpos);
ENDVERBATIM

FUNCTION erand() {
VERBATIM
	if (_p_donotuse) {
		//  printf ("x");
		/*
		:Supports separate independent but reproducible streams for
		: each instance. However, the corresponding hoc Random
		: distribution MUST be set to Random.negexp(1)
		*/
		_lerand = nrn_random_pick(_p_donotuse);
	}else{
		/* only can be used in main thread */
		if (_nt != nrn_threads) {
hoc_execerror("multithread random in NetStim"," only via hoc Random");
		}
ENDVERBATIM
		: the old standby. Cannot use if reproducible parallel sim
		: independent of nhost or which host this instance is on
		: is desired, since each instance on this cpu draws from
		: the same stream
		erand = exprand(1)
VERBATIM
	}
ENDVERBATIM
}

PROCEDURE noiseFromRandom() {
VERBATIM
 {
	void** pv = (void**)(&_p_donotuse);
	if (ifarg(1)) {
		*pv = nrn_random_arg(1);
	}else{
		*pv = (void*)0;
	}
 }
ENDVERBATIM
}

NET_RECEIVE (w) {
      : printf ("NetStimBox: Net_Receive..\n")
 :    if (flag == 0) {            : if activated & external event -> return spike in the interval defined [start:forcestop]
      if (status == 1) {
	: printf ("NetStimBox: I'm active and got an input spike.... I'll answer with a random event....\n")
	FROM i=1 TO nspk {
  	      event = erand()*(forcestop-start)+start
        	: printf ("%f\t%f\t%f\t%f\t%f\n",erand(),erand(),erand(),erand(),erand())
        	if (event < 0) {
        	         event = 0
        	}
		printf ("NetStimBox: Send spike at: %f\n",event)
		net_event(event)
	}
	status = 0			: switch it off 
      } 
 :  }
}


COMMENT

ModelDB file along with publication:
Yim MY, Hanuschkin A, Wolfart J (2015) Hippocampus 25:297-308.
http://onlinelibrary.wiley.com/doi/10.1002/hipo.22373/abstract

ENDCOMMENT

