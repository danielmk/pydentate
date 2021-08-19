#include <stdio.h>
#include "hocdec.h"
#define IMPORT extern __declspec(dllimport)
IMPORT int nrnmpi_myid, nrn_nobanner_;

extern void _CaBK_reg();
extern void _Gfluct2_reg();
extern void _LcaMig_reg();
extern void _bgka_reg();
extern void _ccanl_reg();
extern void _gskch_reg();
extern void _hyperde3_reg();
extern void _ichan2_reg();
extern void _nca_reg();
extern void _netstim125_reg();
extern void _netstimbox_reg();
extern void _tca_reg();
extern void _tmgexp2syn_reg();
extern void _tmgsyn_reg();
extern void _vecevent_reg();

void modl_reg(){
	//nrn_mswindll_stdio(stdin, stdout, stderr);
    if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
	fprintf(stderr, "Additional mechanisms from files\n");

fprintf(stderr," CaBK.mod");
fprintf(stderr," Gfluct2.mod");
fprintf(stderr," LcaMig.mod");
fprintf(stderr," bgka.mod");
fprintf(stderr," ccanl.mod");
fprintf(stderr," gskch.mod");
fprintf(stderr," hyperde3.mod");
fprintf(stderr," ichan2.mod");
fprintf(stderr," nca.mod");
fprintf(stderr," netstim125.mod");
fprintf(stderr," netstimbox.mod");
fprintf(stderr," tca.mod");
fprintf(stderr," tmgexp2syn.mod");
fprintf(stderr," tmgsyn.mod");
fprintf(stderr," vecevent.mod");
fprintf(stderr, "\n");
    }
_CaBK_reg();
_Gfluct2_reg();
_LcaMig_reg();
_bgka_reg();
_ccanl_reg();
_gskch_reg();
_hyperde3_reg();
_ichan2_reg();
_nca_reg();
_netstim125_reg();
_netstimbox_reg();
_tca_reg();
_tmgexp2syn_reg();
_tmgsyn_reg();
_vecevent_reg();
}
