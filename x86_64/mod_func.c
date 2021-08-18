#include <stdio.h>
#include "hocdec.h"
extern int nrnmpi_myid;
extern int nrn_nobanner_;

extern void _bgka_reg(void);
extern void _CaBK_reg(void);
extern void _ccanl_reg(void);
extern void _Gfluct2_reg(void);
extern void _gskch_reg(void);
extern void _hyperde3_reg(void);
extern void _ichan2_reg(void);
extern void _LcaMig_reg(void);
extern void _nca_reg(void);
extern void _netstim125_reg(void);
extern void _netstimbox_reg(void);
extern void _tca_reg(void);
extern void _tmgexp2syn_reg(void);
extern void _tmgsyn_reg(void);
extern void _vecevent_reg(void);

void modl_reg(){
  if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
    fprintf(stderr, "Additional mechanisms from files\n");

    fprintf(stderr," \"mechs/bgka.mod\"");
    fprintf(stderr," \"mechs/CaBK.mod\"");
    fprintf(stderr," \"mechs/ccanl.mod\"");
    fprintf(stderr," \"mechs/Gfluct2.mod\"");
    fprintf(stderr," \"mechs/gskch.mod\"");
    fprintf(stderr," \"mechs/hyperde3.mod\"");
    fprintf(stderr," \"mechs/ichan2.mod\"");
    fprintf(stderr," \"mechs/LcaMig.mod\"");
    fprintf(stderr," \"mechs/nca.mod\"");
    fprintf(stderr," \"mechs/netstim125.mod\"");
    fprintf(stderr," \"mechs/netstimbox.mod\"");
    fprintf(stderr," \"mechs/tca.mod\"");
    fprintf(stderr," \"mechs/tmgexp2syn.mod\"");
    fprintf(stderr," \"mechs/tmgsyn.mod\"");
    fprintf(stderr," \"mechs/vecevent.mod\"");
    fprintf(stderr, "\n");
  }
  _bgka_reg();
  _CaBK_reg();
  _ccanl_reg();
  _Gfluct2_reg();
  _gskch_reg();
  _hyperde3_reg();
  _ichan2_reg();
  _LcaMig_reg();
  _nca_reg();
  _netstim125_reg();
  _netstimbox_reg();
  _tca_reg();
  _tmgexp2syn_reg();
  _tmgsyn_reg();
  _vecevent_reg();
}
