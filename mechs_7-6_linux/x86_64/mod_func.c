#include <stdio.h>
#include "hocdec.h"
extern int nrnmpi_myid;
extern int nrn_nobanner_;

extern void _CaBK_reg(void);
extern void _Gfluct2_reg(void);
extern void _LcaMig_reg(void);
extern void _bgka_reg(void);
extern void _ccanl_reg(void);
extern void _gskch_reg(void);
extern void _hyperde3_reg(void);
extern void _ichan2_reg(void);
extern void _nca_reg(void);
extern void _netstim125_reg(void);
extern void _netstimbox_reg(void);
extern void _tca_reg(void);
extern void _tmgsyn_reg(void);
extern void _vecevent_reg(void);

void modl_reg(){
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
  _tmgsyn_reg();
  _vecevent_reg();
}
