/* Created by Language version: 7.7.0 */
/* VECTORIZED */
#define NRN_VECTORIZED 1
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mech_api.h"
#undef PI
#define nil 0
#include "md1redef.h"
#include "section.h"
#include "nrniv_mf.h"
#include "md2redef.h"
 
#if METHOD3
extern int _method3;
#endif

#if !NRNGPU
#undef exp
#define exp hoc_Exp
extern double hoc_Exp(double);
#endif
 
#define nrn_init _nrn_init__hyperde3
#define _nrn_initial _nrn_initial__hyperde3
#define nrn_cur _nrn_cur__hyperde3
#define _nrn_current _nrn_current__hyperde3
#define nrn_jacob _nrn_jacob__hyperde3
#define nrn_state _nrn_state__hyperde3
#define _net_receive _net_receive__hyperde3 
#define _f_trates _f_trates__hyperde3 
#define rates rates__hyperde3 
#define states states__hyperde3 
#define trates trates__hyperde3 
 
#define _threadargscomma_ _p, _ppvar, _thread, _nt,
#define _threadargsprotocomma_ double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt,
#define _threadargs_ _p, _ppvar, _thread, _nt
#define _threadargsproto_ double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt
 	/*SUPPRESS 761*/
	/*SUPPRESS 762*/
	/*SUPPRESS 763*/
	/*SUPPRESS 765*/
	 extern double *getarg();
 /* Thread safe. No static _p or _ppvar. */
 
#define t _nt->_t
#define dt _nt->_dt
#define ghyfbar _p[0]
#define ghyfbar_columnindex 0
#define ghysbar _p[1]
#define ghysbar_columnindex 1
#define ghyhtfbar _p[2]
#define ghyhtfbar_columnindex 2
#define ghyhtsbar _p[3]
#define ghyhtsbar_columnindex 3
#define ghyf _p[4]
#define ghyf_columnindex 4
#define ghys _p[5]
#define ghys_columnindex 5
#define ghyhtf _p[6]
#define ghyhtf_columnindex 6
#define ghyhts _p[7]
#define ghyhts_columnindex 7
#define ihyf _p[8]
#define ihyf_columnindex 8
#define ihys _p[9]
#define ihys_columnindex 9
#define hyfinf _p[10]
#define hyfinf_columnindex 10
#define hysinf _p[11]
#define hysinf_columnindex 11
#define hyhtfinf _p[12]
#define hyhtfinf_columnindex 12
#define hyhtsinf _p[13]
#define hyhtsinf_columnindex 13
#define hyftau _p[14]
#define hyftau_columnindex 14
#define hystau _p[15]
#define hystau_columnindex 15
#define hyhtftau _p[16]
#define hyhtftau_columnindex 16
#define hyhtstau _p[17]
#define hyhtstau_columnindex 17
#define hyf _p[18]
#define hyf_columnindex 18
#define hys _p[19]
#define hys_columnindex 19
#define hyhtf _p[20]
#define hyhtf_columnindex 20
#define hyhts _p[21]
#define hyhts_columnindex 21
#define ehyf _p[22]
#define ehyf_columnindex 22
#define ehys _p[23]
#define ehys_columnindex 23
#define ehyhtf _p[24]
#define ehyhtf_columnindex 24
#define ehyhts _p[25]
#define ehyhts_columnindex 25
#define Dhyf _p[26]
#define Dhyf_columnindex 26
#define Dhys _p[27]
#define Dhys_columnindex 27
#define Dhyhtf _p[28]
#define Dhyhtf_columnindex 28
#define Dhyhts _p[29]
#define Dhyhts_columnindex 29
#define ihyhtf _p[30]
#define ihyhtf_columnindex 30
#define ihyhts _p[31]
#define ihyhts_columnindex 31
#define hyfexp _p[32]
#define hyfexp_columnindex 32
#define hysexp _p[33]
#define hysexp_columnindex 33
#define hyhtfexp _p[34]
#define hyhtfexp_columnindex 34
#define hyhtsexp _p[35]
#define hyhtsexp_columnindex 35
#define v _p[36]
#define v_columnindex 36
#define _g _p[37]
#define _g_columnindex 37
#define _ion_ehyf	*_ppvar[0]._pval
#define _ion_ihyf	*_ppvar[1]._pval
#define _ion_dihyfdv	*_ppvar[2]._pval
#define _ion_ehys	*_ppvar[3]._pval
#define _ion_ihys	*_ppvar[4]._pval
#define _ion_dihysdv	*_ppvar[5]._pval
#define _ion_ehyhtf	*_ppvar[6]._pval
#define _ion_ihyhtf	*_ppvar[7]._pval
#define _ion_dihyhtfdv	*_ppvar[8]._pval
#define _ion_ehyhts	*_ppvar[9]._pval
#define _ion_ihyhts	*_ppvar[10]._pval
#define _ion_dihyhtsdv	*_ppvar[11]._pval
 
#if MAC
#if !defined(v)
#define v _mlhv
#endif
#if !defined(h)
#define h _mlhh
#endif
#endif
 
#if defined(__cplusplus)
extern "C" {
#endif
 static int hoc_nrnpointerindex =  -1;
 static Datum* _extcall_thread;
 static Prop* _extcall_prop;
 /* external NEURON variables */
 extern double celsius;
 /* declaration of user functions */
 static void _hoc_rates(void);
 static void _hoc_states(void);
 static void _hoc_trates(void);
 static void _hoc_vtrap(void);
 static int _mechtype;
extern void _nrn_cacheloop_reg(int, int);
extern void hoc_register_prop_size(int, int, int);
extern void hoc_register_limits(int, HocParmLimits*);
extern void hoc_register_units(int, HocParmUnits*);
extern void nrn_promote(Prop*, int, int);
extern Memb_func* memb_func;
 
#define NMODL_TEXT 1
#if NMODL_TEXT
static const char* nmodl_file_text;
static const char* nmodl_filename;
extern void hoc_reg_nmodl_text(int, const char*);
extern void hoc_reg_nmodl_filename(int, const char*);
#endif

 extern void _nrn_setdata_reg(int, void(*)(Prop*));
 static void _setdata(Prop* _prop) {
 _extcall_prop = _prop;
 }
 static void _hoc_setdata() {
 Prop *_prop, *hoc_getdata_range(int);
 _prop = hoc_getdata_range(_mechtype);
   _setdata(_prop);
 hoc_retpushx(1.);
}
 /* connect user functions to hoc names */
 static VoidFunc hoc_intfunc[] = {
 "setdata_hyperde3", _hoc_setdata,
 "rates_hyperde3", _hoc_rates,
 "states_hyperde3", _hoc_states,
 "trates_hyperde3", _hoc_trates,
 "vtrap_hyperde3", _hoc_vtrap,
 0, 0
};
#define vtrap vtrap_hyperde3
 extern double vtrap( _threadargsprotocomma_ double , double );
 
static void _check_trates(double*, Datum*, Datum*, NrnThread*); 
static void _check_table_thread(double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt, int _type) {
   _check_trates(_p, _ppvar, _thread, _nt);
 }
 #define _zq10 _thread[0]._pval[0]
 /* declare global and static user variables */
#define usetable usetable_hyperde3
 double usetable = 1;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 "usetable_hyperde3", 0, 1,
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "ghyfbar_hyperde3", "mho/cm2",
 "ghysbar_hyperde3", "mho/cm2",
 "ghyhtfbar_hyperde3", "mho/cm2",
 "ghyhtsbar_hyperde3", "mho/cm2",
 "ghyf_hyperde3", "mho/cm2",
 "ghys_hyperde3", "mho/cm2",
 "ghyhtf_hyperde3", "mho/cm2",
 "ghyhts_hyperde3", "mho/cm2",
 "ihyf_hyperde3", "mA/cm2",
 "ihys_hyperde3", "mA/cm2",
 "hyftau_hyperde3", "ms",
 "hystau_hyperde3", "ms",
 "hyhtftau_hyperde3", "ms",
 "hyhtstau_hyperde3", "ms",
 0,0
};
 static double delta_t = 1;
 static double hyhts0 = 0;
 static double hyhtf0 = 0;
 static double hys0 = 0;
 static double hyf0 = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "usetable_hyperde3", &usetable_hyperde3,
 0,0
};
 static DoubVec hoc_vdoub[] = {
 0,0,0
};
 static double _sav_indep;
 static void nrn_alloc(Prop*);
static void  nrn_init(NrnThread*, _Memb_list*, int);
static void nrn_state(NrnThread*, _Memb_list*, int);
 static void nrn_cur(NrnThread*, _Memb_list*, int);
static void  nrn_jacob(NrnThread*, _Memb_list*, int);
 
static int _ode_count(int);
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"hyperde3",
 "ghyfbar_hyperde3",
 "ghysbar_hyperde3",
 "ghyhtfbar_hyperde3",
 "ghyhtsbar_hyperde3",
 0,
 "ghyf_hyperde3",
 "ghys_hyperde3",
 "ghyhtf_hyperde3",
 "ghyhts_hyperde3",
 "ihyf_hyperde3",
 "ihys_hyperde3",
 "hyfinf_hyperde3",
 "hysinf_hyperde3",
 "hyhtfinf_hyperde3",
 "hyhtsinf_hyperde3",
 "hyftau_hyperde3",
 "hystau_hyperde3",
 "hyhtftau_hyperde3",
 "hyhtstau_hyperde3",
 0,
 "hyf_hyperde3",
 "hys_hyperde3",
 "hyhtf_hyperde3",
 "hyhts_hyperde3",
 0,
 0};
 static Symbol* _hyf_sym;
 static Symbol* _hys_sym;
 static Symbol* _hyhtf_sym;
 static Symbol* _hyhts_sym;
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 38, _prop);
 	/*initialize range parameters*/
 	ghyfbar = 0;
 	ghysbar = 0;
 	ghyhtfbar = 0;
 	ghyhtsbar = 0;
 	_prop->param = _p;
 	_prop->param_size = 38;
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 12, _prop);
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 prop_ion = need_memb(_hyf_sym);
 nrn_promote(prop_ion, 0, 1);
 	_ppvar[0]._pval = &prop_ion->param[0]; /* ehyf */
 	_ppvar[1]._pval = &prop_ion->param[3]; /* ihyf */
 	_ppvar[2]._pval = &prop_ion->param[4]; /* _ion_dihyfdv */
 prop_ion = need_memb(_hys_sym);
 nrn_promote(prop_ion, 0, 1);
 	_ppvar[3]._pval = &prop_ion->param[0]; /* ehys */
 	_ppvar[4]._pval = &prop_ion->param[3]; /* ihys */
 	_ppvar[5]._pval = &prop_ion->param[4]; /* _ion_dihysdv */
 prop_ion = need_memb(_hyhtf_sym);
 nrn_promote(prop_ion, 0, 1);
 	_ppvar[6]._pval = &prop_ion->param[0]; /* ehyhtf */
 	_ppvar[7]._pval = &prop_ion->param[3]; /* ihyhtf */
 	_ppvar[8]._pval = &prop_ion->param[4]; /* _ion_dihyhtfdv */
 prop_ion = need_memb(_hyhts_sym);
 nrn_promote(prop_ion, 0, 1);
 	_ppvar[9]._pval = &prop_ion->param[0]; /* ehyhts */
 	_ppvar[10]._pval = &prop_ion->param[3]; /* ihyhts */
 	_ppvar[11]._pval = &prop_ion->param[4]; /* _ion_dihyhtsdv */
 
}
 static void _initlists();
 static void _thread_mem_init(Datum*);
 static void _thread_cleanup(Datum*);
 static void _update_ion_pointer(Datum*);
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
extern void _nrn_thread_table_reg(int, void(*)(double*, Datum*, Datum*, NrnThread*, int));
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 void _hyperde3_reg() {
	int _vectorized = 1;
  _initlists();
 	ion_reg("hyf", 1.0);
 	ion_reg("hys", 1.0);
 	ion_reg("hyhtf", 1.0);
 	ion_reg("hyhts", 1.0);
 	_hyf_sym = hoc_lookup("hyf_ion");
 	_hys_sym = hoc_lookup("hys_ion");
 	_hyhtf_sym = hoc_lookup("hyhtf_ion");
 	_hyhts_sym = hoc_lookup("hyhts_ion");
 	register_mech(_mechanism, nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init, hoc_nrnpointerindex, 2);
  _extcall_thread = (Datum*)ecalloc(1, sizeof(Datum));
  _thread_mem_init(_extcall_thread);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
     _nrn_thread_reg(_mechtype, 1, _thread_mem_init);
     _nrn_thread_reg(_mechtype, 0, _thread_cleanup);
     _nrn_thread_reg(_mechtype, 2, _update_ion_pointer);
     _nrn_thread_table_reg(_mechtype, _check_table_thread);
 #if NMODL_TEXT
  hoc_reg_nmodl_text(_mechtype, nmodl_file_text);
  hoc_reg_nmodl_filename(_mechtype, nmodl_filename);
#endif
  hoc_register_prop_size(_mechtype, 38, 12);
  hoc_register_dparam_semantics(_mechtype, 0, "hyf_ion");
  hoc_register_dparam_semantics(_mechtype, 1, "hyf_ion");
  hoc_register_dparam_semantics(_mechtype, 2, "hyf_ion");
  hoc_register_dparam_semantics(_mechtype, 3, "hys_ion");
  hoc_register_dparam_semantics(_mechtype, 4, "hys_ion");
  hoc_register_dparam_semantics(_mechtype, 5, "hys_ion");
  hoc_register_dparam_semantics(_mechtype, 6, "hyhtf_ion");
  hoc_register_dparam_semantics(_mechtype, 7, "hyhtf_ion");
  hoc_register_dparam_semantics(_mechtype, 8, "hyhtf_ion");
  hoc_register_dparam_semantics(_mechtype, 9, "hyhts_ion");
  hoc_register_dparam_semantics(_mechtype, 10, "hyhts_ion");
  hoc_register_dparam_semantics(_mechtype, 11, "hyhts_ion");
 	hoc_register_cvode(_mechtype, _ode_count, 0, 0, 0);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 hyperde3 /Users/temma/ghq/pydentate/mechs/hyperde3.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
 static double FARADAY = 96520.0;
 static double R = 8.3134;
 /*Top LOCAL _zq10 */
 static double *_t_hyfinf;
 static double *_t_hyhtfinf;
 static double *_t_hyfexp;
 static double *_t_hyhtfexp;
 static double *_t_hyftau;
 static double *_t_hyhtftau;
 static double *_t_hysinf;
 static double *_t_hyhtsinf;
 static double *_t_hysexp;
 static double *_t_hyhtsexp;
 static double *_t_hystau;
 static double *_t_hyhtstau;
static int _reset;
static char *modelname = "hyperde3.mod  ";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static int _f_trates(_threadargsprotocomma_ double);
static int rates(_threadargsprotocomma_ double);
static int states(_threadargsproto_);
static int trates(_threadargsprotocomma_ double);
 static void _n_trates(_threadargsprotocomma_ double _lv);
 
static int  states ( _threadargsproto_ ) {
   trates ( _threadargscomma_ v ) ;
   hyf = hyf + hyfexp * ( hyfinf - hyf ) ;
   hys = hys + hysexp * ( hysinf - hys ) ;
   hyhtf = hyhtf + hyhtfexp * ( hyhtfinf - hyhtf ) ;
   hyhts = hyhts + hyhtsexp * ( hyhtsinf - hyhts ) ;
    return 0; }
 
static void _hoc_states(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r = 1.;
 states ( _p, _ppvar, _thread, _nt );
 hoc_retpushx(_r);
}
 
static int  rates ( _threadargsprotocomma_ double _lv ) {
   double _lalpha , _lbeta , _lsum ;
 _zq10 = pow( 3.0 , ( ( celsius - 6.3 ) / 10.0 ) ) ;
   hyfinf = 1.0 / ( 1.0 + exp ( ( _lv + 91.0 ) / 10.0 ) ) ;
   hyftau = 14.9 + 14.1 / ( 1.0 + exp ( - ( _lv + 95.2 ) / 0.5 ) ) ;
   hysinf = 1.0 / ( 1.0 + exp ( ( _lv + 91.0 ) / 10.0 ) ) ;
   hystau = 80.0 + 172.7 / ( 1.0 + exp ( - ( _lv + 59.3 ) / - 0.83 ) ) ;
   hyhtfinf = 1.0 / ( 1.0 + exp ( ( _lv + 87.0 ) / 10.0 ) ) ;
   hyhtftau = 23.2 + 16.1 / ( 1.0 + exp ( - ( _lv + 91.2 ) / 0.83 ) ) ;
   hyhtsinf = 1.0 / ( 1.0 + exp ( ( _lv + 87.0 ) / 10.0 ) ) ;
   hyhtstau = 227.3 + 170.7 * exp ( - 0.5 * pow( ( ( _lv + 80.4 ) / 11.0 ) , 2.0 ) ) ;
    return 0; }
 
static void _hoc_rates(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r = 1.;
 rates ( _p, _ppvar, _thread, _nt, *getarg(1) );
 hoc_retpushx(_r);
}
 static double _mfac_trates, _tmin_trates;
  static void _check_trates(double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt) {
  static int _maktable=1; int _i, _j, _ix = 0;
  double _xi, _tmax;
  static double _sav_dt;
  static double _sav_celsius;
  if (!usetable) {return;}
  if (_sav_dt != dt) { _maktable = 1;}
  if (_sav_celsius != celsius) { _maktable = 1;}
  if (_maktable) { double _x, _dx; _maktable=0;
   _tmin_trates =  - 120.0 ;
   _tmax =  100.0 ;
   _dx = (_tmax - _tmin_trates)/220.; _mfac_trates = 1./_dx;
   for (_i=0, _x=_tmin_trates; _i < 221; _x += _dx, _i++) {
    _f_trates(_p, _ppvar, _thread, _nt, _x);
    _t_hyfinf[_i] = hyfinf;
    _t_hyhtfinf[_i] = hyhtfinf;
    _t_hyfexp[_i] = hyfexp;
    _t_hyhtfexp[_i] = hyhtfexp;
    _t_hyftau[_i] = hyftau;
    _t_hyhtftau[_i] = hyhtftau;
    _t_hysinf[_i] = hysinf;
    _t_hyhtsinf[_i] = hyhtsinf;
    _t_hysexp[_i] = hysexp;
    _t_hyhtsexp[_i] = hyhtsexp;
    _t_hystau[_i] = hystau;
    _t_hyhtstau[_i] = hyhtstau;
   }
   _sav_dt = dt;
   _sav_celsius = celsius;
  }
 }

 static int trates(double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt, double _lv) { 
#if 0
_check_trates(_p, _ppvar, _thread, _nt);
#endif
 _n_trates(_p, _ppvar, _thread, _nt, _lv);
 return 0;
 }

 static void _n_trates(double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt, double _lv){ int _i, _j;
 double _xi, _theta;
 if (!usetable) {
 _f_trates(_p, _ppvar, _thread, _nt, _lv); return; 
}
 _xi = _mfac_trates * (_lv - _tmin_trates);
 if (isnan(_xi)) {
  hyfinf = _xi;
  hyhtfinf = _xi;
  hyfexp = _xi;
  hyhtfexp = _xi;
  hyftau = _xi;
  hyhtftau = _xi;
  hysinf = _xi;
  hyhtsinf = _xi;
  hysexp = _xi;
  hyhtsexp = _xi;
  hystau = _xi;
  hyhtstau = _xi;
  return;
 }
 if (_xi <= 0.) {
 hyfinf = _t_hyfinf[0];
 hyhtfinf = _t_hyhtfinf[0];
 hyfexp = _t_hyfexp[0];
 hyhtfexp = _t_hyhtfexp[0];
 hyftau = _t_hyftau[0];
 hyhtftau = _t_hyhtftau[0];
 hysinf = _t_hysinf[0];
 hyhtsinf = _t_hyhtsinf[0];
 hysexp = _t_hysexp[0];
 hyhtsexp = _t_hyhtsexp[0];
 hystau = _t_hystau[0];
 hyhtstau = _t_hyhtstau[0];
 return; }
 if (_xi >= 220.) {
 hyfinf = _t_hyfinf[220];
 hyhtfinf = _t_hyhtfinf[220];
 hyfexp = _t_hyfexp[220];
 hyhtfexp = _t_hyhtfexp[220];
 hyftau = _t_hyftau[220];
 hyhtftau = _t_hyhtftau[220];
 hysinf = _t_hysinf[220];
 hyhtsinf = _t_hyhtsinf[220];
 hysexp = _t_hysexp[220];
 hyhtsexp = _t_hyhtsexp[220];
 hystau = _t_hystau[220];
 hyhtstau = _t_hyhtstau[220];
 return; }
 _i = (int) _xi;
 _theta = _xi - (double)_i;
 hyfinf = _t_hyfinf[_i] + _theta*(_t_hyfinf[_i+1] - _t_hyfinf[_i]);
 hyhtfinf = _t_hyhtfinf[_i] + _theta*(_t_hyhtfinf[_i+1] - _t_hyhtfinf[_i]);
 hyfexp = _t_hyfexp[_i] + _theta*(_t_hyfexp[_i+1] - _t_hyfexp[_i]);
 hyhtfexp = _t_hyhtfexp[_i] + _theta*(_t_hyhtfexp[_i+1] - _t_hyhtfexp[_i]);
 hyftau = _t_hyftau[_i] + _theta*(_t_hyftau[_i+1] - _t_hyftau[_i]);
 hyhtftau = _t_hyhtftau[_i] + _theta*(_t_hyhtftau[_i+1] - _t_hyhtftau[_i]);
 hysinf = _t_hysinf[_i] + _theta*(_t_hysinf[_i+1] - _t_hysinf[_i]);
 hyhtsinf = _t_hyhtsinf[_i] + _theta*(_t_hyhtsinf[_i+1] - _t_hyhtsinf[_i]);
 hysexp = _t_hysexp[_i] + _theta*(_t_hysexp[_i+1] - _t_hysexp[_i]);
 hyhtsexp = _t_hyhtsexp[_i] + _theta*(_t_hyhtsexp[_i+1] - _t_hyhtsexp[_i]);
 hystau = _t_hystau[_i] + _theta*(_t_hystau[_i+1] - _t_hystau[_i]);
 hyhtstau = _t_hyhtstau[_i] + _theta*(_t_hyhtstau[_i+1] - _t_hyhtstau[_i]);
 }

 
static int  _f_trates ( _threadargsprotocomma_ double _lv ) {
   double _ltinc ;
 rates ( _threadargscomma_ _lv ) ;
   _ltinc = - dt * _zq10 ;
   hyfexp = 1.0 - exp ( _ltinc / hyftau ) ;
   hysexp = 1.0 - exp ( _ltinc / hystau ) ;
   hyhtfexp = 1.0 - exp ( _ltinc / hyhtftau ) ;
   hyhtsexp = 1.0 - exp ( _ltinc / hyhtstau ) ;
    return 0; }
 
static void _hoc_trates(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 
#if 1
 _check_trates(_p, _ppvar, _thread, _nt);
#endif
 _r = 1.;
 trates ( _p, _ppvar, _thread, _nt, *getarg(1) );
 hoc_retpushx(_r);
}
 
double vtrap ( _threadargsprotocomma_ double _lx , double _ly ) {
   double _lvtrap;
 if ( fabs ( _lx / _ly ) < 1e-6 ) {
     _lvtrap = _ly * ( 1.0 - _lx / _ly / 2.0 ) ;
     }
   else {
     _lvtrap = _lx / ( exp ( _lx / _ly ) - 1.0 ) ;
     }
   
return _lvtrap;
 }
 
static void _hoc_vtrap(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  vtrap ( _p, _ppvar, _thread, _nt, *getarg(1) , *getarg(2) );
 hoc_retpushx(_r);
}
 
static int _ode_count(int _type){ hoc_execerror("hyperde3", "cannot be used with CVODE"); return 0;}
 
static void _thread_mem_init(Datum* _thread) {
   _thread[0]._pval = (double*)ecalloc(1, sizeof(double));
 }
 
static void _thread_cleanup(Datum* _thread) {
   free((void*)(_thread[0]._pval));
 }
 extern void nrn_update_ion_pointer(Symbol*, Datum*, int, int);
 static void _update_ion_pointer(Datum* _ppvar) {
   nrn_update_ion_pointer(_hyf_sym, _ppvar, 0, 0);
   nrn_update_ion_pointer(_hyf_sym, _ppvar, 1, 3);
   nrn_update_ion_pointer(_hyf_sym, _ppvar, 2, 4);
   nrn_update_ion_pointer(_hys_sym, _ppvar, 3, 0);
   nrn_update_ion_pointer(_hys_sym, _ppvar, 4, 3);
   nrn_update_ion_pointer(_hys_sym, _ppvar, 5, 4);
   nrn_update_ion_pointer(_hyhtf_sym, _ppvar, 6, 0);
   nrn_update_ion_pointer(_hyhtf_sym, _ppvar, 7, 3);
   nrn_update_ion_pointer(_hyhtf_sym, _ppvar, 8, 4);
   nrn_update_ion_pointer(_hyhts_sym, _ppvar, 9, 0);
   nrn_update_ion_pointer(_hyhts_sym, _ppvar, 10, 3);
   nrn_update_ion_pointer(_hyhts_sym, _ppvar, 11, 4);
 }

static void initmodel(double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt) {
  int _i; double _save;{
  hyhts = hyhts0;
  hyhtf = hyhtf0;
  hys = hys0;
  hyf = hyf0;
 {
   trates ( _threadargscomma_ v ) ;
   hyf = hyfinf ;
   hys = hysinf ;
   hyhtf = hyhtfinf ;
   hyhts = hyhtsinf ;
   }
 
}
}

static void nrn_init(NrnThread* _nt, _Memb_list* _ml, int _type){
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; double _v; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];

#if 0
 _check_trates(_p, _ppvar, _thread, _nt);
#endif
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v = _v;
  ehyf = _ion_ehyf;
  ehys = _ion_ehys;
  ehyhtf = _ion_ehyhtf;
  ehyhts = _ion_ehyhts;
 initmodel(_p, _ppvar, _thread, _nt);
    }
}

static double _nrn_current(double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt, double _v){double _current=0.;v=_v;{ {
   ghyf = ghyfbar * hyf * hyf ;
   ihyf = ghyf * ( v - ehyf ) ;
   ghys = ghysbar * hys * hys ;
   ihys = ghys * ( v - ehys ) ;
   ghyhtf = ghyhtfbar * hyhtf * hyhtf ;
   ihyhtf = ghyhtf * ( v - ehyhtf ) ;
   ghyhts = ghyhtsbar * hyhts * hyhts ;
   ihyhts = ghyhts * ( v - ehyhts ) ;
   }
 _current += ihyf;
 _current += ihys;
 _current += ihyhtf;
 _current += ihyhts;

} return _current;
}

static void nrn_cur(NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; int* _ni; double _rhs, _v; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
  ehyf = _ion_ehyf;
  ehys = _ion_ehys;
  ehyhtf = _ion_ehyhtf;
  ehyhts = _ion_ehyhts;
 _g = _nrn_current(_p, _ppvar, _thread, _nt, _v + .001);
 	{ double _dihyhts;
 double _dihyhtf;
 double _dihys;
 double _dihyf;
  _dihyf = ihyf;
  _dihys = ihys;
  _dihyhtf = ihyhtf;
  _dihyhts = ihyhts;
 _rhs = _nrn_current(_p, _ppvar, _thread, _nt, _v);
  _ion_dihyfdv += (_dihyf - ihyf)/.001 ;
  _ion_dihysdv += (_dihys - ihys)/.001 ;
  _ion_dihyhtfdv += (_dihyhtf - ihyhtf)/.001 ;
  _ion_dihyhtsdv += (_dihyhts - ihyhts)/.001 ;
 	}
 _g = (_g - _rhs)/.001;
  _ion_ihyf += ihyf ;
  _ion_ihys += ihys ;
  _ion_ihyhtf += ihyhtf ;
  _ion_ihyhts += ihyhts ;
#if CACHEVEC
  if (use_cachevec) {
	VEC_RHS(_ni[_iml]) -= _rhs;
  }else
#endif
  {
	NODERHS(_nd) -= _rhs;
  }
 
}
 
}

static void nrn_jacob(NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml];
#if CACHEVEC
  if (use_cachevec) {
	VEC_D(_ni[_iml]) += _g;
  }else
#endif
  {
     _nd = _ml->_nodelist[_iml];
	NODED(_nd) += _g;
  }
 
}
 
}

static void nrn_state(NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; double _v = 0.0; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
 _nd = _ml->_nodelist[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v=_v;
{
  ehyf = _ion_ehyf;
  ehys = _ion_ehys;
  ehyhtf = _ion_ehyhtf;
  ehyhts = _ion_ehyhts;
 {  { states(_p, _ppvar, _thread, _nt); }
  }    }}

}

static void terminal(){}

static void _initlists(){
 double _x; double* _p = &_x;
 int _i; static int _first = 1;
  if (!_first) return;
   _t_hyfinf = makevector(221*sizeof(double));
   _t_hyhtfinf = makevector(221*sizeof(double));
   _t_hyfexp = makevector(221*sizeof(double));
   _t_hyhtfexp = makevector(221*sizeof(double));
   _t_hyftau = makevector(221*sizeof(double));
   _t_hyhtftau = makevector(221*sizeof(double));
   _t_hysinf = makevector(221*sizeof(double));
   _t_hyhtsinf = makevector(221*sizeof(double));
   _t_hysexp = makevector(221*sizeof(double));
   _t_hyhtsexp = makevector(221*sizeof(double));
   _t_hystau = makevector(221*sizeof(double));
   _t_hyhtstau = makevector(221*sizeof(double));
_first = 0;
}

#if defined(__cplusplus)
} /* extern "C" */
#endif

#if NMODL_TEXT
static const char* nmodl_filename = "/Users/temma/ghq/pydentate/mechs/hyperde3.mod";
static const char* nmodl_file_text = 
  "TITLE hyperde3.mod  \n"
  " \n"
  "COMMENT\n"
  "Chen K, Aradi I, Thon N, Eghbal-Ahmadi M, Baram TZ, Soltesz I: Persistently\n"
  "modified\n"
  "h-channels after complex febrile seizures convert the seizure-induced\n"
  "enhancement of\n"
  "inhibition to hyperexcitability. Nature Medicine, 7(3) pp. 331-337, 2001.\n"
  "(modeling by Ildiko Aradi, iaradi@uci.edu)\n"
  "distal dendritic Ih channel kinetics for both HT and Control anlimals\n"
  "ENDCOMMENT\n"
  " \n"
  "UNITS {\n"
  "        (mA) =(milliamp)\n"
  "        (mV) =(millivolt)\n"
  "        (uF) = (microfarad)\n"
  "	(molar) = (1/liter)\n"
  "	(nA) = (nanoamp)\n"
  "	(mM) = (millimolar)\n"
  "	(um) = (micron)\n"
  "	FARADAY = 96520 (coul)\n"
  "	R = 8.3134	(joule/degC)\n"
  "}\n"
  " \n"
  "? interface \n"
  "NEURON { \n"
  "SUFFIX hyperde3 \n"
  "USEION hyf READ ehyf WRITE ihyf VALENCE 1\n"
  "USEION hys READ ehys WRITE ihys VALENCE 1\n"
  "USEION hyhtf READ ehyhtf WRITE ihyhtf VALENCE 1\n"
  "USEION hyhts READ ehyhts WRITE ihyhts VALENCE 1\n"
  "RANGE  ghyf, ghys, ghyhtf, ghyhts\n"
  "RANGE ghyfbar, ghysbar, ghyhtfbar, ghyhtsbar\n"
  "RANGE hyfinf, hysinf, hyftau, hystau\n"
  "RANGE hyhtfinf, hyhtsinf, hyhtftau, hyhtstau, ihyf, ihys\n"
  "}\n"
  " \n"
  "INDEPENDENT {t FROM 0 TO 100 WITH 100 (ms)}\n"
  " \n"
  "PARAMETER {\n"
  "      v (mV) \n"
  "      celsius = 6.3 (degC)\n"
  "      dt (ms) \n"
  "\n"
  "	ghyfbar (mho/cm2)\n"
  "	ghysbar (mho/cm2)\n"
  "	ehyf (mV)\n"
  "	ehys (mV)\n"
  "	ghyhtfbar (mho/cm2)\n"
  "	ghyhtsbar (mho/cm2)\n"
  "	ehyhtf (mV)\n"
  "	ehyhts (mV)\n"
  "}\n"
  " \n"
  "STATE {\n"
  "	hyf hys hyhtf hyhts\n"
  "}\n"
  " \n"
  "ASSIGNED {\n"
  "         \n"
  "  \n"
  "	ghyf (mho/cm2)\n"
  " 	ghys (mho/cm2)\n"
  "\n"
  "	ghyhtf (mho/cm2)\n"
  "	ghyhts (mho/cm2)\n"
  "\n"
  "  \n"
  "	ihyf (mA/cm2)\n"
  "	ihys (mA/cm2)\n"
  "	ihyhtf (mA/cm2)\n"
  "	ihyhts (mA/cm2)\n"
  "\n"
  "	hyfinf hysinf hyhtfinf hyhtsinf\n"
  " 	hyftau (ms) hystau (ms) hyhtftau (ms) hyhtstau (ms)\n"
  "	hyfexp hysexp hyhtfexp hyhtsexp     \n"
  "} \n"
  "\n"
  "? currents\n"
  "BREAKPOINT {\n"
  "\n"
  "	SOLVE states\n"
  "\n"
  "	ghyf = ghyfbar * hyf*hyf\n"
  "	ihyf = ghyf * (v-ehyf)\n"
  "	ghys = ghysbar * hys*hys\n"
  "	ihys = ghys * (v-ehys)\n"
  "\n"
  "	ghyhtf = ghyhtfbar * hyhtf* hyhtf\n"
  "	ihyhtf = ghyhtf * (v-ehyhtf)\n"
  "	ghyhts = ghyhtsbar * hyhts* hyhts\n"
  "	ihyhts = ghyhts * (v-ehyhts)\n"
  "		\n"
  "		}\n"
  " \n"
  "UNITSOFF\n"
  " \n"
  "INITIAL {\n"
  "	trates(v)\n"
  "	\n"
  "	hyf = hyfinf\n"
  "      hys = hysinf\n"
  "	hyhtf = hyhtfinf\n"
  "	hyhts = hyhtsinf\n"
  "}\n"
  "\n"
  "? states\n"
  "PROCEDURE states() {	:Computes state variables m, h, and n \n"
  "        trates(v)	:      at the current v and dt.\n"
  "        \n"
  "        hyf = hyf + hyfexp*(hyfinf-hyf)\n"
  "        hys = hys + hysexp*(hysinf-hys)\n"
  "	  hyhtf = hyhtf + hyhtfexp*(hyhtfinf-hyhtf)\n"
  "	  hyhts = hyhts + hyhtsexp*(hyhtsinf-hyhts)\n"
  "\n"
  "}\n"
  " \n"
  "LOCAL q10\n"
  "\n"
  "? rates\n"
  "PROCEDURE rates(v) {  :Computes rate and other constants at current v.\n"
  "                      :Call once from HOC to initialize inf at resting v.\n"
  "        LOCAL  alpha, beta, sum\n"
  "       q10 = 3^((celsius - 6.3)/10)\n"
  "       \n"
  "	:\"hyf\" FAST CONTROL Hype activation system\n"
  "	hyfinf =  1 / (1 + exp( (v+91)/10 ))\n"
  "	hyftau = 14.9 + 14.1 / (1+exp(-(v+95.2)/0.5))\n"
  "\n"
  "	:\"hys\" SLOW CONTROL Hype activation system\n"
  "	hysinf =  1 / (1 + exp( (v+91)/10 ))\n"
  "	hystau = 80 + 172.7 / (1+exp(-(v+59.3)/-0.83))\n"
  "\n"
  "		:\"hyhtf\" FAST HT Hypeht activation system\n"
  "	hyhtfinf =  1 / (1 + exp( (v+87)/10 ))\n"
  "	hyhtftau = 23.2 + 16.1 / (1+exp(-(v+91.2)/0.83))\n"
  "\n"
  "		:\"hyhts\" SLOW HT Hypeht activation system\n"
  "	hyhtsinf =  1 / (1 + exp( (v+87)/10 ))\n"
  "	hyhtstau = 227.3 + 170.7*exp(-0.5*((v+80.4)/11)^2)\n"
  "}\n"
  " \n"
  "PROCEDURE trates(v) {  :Computes rate and other constants at current v.\n"
  "                      :Call once from HOC to initialize inf at resting v.\n"
  "	LOCAL tinc\n"
  "      TABLE hyfinf, hyhtfinf, hyfexp, hyhtfexp, hyftau, hyhtftau, \n"
  "		hysinf, hyhtsinf, hysexp, hyhtsexp, hystau, hyhtstau	\n"
  "	DEPEND dt, celsius FROM -120 TO 100 WITH 220\n"
  "                           \n"
  "	rates(v)	: not consistently executed from here if usetable_hh == 1\n"
  "		: so don't expect the tau values to be tracking along with\n"
  "		: the inf values in hoc\n"
  "\n"
  "	       tinc = -dt * q10\n"
  "        \n"
  "        hyfexp = 1 - exp(tinc/hyftau)\n"
  "	  hysexp = 1 - exp(tinc/hystau)\n"
  "	  hyhtfexp = 1 - exp(tinc/hyhtftau)\n"
  "	  hyhtsexp = 1 - exp(tinc/hyhtstau)\n"
  "}\n"
  " \n"
  "FUNCTION vtrap(x,y) {  :Traps for 0 in denominator of rate eqns.\n"
  "        if (fabs(x/y) < 1e-6) {\n"
  "                vtrap = y*(1 - x/y/2)\n"
  "        }else{  \n"
  "                vtrap = x/(exp(x/y) - 1)\n"
  "        }\n"
  "}\n"
  " \n"
  "UNITSON\n"
  "\n"
  ;
#endif
