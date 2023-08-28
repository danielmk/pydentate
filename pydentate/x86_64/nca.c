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
 
#define nrn_init _nrn_init__nca
#define _nrn_initial _nrn_initial__nca
#define nrn_cur _nrn_cur__nca
#define _nrn_current _nrn_current__nca
#define nrn_jacob _nrn_jacob__nca
#define nrn_state _nrn_state__nca
#define _net_receive _net_receive__nca 
#define _f_trates _f_trates__nca 
#define rates rates__nca 
#define states states__nca 
#define trates trates__nca 
 
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
#define gncabar _p[0]
#define gncabar_columnindex 0
#define gnca _p[1]
#define gnca_columnindex 1
#define inca _p[2]
#define inca_columnindex 2
#define cinf _p[3]
#define cinf_columnindex 3
#define dinf _p[4]
#define dinf_columnindex 4
#define ctau _p[5]
#define ctau_columnindex 5
#define dtau _p[6]
#define dtau_columnindex 6
#define c _p[7]
#define c_columnindex 7
#define d _p[8]
#define d_columnindex 8
#define Dc _p[9]
#define Dc_columnindex 9
#define Dd _p[10]
#define Dd_columnindex 10
#define enca _p[11]
#define enca_columnindex 11
#define cexp _p[12]
#define cexp_columnindex 12
#define dexp _p[13]
#define dexp_columnindex 13
#define v _p[14]
#define v_columnindex 14
#define _g _p[15]
#define _g_columnindex 15
#define _ion_enca	*_ppvar[0]._pval
#define _ion_inca	*_ppvar[1]._pval
#define _ion_dincadv	*_ppvar[2]._pval
 
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
 "setdata_nca", _hoc_setdata,
 "rates_nca", _hoc_rates,
 "states_nca", _hoc_states,
 "trates_nca", _hoc_trates,
 "vtrap_nca", _hoc_vtrap,
 0, 0
};
#define vtrap vtrap_nca
 extern double vtrap( _threadargsprotocomma_ double , double );
 
static void _check_trates(double*, Datum*, Datum*, NrnThread*); 
static void _check_table_thread(double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt, int _type) {
   _check_trates(_p, _ppvar, _thread, _nt);
 }
 #define _zq10 _thread[0]._pval[0]
 /* declare global and static user variables */
#define usetable usetable_nca
 double usetable = 1;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 "usetable_nca", 0, 1,
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "gncabar_nca", "mho/cm2",
 "gnca_nca", "mho/cm2",
 "inca_nca", "mA/cm2",
 "ctau_nca", "ms",
 "dtau_nca", "ms",
 0,0
};
 static double c0 = 0;
 static double delta_t = 1;
 static double d0 = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "usetable_nca", &usetable_nca,
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
"nca",
 "gncabar_nca",
 0,
 "gnca_nca",
 "inca_nca",
 "cinf_nca",
 "dinf_nca",
 "ctau_nca",
 "dtau_nca",
 0,
 "c_nca",
 "d_nca",
 0,
 0};
 static Symbol* _nca_sym;
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 16, _prop);
 	/*initialize range parameters*/
 	gncabar = 0;
 	_prop->param = _p;
 	_prop->param_size = 16;
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 3, _prop);
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 prop_ion = need_memb(_nca_sym);
 nrn_promote(prop_ion, 0, 1);
 	_ppvar[0]._pval = &prop_ion->param[0]; /* enca */
 	_ppvar[1]._pval = &prop_ion->param[3]; /* inca */
 	_ppvar[2]._pval = &prop_ion->param[4]; /* _ion_dincadv */
 
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

 void _nca_reg() {
	int _vectorized = 1;
  _initlists();
 	ion_reg("nca", 2.0);
 	_nca_sym = hoc_lookup("nca_ion");
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
  hoc_register_prop_size(_mechtype, 16, 3);
  hoc_register_dparam_semantics(_mechtype, 0, "nca_ion");
  hoc_register_dparam_semantics(_mechtype, 1, "nca_ion");
  hoc_register_dparam_semantics(_mechtype, 2, "nca_ion");
 	hoc_register_cvode(_mechtype, _ode_count, 0, 0, 0);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 nca /Users/temma/ghq/pydentate/mechs/nca.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
 static double FARADAY = 96520.0;
 static double R = 8.3134;
 /*Top LOCAL _zq10 */
 static double *_t_cinf;
 static double *_t_cexp;
 static double *_t_dinf;
 static double *_t_dexp;
 static double *_t_ctau;
 static double *_t_dtau;
static int _reset;
static char *modelname = "nca.mod  ";

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
   c = c + cexp * ( cinf - c ) ;
   d = d + dexp * ( dinf - d ) ;
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
   _lalpha = - 0.19 * vtrap ( _threadargscomma_ _lv - 19.88 , - 10.0 ) ;
   _lbeta = 0.046 * exp ( - _lv / 20.73 ) ;
   _lsum = _lalpha + _lbeta ;
   ctau = 1.0 / _lsum ;
   cinf = _lalpha / _lsum ;
   _lalpha = 0.00016 / exp ( - _lv / 48.4 ) ;
   _lbeta = 1.0 / ( exp ( ( - _lv + 39.0 ) / 10.0 ) + 1.0 ) ;
   _lsum = _lalpha + _lbeta ;
   dtau = 1.0 / _lsum ;
   dinf = _lalpha / _lsum ;
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
   _tmin_trates =  - 100.0 ;
   _tmax =  100.0 ;
   _dx = (_tmax - _tmin_trates)/200.; _mfac_trates = 1./_dx;
   for (_i=0, _x=_tmin_trates; _i < 201; _x += _dx, _i++) {
    _f_trates(_p, _ppvar, _thread, _nt, _x);
    _t_cinf[_i] = cinf;
    _t_cexp[_i] = cexp;
    _t_dinf[_i] = dinf;
    _t_dexp[_i] = dexp;
    _t_ctau[_i] = ctau;
    _t_dtau[_i] = dtau;
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
  cinf = _xi;
  cexp = _xi;
  dinf = _xi;
  dexp = _xi;
  ctau = _xi;
  dtau = _xi;
  return;
 }
 if (_xi <= 0.) {
 cinf = _t_cinf[0];
 cexp = _t_cexp[0];
 dinf = _t_dinf[0];
 dexp = _t_dexp[0];
 ctau = _t_ctau[0];
 dtau = _t_dtau[0];
 return; }
 if (_xi >= 200.) {
 cinf = _t_cinf[200];
 cexp = _t_cexp[200];
 dinf = _t_dinf[200];
 dexp = _t_dexp[200];
 ctau = _t_ctau[200];
 dtau = _t_dtau[200];
 return; }
 _i = (int) _xi;
 _theta = _xi - (double)_i;
 cinf = _t_cinf[_i] + _theta*(_t_cinf[_i+1] - _t_cinf[_i]);
 cexp = _t_cexp[_i] + _theta*(_t_cexp[_i+1] - _t_cexp[_i]);
 dinf = _t_dinf[_i] + _theta*(_t_dinf[_i+1] - _t_dinf[_i]);
 dexp = _t_dexp[_i] + _theta*(_t_dexp[_i+1] - _t_dexp[_i]);
 ctau = _t_ctau[_i] + _theta*(_t_ctau[_i+1] - _t_ctau[_i]);
 dtau = _t_dtau[_i] + _theta*(_t_dtau[_i+1] - _t_dtau[_i]);
 }

 
static int  _f_trates ( _threadargsprotocomma_ double _lv ) {
   double _ltinc ;
 rates ( _threadargscomma_ _lv ) ;
   _ltinc = - dt * _zq10 ;
   cexp = 1.0 - exp ( _ltinc / ctau ) ;
   dexp = 1.0 - exp ( _ltinc / dtau ) ;
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
 
static int _ode_count(int _type){ hoc_execerror("nca", "cannot be used with CVODE"); return 0;}
 
static void _thread_mem_init(Datum* _thread) {
   _thread[0]._pval = (double*)ecalloc(1, sizeof(double));
 }
 
static void _thread_cleanup(Datum* _thread) {
   free((void*)(_thread[0]._pval));
 }
 extern void nrn_update_ion_pointer(Symbol*, Datum*, int, int);
 static void _update_ion_pointer(Datum* _ppvar) {
   nrn_update_ion_pointer(_nca_sym, _ppvar, 0, 0);
   nrn_update_ion_pointer(_nca_sym, _ppvar, 1, 3);
   nrn_update_ion_pointer(_nca_sym, _ppvar, 2, 4);
 }

static void initmodel(double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt) {
  int _i; double _save;{
  c = c0;
  d = d0;
 {
   trates ( _threadargscomma_ v ) ;
   c = cinf ;
   d = dinf ;
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
  enca = _ion_enca;
 initmodel(_p, _ppvar, _thread, _nt);
 }
}

static double _nrn_current(double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt, double _v){double _current=0.;v=_v;{ {
   gnca = gncabar * c * c * d ;
   inca = gnca * ( v - enca ) ;
   }
 _current += inca;

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
  enca = _ion_enca;
 _g = _nrn_current(_p, _ppvar, _thread, _nt, _v + .001);
 	{ double _dinca;
  _dinca = inca;
 _rhs = _nrn_current(_p, _ppvar, _thread, _nt, _v);
  _ion_dincadv += (_dinca - inca)/.001 ;
 	}
 _g = (_g - _rhs)/.001;
  _ion_inca += inca ;
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
  enca = _ion_enca;
 {  { states(_p, _ppvar, _thread, _nt); }
  } }}

}

static void terminal(){}

static void _initlists(){
 double _x; double* _p = &_x;
 int _i; static int _first = 1;
  if (!_first) return;
   _t_cinf = makevector(201*sizeof(double));
   _t_cexp = makevector(201*sizeof(double));
   _t_dinf = makevector(201*sizeof(double));
   _t_dexp = makevector(201*sizeof(double));
   _t_ctau = makevector(201*sizeof(double));
   _t_dtau = makevector(201*sizeof(double));
_first = 0;
}

#if defined(__cplusplus)
} /* extern "C" */
#endif

#if NMODL_TEXT
static const char* nmodl_filename = "/Users/temma/ghq/pydentate/mechs/nca.mod";
static const char* nmodl_file_text = 
  "TITLE nca.mod  \n"
  " \n"
  "COMMENT\n"
  "konduktivitas valtozas hatasa- somaban \n"
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
  "SUFFIX nca\n"
  "USEION nca READ enca WRITE inca VALENCE 2 \n"
  "RANGE  gnca\n"
  "RANGE gncabar\n"
  "RANGE cinf, ctau, dinf, dtau, inca\n"
  "}\n"
  " \n"
  "INDEPENDENT {t FROM 0 TO 100 WITH 100 (ms)}\n"
  " \n"
  "PARAMETER {\n"
  "        v (mV) \n"
  "        celsius = 6.3 (degC)\n"
  "        dt (ms) \n"
  "	gncabar (mho/cm2)\n"
  "}\n"
  " \n"
  "STATE {\n"
  "	c d\n"
  "}\n"
  " \n"
  "ASSIGNED {\n"
  "	  gnca (mho/cm2)\n"
  "	inca (mA/cm2)\n"
  "	enca (mV)\n"
  "\n"
  "	cinf dinf\n"
  "	ctau (ms) dtau (ms) \n"
  "	cexp dexp      \n"
  "} \n"
  "\n"
  "? currents\n"
  "BREAKPOINT {\n"
  "	SOLVE states\n"
  "        gnca = gncabar*c*c*d\n"
  "	inca = gnca*(v-enca)\n"
  "}\n"
  " \n"
  "UNITSOFF\n"
  " \n"
  "INITIAL {\n"
  "	trates(v)\n"
  "	c = cinf\n"
  "	d = dinf\n"
  "}\n"
  "\n"
  "? states\n"
  "PROCEDURE states() {	:Computes state variables m, h, and n \n"
  "        trates(v)	:      at the current v and dt.\n"
  "	c = c + cexp*(cinf-c)\n"
  "	d = d + dexp*(dinf-d)\n"
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
  "                :\"c\" NCa activation system\n"
  "        alpha = -0.19*vtrap(v-19.88,-10)\n"
  "	beta = 0.046*exp(-v/20.73)\n"
  "	sum = alpha+beta        \n"
  "	ctau = 1/sum      cinf = alpha/sum\n"
  "                :\"d\" NCa inactivation system\n"
  "	alpha = 0.00016/exp(-v/48.4)\n"
  "	beta = 1/(exp((-v+39)/10)+1)\n"
  "	sum = alpha+beta        \n"
  "	dtau = 1/sum      dinf = alpha/sum\n"
  "}\n"
  " \n"
  "PROCEDURE trates(v) {  :Computes rate and other constants at current v.\n"
  "                      :Call once from HOC to initialize inf at resting v.\n"
  "	LOCAL tinc\n"
  "        TABLE  cinf, cexp, dinf, dexp, ctau, dtau\n"
  "	DEPEND dt, celsius FROM -100 TO 100 WITH 200\n"
  "                           \n"
  "	rates(v)	: not consistently executed from here if usetable_hh == 1\n"
  "		: so don't expect the tau values to be tracking along with\n"
  "		: the inf values in hoc\n"
  "\n"
  "	       tinc = -dt * q10\n"
  "	cexp = 1 - exp(tinc/ctau)\n"
  "	dexp = 1 - exp(tinc/dtau)\n"
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
