/* Created by Language version: 7.7.0 */
/* NOT VECTORIZED */
#define NRN_VECTORIZED 0
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "scoplib_ansi.h"
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
 
#define nrn_init _nrn_init__ichan2
#define _nrn_initial _nrn_initial__ichan2
#define nrn_cur _nrn_cur__ichan2
#define _nrn_current _nrn_current__ichan2
#define nrn_jacob _nrn_jacob__ichan2
#define nrn_state _nrn_state__ichan2
#define _net_receive _net_receive__ichan2 
#define _f_trates _f_trates__ichan2 
#define rates rates__ichan2 
#define states states__ichan2 
#define trates trates__ichan2 
 
#define _threadargscomma_ /**/
#define _threadargsprotocomma_ /**/
#define _threadargs_ /**/
#define _threadargsproto_ /**/
 	/*SUPPRESS 761*/
	/*SUPPRESS 762*/
	/*SUPPRESS 763*/
	/*SUPPRESS 765*/
	 extern double *getarg();
 static double *_p; static Datum *_ppvar;
 
#define t nrn_threads->_t
#define dt nrn_threads->_dt
#define gnatbar _p[0]
#define gkfbar _p[1]
#define gksbar _p[2]
#define gl _p[3]
#define el _p[4]
#define gnat _p[5]
#define gkf _p[6]
#define gks _p[7]
#define inat _p[8]
#define ikf _p[9]
#define iks _p[10]
#define il _p[11]
#define minf _p[12]
#define hinf _p[13]
#define nfinf _p[14]
#define nsinf _p[15]
#define mtau _p[16]
#define htau _p[17]
#define nftau _p[18]
#define nstau _p[19]
#define m _p[20]
#define h _p[21]
#define nf _p[22]
#define ns _p[23]
#define enat _p[24]
#define ekf _p[25]
#define eks _p[26]
#define Dm _p[27]
#define Dh _p[28]
#define Dnf _p[29]
#define Dns _p[30]
#define mexp _p[31]
#define hexp _p[32]
#define nfexp _p[33]
#define nsexp _p[34]
#define _g _p[35]
#define _ion_enat	*_ppvar[0]._pval
#define _ion_inat	*_ppvar[1]._pval
#define _ion_dinatdv	*_ppvar[2]._pval
#define _ion_ekf	*_ppvar[3]._pval
#define _ion_ikf	*_ppvar[4]._pval
#define _ion_dikfdv	*_ppvar[5]._pval
#define _ion_eks	*_ppvar[6]._pval
#define _ion_iks	*_ppvar[7]._pval
#define _ion_diksdv	*_ppvar[8]._pval
 
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
 _p = _prop->param; _ppvar = _prop->dparam;
 }
 static void _hoc_setdata() {
 Prop *_prop, *hoc_getdata_range(int);
 _prop = hoc_getdata_range(_mechtype);
   _setdata(_prop);
 hoc_retpushx(1.);
}
 /* connect user functions to hoc names */
 static VoidFunc hoc_intfunc[] = {
 "setdata_ichan2", _hoc_setdata,
 "rates_ichan2", _hoc_rates,
 "states_ichan2", _hoc_states,
 "trates_ichan2", _hoc_trates,
 "vtrap_ichan2", _hoc_vtrap,
 0, 0
};
#define vtrap vtrap_ichan2
 extern double vtrap( double , double );
 /* declare global and static user variables */
#define usetable usetable_ichan2
 double usetable = 1;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 "usetable_ichan2", 0, 1,
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "gnatbar_ichan2", "mho/cm2",
 "gkfbar_ichan2", "mho/cm2",
 "gksbar_ichan2", "mho/cm2",
 "gl_ichan2", "mho/cm2",
 "el_ichan2", "mV",
 "gnat_ichan2", "mho/cm2",
 "gkf_ichan2", "mho/cm2",
 "gks_ichan2", "mho/cm2",
 "inat_ichan2", "mA/cm2",
 "ikf_ichan2", "mA/cm2",
 "iks_ichan2", "mA/cm2",
 "il_ichan2", "mA/cm2",
 "mtau_ichan2", "ms",
 "htau_ichan2", "ms",
 "nftau_ichan2", "ms",
 "nstau_ichan2", "ms",
 0,0
};
 static double delta_t = 1;
 static double h0 = 0;
 static double m0 = 0;
 static double ns0 = 0;
 static double nf0 = 0;
 static double v = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "usetable_ichan2", &usetable_ichan2,
 0,0
};
 static DoubVec hoc_vdoub[] = {
 0,0,0
};
 static double _sav_indep;
 static void nrn_alloc(Prop*);
static void  nrn_init(_NrnThread*, _Memb_list*, int);
static void nrn_state(_NrnThread*, _Memb_list*, int);
 static void nrn_cur(_NrnThread*, _Memb_list*, int);
static void  nrn_jacob(_NrnThread*, _Memb_list*, int);
 
static int _ode_count(int);
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"ichan2",
 "gnatbar_ichan2",
 "gkfbar_ichan2",
 "gksbar_ichan2",
 "gl_ichan2",
 "el_ichan2",
 0,
 "gnat_ichan2",
 "gkf_ichan2",
 "gks_ichan2",
 "inat_ichan2",
 "ikf_ichan2",
 "iks_ichan2",
 "il_ichan2",
 "minf_ichan2",
 "hinf_ichan2",
 "nfinf_ichan2",
 "nsinf_ichan2",
 "mtau_ichan2",
 "htau_ichan2",
 "nftau_ichan2",
 "nstau_ichan2",
 0,
 "m_ichan2",
 "h_ichan2",
 "nf_ichan2",
 "ns_ichan2",
 0,
 0};
 static Symbol* _nat_sym;
 static Symbol* _kf_sym;
 static Symbol* _ks_sym;
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 36, _prop);
 	/*initialize range parameters*/
 	gnatbar = 0;
 	gkfbar = 0;
 	gksbar = 0;
 	gl = 0;
 	el = 0;
 	_prop->param = _p;
 	_prop->param_size = 36;
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 9, _prop);
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 prop_ion = need_memb(_nat_sym);
 nrn_promote(prop_ion, 0, 1);
 	_ppvar[0]._pval = &prop_ion->param[0]; /* enat */
 	_ppvar[1]._pval = &prop_ion->param[3]; /* inat */
 	_ppvar[2]._pval = &prop_ion->param[4]; /* _ion_dinatdv */
 prop_ion = need_memb(_kf_sym);
 nrn_promote(prop_ion, 0, 1);
 	_ppvar[3]._pval = &prop_ion->param[0]; /* ekf */
 	_ppvar[4]._pval = &prop_ion->param[3]; /* ikf */
 	_ppvar[5]._pval = &prop_ion->param[4]; /* _ion_dikfdv */
 prop_ion = need_memb(_ks_sym);
 nrn_promote(prop_ion, 0, 1);
 	_ppvar[6]._pval = &prop_ion->param[0]; /* eks */
 	_ppvar[7]._pval = &prop_ion->param[3]; /* iks */
 	_ppvar[8]._pval = &prop_ion->param[4]; /* _ion_diksdv */
 
}
 static void _initlists();
 static void _update_ion_pointer(Datum*);
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
extern void _nrn_thread_table_reg(int, void(*)(double*, Datum*, Datum*, _NrnThread*, int));
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 void _ichan2_reg() {
	int _vectorized = 0;
  _initlists();
 	ion_reg("nat", 1.0);
 	ion_reg("kf", 1.0);
 	ion_reg("ks", 1.0);
 	_nat_sym = hoc_lookup("nat_ion");
 	_kf_sym = hoc_lookup("kf_ion");
 	_ks_sym = hoc_lookup("ks_ion");
 	register_mech(_mechanism, nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init, hoc_nrnpointerindex, 0);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
     _nrn_thread_reg(_mechtype, 2, _update_ion_pointer);
 #if NMODL_TEXT
  hoc_reg_nmodl_text(_mechtype, nmodl_file_text);
  hoc_reg_nmodl_filename(_mechtype, nmodl_filename);
#endif
  hoc_register_prop_size(_mechtype, 36, 9);
  hoc_register_dparam_semantics(_mechtype, 0, "nat_ion");
  hoc_register_dparam_semantics(_mechtype, 1, "nat_ion");
  hoc_register_dparam_semantics(_mechtype, 2, "nat_ion");
  hoc_register_dparam_semantics(_mechtype, 3, "kf_ion");
  hoc_register_dparam_semantics(_mechtype, 4, "kf_ion");
  hoc_register_dparam_semantics(_mechtype, 5, "kf_ion");
  hoc_register_dparam_semantics(_mechtype, 6, "ks_ion");
  hoc_register_dparam_semantics(_mechtype, 7, "ks_ion");
  hoc_register_dparam_semantics(_mechtype, 8, "ks_ion");
 	hoc_register_cvode(_mechtype, _ode_count, 0, 0, 0);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 ichan2 /home/danielmk/repos/pyDentate/mechs/ichan2.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
 static double FARADAY = 96520.0;
 static double R = 8.3134;
 static double _zq10 ;
 static double *_t_minf;
 static double *_t_mexp;
 static double *_t_hinf;
 static double *_t_hexp;
 static double *_t_nfinf;
 static double *_t_nfexp;
 static double *_t_nsinf;
 static double *_t_nsexp;
 static double *_t_mtau;
 static double *_t_htau;
 static double *_t_nftau;
 static double *_t_nstau;
static int _reset;
static char *modelname = "ichan2.mod  ";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static int _f_trates(double);
static int rates(double);
static int states();
static int trates(double);
 static void _n_trates(double);
 
static int  states (  ) {
   trates ( _threadargscomma_ v ) ;
   m = m + mexp * ( minf - m ) ;
   h = h + hexp * ( hinf - h ) ;
   nf = nf + nfexp * ( nfinf - nf ) ;
   ns = ns + nsexp * ( nsinf - ns ) ;
   
/*VERBATIM*/
        return 0;
  return 0; }
 
static void _hoc_states(void) {
  double _r;
   _r = 1.;
 states (  );
 hoc_retpushx(_r);
}
 
static int  rates (  double _lv ) {
   double _lalpha , _lbeta , _lsum ;
 _zq10 = pow( 3.0 , ( ( celsius - 6.3 ) / 10.0 ) ) ;
   _lalpha = - 0.3 * vtrap ( _threadargscomma_ ( _lv + 60.0 - 17.0 ) , - 5.0 ) ;
   _lbeta = 0.3 * vtrap ( _threadargscomma_ ( _lv + 60.0 - 45.0 ) , 5.0 ) ;
   _lsum = _lalpha + _lbeta ;
   mtau = 1.0 / _lsum ;
   minf = _lalpha / _lsum ;
   _lalpha = 0.23 / exp ( ( _lv + 60.0 + 5.0 ) / 20.0 ) ;
   _lbeta = 3.33 / ( 1.0 + exp ( ( _lv + 60.0 - 47.5 ) / - 10.0 ) ) ;
   _lsum = _lalpha + _lbeta ;
   htau = 1.0 / _lsum ;
   hinf = _lalpha / _lsum ;
   _lalpha = - 0.028 * vtrap ( _threadargscomma_ ( _lv + 65.0 - 35.0 ) , - 6.0 ) ;
   _lbeta = 0.1056 / exp ( ( _lv + 65.0 - 10.0 ) / 40.0 ) ;
   _lsum = _lalpha + _lbeta ;
   nstau = 1.0 / _lsum ;
   nsinf = _lalpha / _lsum ;
   _lalpha = - 0.07 * vtrap ( _threadargscomma_ ( _lv + 65.0 - 47.0 ) , - 6.0 ) ;
   _lbeta = 0.264 / exp ( ( _lv + 65.0 - 22.0 ) / 40.0 ) ;
   _lsum = _lalpha + _lbeta ;
   nftau = 1.0 / _lsum ;
   nfinf = _lalpha / _lsum ;
    return 0; }
 
static void _hoc_rates(void) {
  double _r;
   _r = 1.;
 rates (  *getarg(1) );
 hoc_retpushx(_r);
}
 static double _mfac_trates, _tmin_trates;
 static void _check_trates();
 static void _check_trates() {
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
    _f_trates(_x);
    _t_minf[_i] = minf;
    _t_mexp[_i] = mexp;
    _t_hinf[_i] = hinf;
    _t_hexp[_i] = hexp;
    _t_nfinf[_i] = nfinf;
    _t_nfexp[_i] = nfexp;
    _t_nsinf[_i] = nsinf;
    _t_nsexp[_i] = nsexp;
    _t_mtau[_i] = mtau;
    _t_htau[_i] = htau;
    _t_nftau[_i] = nftau;
    _t_nstau[_i] = nstau;
   }
   _sav_dt = dt;
   _sav_celsius = celsius;
  }
 }

 static int trates(double _lv){ _check_trates();
 _n_trates(_lv);
 return 0;
 }

 static void _n_trates(double _lv){ int _i, _j;
 double _xi, _theta;
 if (!usetable) {
 _f_trates(_lv); return; 
}
 _xi = _mfac_trates * (_lv - _tmin_trates);
 if (isnan(_xi)) {
  minf = _xi;
  mexp = _xi;
  hinf = _xi;
  hexp = _xi;
  nfinf = _xi;
  nfexp = _xi;
  nsinf = _xi;
  nsexp = _xi;
  mtau = _xi;
  htau = _xi;
  nftau = _xi;
  nstau = _xi;
  return;
 }
 if (_xi <= 0.) {
 minf = _t_minf[0];
 mexp = _t_mexp[0];
 hinf = _t_hinf[0];
 hexp = _t_hexp[0];
 nfinf = _t_nfinf[0];
 nfexp = _t_nfexp[0];
 nsinf = _t_nsinf[0];
 nsexp = _t_nsexp[0];
 mtau = _t_mtau[0];
 htau = _t_htau[0];
 nftau = _t_nftau[0];
 nstau = _t_nstau[0];
 return; }
 if (_xi >= 200.) {
 minf = _t_minf[200];
 mexp = _t_mexp[200];
 hinf = _t_hinf[200];
 hexp = _t_hexp[200];
 nfinf = _t_nfinf[200];
 nfexp = _t_nfexp[200];
 nsinf = _t_nsinf[200];
 nsexp = _t_nsexp[200];
 mtau = _t_mtau[200];
 htau = _t_htau[200];
 nftau = _t_nftau[200];
 nstau = _t_nstau[200];
 return; }
 _i = (int) _xi;
 _theta = _xi - (double)_i;
 minf = _t_minf[_i] + _theta*(_t_minf[_i+1] - _t_minf[_i]);
 mexp = _t_mexp[_i] + _theta*(_t_mexp[_i+1] - _t_mexp[_i]);
 hinf = _t_hinf[_i] + _theta*(_t_hinf[_i+1] - _t_hinf[_i]);
 hexp = _t_hexp[_i] + _theta*(_t_hexp[_i+1] - _t_hexp[_i]);
 nfinf = _t_nfinf[_i] + _theta*(_t_nfinf[_i+1] - _t_nfinf[_i]);
 nfexp = _t_nfexp[_i] + _theta*(_t_nfexp[_i+1] - _t_nfexp[_i]);
 nsinf = _t_nsinf[_i] + _theta*(_t_nsinf[_i+1] - _t_nsinf[_i]);
 nsexp = _t_nsexp[_i] + _theta*(_t_nsexp[_i+1] - _t_nsexp[_i]);
 mtau = _t_mtau[_i] + _theta*(_t_mtau[_i+1] - _t_mtau[_i]);
 htau = _t_htau[_i] + _theta*(_t_htau[_i+1] - _t_htau[_i]);
 nftau = _t_nftau[_i] + _theta*(_t_nftau[_i+1] - _t_nftau[_i]);
 nstau = _t_nstau[_i] + _theta*(_t_nstau[_i+1] - _t_nstau[_i]);
 }

 
static int  _f_trates (  double _lv ) {
   double _ltinc ;
 rates ( _threadargscomma_ _lv ) ;
   _ltinc = - dt * _zq10 ;
   mexp = 1.0 - exp ( _ltinc / mtau ) ;
   hexp = 1.0 - exp ( _ltinc / htau ) ;
   nfexp = 1.0 - exp ( _ltinc / nftau ) ;
   nsexp = 1.0 - exp ( _ltinc / nstau ) ;
    return 0; }
 
static void _hoc_trates(void) {
  double _r;
    _r = 1.;
 trates (  *getarg(1) );
 hoc_retpushx(_r);
}
 
double vtrap (  double _lx , double _ly ) {
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
   _r =  vtrap (  *getarg(1) , *getarg(2) );
 hoc_retpushx(_r);
}
 
static int _ode_count(int _type){ hoc_execerror("ichan2", "cannot be used with CVODE"); return 0;}
 extern void nrn_update_ion_pointer(Symbol*, Datum*, int, int);
 static void _update_ion_pointer(Datum* _ppvar) {
   nrn_update_ion_pointer(_nat_sym, _ppvar, 0, 0);
   nrn_update_ion_pointer(_nat_sym, _ppvar, 1, 3);
   nrn_update_ion_pointer(_nat_sym, _ppvar, 2, 4);
   nrn_update_ion_pointer(_kf_sym, _ppvar, 3, 0);
   nrn_update_ion_pointer(_kf_sym, _ppvar, 4, 3);
   nrn_update_ion_pointer(_kf_sym, _ppvar, 5, 4);
   nrn_update_ion_pointer(_ks_sym, _ppvar, 6, 0);
   nrn_update_ion_pointer(_ks_sym, _ppvar, 7, 3);
   nrn_update_ion_pointer(_ks_sym, _ppvar, 8, 4);
 }

static void initmodel() {
  int _i; double _save;_ninits++;
 _save = t;
 t = 0.0;
{
  h = h0;
  m = m0;
  ns = ns0;
  nf = nf0;
 {
   trates ( _threadargscomma_ v ) ;
   m = minf ;
   h = hinf ;
   nf = nfinf ;
   ns = nsinf ;
   
/*VERBATIM*/
	return 0;
 }
  _sav_indep = t; t = _save;

}
}

static void nrn_init(_NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; double _v; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
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
 v = _v;
  enat = _ion_enat;
  ekf = _ion_ekf;
  eks = _ion_eks;
 initmodel();
   }}

static double _nrn_current(double _v){double _current=0.;v=_v;{ {
   gnat = gnatbar * m * m * m * h ;
   inat = gnat * ( v - enat ) ;
   gkf = gkfbar * nf * nf * nf * nf ;
   ikf = gkf * ( v - ekf ) ;
   gks = gksbar * ns * ns * ns * ns ;
   iks = gks * ( v - eks ) ;
   il = gl * ( v - el ) ;
   }
 _current += inat;
 _current += ikf;
 _current += iks;
 _current += il;

} return _current;
}

static void nrn_cur(_NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; int* _ni; double _rhs, _v; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
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
  enat = _ion_enat;
  ekf = _ion_ekf;
  eks = _ion_eks;
 _g = _nrn_current(_v + .001);
 	{ double _diks;
 double _dikf;
 double _dinat;
  _dinat = inat;
  _dikf = ikf;
  _diks = iks;
 _rhs = _nrn_current(_v);
  _ion_dinatdv += (_dinat - inat)/.001 ;
  _ion_dikfdv += (_dikf - ikf)/.001 ;
  _ion_diksdv += (_diks - iks)/.001 ;
 	}
 _g = (_g - _rhs)/.001;
  _ion_inat += inat ;
  _ion_ikf += ikf ;
  _ion_iks += iks ;
#if CACHEVEC
  if (use_cachevec) {
	VEC_RHS(_ni[_iml]) -= _rhs;
  }else
#endif
  {
	NODERHS(_nd) -= _rhs;
  }
 
}}

static void nrn_jacob(_NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
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
 
}}

static void nrn_state(_NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; double _v = 0.0; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
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
  enat = _ion_enat;
  ekf = _ion_ekf;
  eks = _ion_eks;
 { error =  states();
 if(error){fprintf(stderr,"at line 72 in file ichan2.mod:\n	SOLVE states\n"); nrn_complain(_p); abort_run(error);}
 }   }}

}

static void terminal(){}

static void _initlists() {
 int _i; static int _first = 1;
  if (!_first) return;
   _t_minf = makevector(201*sizeof(double));
   _t_mexp = makevector(201*sizeof(double));
   _t_hinf = makevector(201*sizeof(double));
   _t_hexp = makevector(201*sizeof(double));
   _t_nfinf = makevector(201*sizeof(double));
   _t_nfexp = makevector(201*sizeof(double));
   _t_nsinf = makevector(201*sizeof(double));
   _t_nsexp = makevector(201*sizeof(double));
   _t_mtau = makevector(201*sizeof(double));
   _t_htau = makevector(201*sizeof(double));
   _t_nftau = makevector(201*sizeof(double));
   _t_nstau = makevector(201*sizeof(double));
_first = 0;
}

#if NMODL_TEXT
static const char* nmodl_filename = "/home/danielmk/repos/pyDentate/mechs/ichan2.mod";
static const char* nmodl_file_text = 
  "TITLE ichan2.mod  \n"
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
  "SUFFIX ichan2 \n"
  "USEION nat READ enat WRITE inat VALENCE 1\n"
  "USEION kf READ ekf WRITE ikf  VALENCE 1\n"
  "USEION ks READ eks WRITE iks  VALENCE 1\n"
  "NONSPECIFIC_CURRENT il \n"
  "RANGE  gnat, gkf, gks\n"
  "RANGE gnatbar, gkfbar, gksbar\n"
  "RANGE gl, el\n"
  "RANGE minf, mtau, hinf, htau, nfinf, nftau, inat, ikf, nsinf, nstau, iks\n"
  "}\n"
  " \n"
  "INDEPENDENT {t FROM 0 TO 100 WITH 100 (ms)}\n"
  " \n"
  "PARAMETER {\n"
  "        v (mV) \n"
  "        celsius = 6.3 (degC)\n"
  "        dt (ms) \n"
  "        enat  (mV)\n"
  "	gnatbar (mho/cm2)   \n"
  "        ekf  (mV)\n"
  "	gkfbar (mho/cm2)\n"
  "        eks  (mV)\n"
  "	gksbar (mho/cm2)\n"
  "	gl (mho/cm2)    \n"
  " 	el (mV)\n"
  "}\n"
  " \n"
  "STATE {\n"
  "	m h nf ns\n"
  "}\n"
  " \n"
  "ASSIGNED {\n"
  "         \n"
  "        gnat (mho/cm2) \n"
  "        gkf (mho/cm2)\n"
  "        gks (mho/cm2)\n"
  "\n"
  "        inat (mA/cm2)\n"
  "        ikf (mA/cm2)\n"
  "        iks (mA/cm2)\n"
  "\n"
  "\n"
  "	il (mA/cm2)\n"
  "\n"
  "	minf hinf nfinf nsinf\n"
  " 	mtau (ms) htau (ms) nftau (ms) nstau (ms)\n"
  "	mexp hexp nfexp nsexp\n"
  "} \n"
  "\n"
  "? currents\n"
  "BREAKPOINT {\n"
  "	SOLVE states\n"
  "        gnat = gnatbar*m*m*m*h  \n"
  "        inat = gnat*(v - enat)\n"
  "        gkf = gkfbar*nf*nf*nf*nf\n"
  "        ikf = gkf*(v-ekf)\n"
  "        gks = gksbar*ns*ns*ns*ns\n"
  "        iks = gks*(v-eks)\n"
  "\n"
  "	il = gl*(v-el)\n"
  "}\n"
  " \n"
  "UNITSOFF\n"
  " \n"
  "INITIAL {\n"
  "	trates(v)\n"
  "	\n"
  "	m = minf\n"
  "	h = hinf\n"
  "      nf = nfinf\n"
  "      ns = nsinf\n"
  "	\n"
  "	VERBATIM\n"
  "	return 0;\n"
  "	ENDVERBATIM\n"
  "}\n"
  "\n"
  "? states\n"
  "PROCEDURE states() {	:Computes state variables m, h, and n \n"
  "        trates(v)	:      at the current v and dt.\n"
  "        m = m + mexp*(minf-m)\n"
  "        h = h + hexp*(hinf-h)\n"
  "        nf = nf + nfexp*(nfinf-nf)\n"
  "        ns = ns + nsexp*(nsinf-ns)\n"
  "        VERBATIM\n"
  "        return 0;\n"
  "        ENDVERBATIM\n"
  "}\n"
  " \n"
  "LOCAL q10\n"
  "\n"
  "? rates\n"
  "PROCEDURE rates(v) {  :Computes rate and other constants at current v.\n"
  "                      :Call once from HOC to initialize inf at resting v.\n"
  "        LOCAL  alpha, beta, sum\n"
  "       q10 = 3^((celsius - 6.3)/10)\n"
  "                :\"m\" sodium activation system - act and inact cross at -40\n"
  "	alpha = -0.3*vtrap((v+60-17),-5)\n"
  "	beta = 0.3*vtrap((v+60-45),5)\n"
  "	sum = alpha+beta        \n"
  "	mtau = 1/sum      minf = alpha/sum\n"
  "                :\"h\" sodium inactivation system\n"
  "	alpha = 0.23/exp((v+60+5)/20)\n"
  "	beta = 3.33/(1+exp((v+60-47.5)/-10))\n"
  "	sum = alpha+beta\n"
  "	htau = 1/sum \n"
  "        hinf = alpha/sum \n"
  "             :\"ns\" sKDR activation system\n"
  "        alpha = -0.028*vtrap((v+65-35),-6)\n"
  "	beta = 0.1056/exp((v+65-10)/40)\n"
  "	sum = alpha+beta        \n"
  "	nstau = 1/sum      nsinf = alpha/sum\n"
  "            :\"nf\" fKDR activation system\n"
  "        alpha = -0.07*vtrap((v+65-47),-6)\n"
  "	beta = 0.264/exp((v+65-22)/40)\n"
  "	sum = alpha+beta        \n"
  "	nftau = 1/sum      nfinf = alpha/sum\n"
  "	\n"
  "}\n"
  " \n"
  "PROCEDURE trates(v) {  :Computes rate and other constants at current v.\n"
  "                      :Call once from HOC to initialize inf at resting v.\n"
  "	LOCAL tinc\n"
  "        TABLE minf, mexp, hinf, hexp, nfinf, nfexp, nsinf, nsexp, mtau, htau, nftau, nstau\n"
  "	DEPEND dt, celsius FROM -100 TO 100 WITH 200\n"
  "                           \n"
  "	rates(v)	: not consistently executed from here if usetable_hh == 1\n"
  "		: so don't expect the tau values to be tracking along with\n"
  "		: the inf values in hoc\n"
  "\n"
  "	       tinc = -dt * q10\n"
  "        mexp = 1 - exp(tinc/mtau)\n"
  "        hexp = 1 - exp(tinc/htau)\n"
  "	nfexp = 1 - exp(tinc/nftau)\n"
  "	nsexp = 1 - exp(tinc/nstau)\n"
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
