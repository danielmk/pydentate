/* Created by Language version: 7.5.0 */
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
 
#define nrn_init _nrn_init__cat
#define _nrn_initial _nrn_initial__cat
#define nrn_cur _nrn_cur__cat
#define _nrn_current _nrn_current__cat
#define nrn_jacob _nrn_jacob__cat
#define nrn_state _nrn_state__cat
#define _net_receive _net_receive__cat 
#define states states__cat 
 
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
#define gcatbar _p[0]
#define itca _p[1]
#define m _p[2]
#define h _p[3]
#define cai _p[4]
#define cao _p[5]
#define Dm _p[6]
#define Dh _p[7]
#define gcat _p[8]
#define etca _p[9]
#define _g _p[10]
#define _ion_etca	*_ppvar[0]._pval
#define _ion_itca	*_ppvar[1]._pval
#define _ion_ditcadv	*_ppvar[2]._pval
#define _ion_cai	*_ppvar[3]._pval
#define _ion_cao	*_ppvar[4]._pval
 
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
 static void _hoc_KTF(void);
 static void _hoc_efun(void);
 static void _hoc_ghk(void);
 static void _hoc_h_tau(void);
 static void _hoc_hinf(void);
 static void _hoc_m_tau(void);
 static void _hoc_minf(void);
 static int _mechtype;
extern void _nrn_cacheloop_reg(int, int);
extern void hoc_register_prop_size(int, int, int);
extern void hoc_register_limits(int, HocParmLimits*);
extern void hoc_register_units(int, HocParmUnits*);
extern void nrn_promote(Prop*, int, int);
extern Memb_func* memb_func;
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
 "setdata_cat", _hoc_setdata,
 "KTF_cat", _hoc_KTF,
 "efun_cat", _hoc_efun,
 "ghk_cat", _hoc_ghk,
 "h_tau_cat", _hoc_h_tau,
 "hinf_cat", _hoc_hinf,
 "m_tau_cat", _hoc_m_tau,
 "minf_cat", _hoc_minf,
 0, 0
};
#define KTF KTF_cat
#define _f_h_tau _f_h_tau_cat
#define _f_m_tau _f_m_tau_cat
#define _f_minf _f_minf_cat
#define _f_hinf _f_hinf_cat
#define efun efun_cat
#define ghk ghk_cat
#define h_tau h_tau_cat
#define hinf hinf_cat
#define m_tau m_tau_cat
#define minf minf_cat
 extern double KTF( double );
 extern double _f_h_tau( double );
 extern double _f_m_tau( double );
 extern double _f_minf( double );
 extern double _f_hinf( double );
 extern double efun( double );
 extern double ghk( double , double , double );
 extern double h_tau( double );
 extern double hinf( double );
 extern double m_tau( double );
 extern double minf( double );
 /* declare global and static user variables */
#define usetable usetable_cat
 double usetable = 1;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 "usetable_cat", 0, 1,
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "gcatbar_cat", "mho/cm2",
 "itca_cat", "mA/cm2",
 0,0
};
 static double delta_t = 0.01;
 static double h0 = 0;
 static double m0 = 0;
 static double v = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "usetable_cat", &usetable_cat,
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
static void _ode_map(int, double**, double**, double*, Datum*, double*, int);
static void _ode_spec(_NrnThread*, _Memb_list*, int);
static void _ode_matsol(_NrnThread*, _Memb_list*, int);
 
#define _cvode_ieq _ppvar[5]._i
 static void _ode_matsol_instance1(_threadargsproto_);
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.5.0",
"cat",
 "gcatbar_cat",
 0,
 "itca_cat",
 0,
 "m_cat",
 "h_cat",
 0,
 0};
 static Symbol* _tca_sym;
 static Symbol* _ca_sym;
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 11, _prop);
 	/*initialize range parameters*/
 	gcatbar = 0.003;
 	_prop->param = _p;
 	_prop->param_size = 11;
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 6, _prop);
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 prop_ion = need_memb(_tca_sym);
 nrn_promote(prop_ion, 0, 1);
 	_ppvar[0]._pval = &prop_ion->param[0]; /* etca */
 	_ppvar[1]._pval = &prop_ion->param[3]; /* itca */
 	_ppvar[2]._pval = &prop_ion->param[4]; /* _ion_ditcadv */
 prop_ion = need_memb(_ca_sym);
 nrn_promote(prop_ion, 1, 0);
 	_ppvar[3]._pval = &prop_ion->param[1]; /* cai */
 	_ppvar[4]._pval = &prop_ion->param[2]; /* cao */
 
}
 static void _initlists();
  /* some states have an absolute tolerance */
 static Symbol** _atollist;
 static HocStateTolerance _hoc_state_tol[] = {
 0,0
};
 static void _update_ion_pointer(Datum*);
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
extern void _nrn_thread_table_reg(int, void(*)(double*, Datum*, Datum*, _NrnThread*, int));
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 void _tca_reg() {
	int _vectorized = 0;
  _initlists();
 	ion_reg("tca", 2.0);
 	ion_reg("ca", 2.0);
 	_tca_sym = hoc_lookup("tca_ion");
 	_ca_sym = hoc_lookup("ca_ion");
 	register_mech(_mechanism, nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init, hoc_nrnpointerindex, 0);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
     _nrn_thread_reg(_mechtype, 2, _update_ion_pointer);
  hoc_register_prop_size(_mechtype, 11, 6);
  hoc_register_dparam_semantics(_mechtype, 0, "tca_ion");
  hoc_register_dparam_semantics(_mechtype, 1, "tca_ion");
  hoc_register_dparam_semantics(_mechtype, 2, "tca_ion");
  hoc_register_dparam_semantics(_mechtype, 3, "ca_ion");
  hoc_register_dparam_semantics(_mechtype, 4, "ca_ion");
  hoc_register_dparam_semantics(_mechtype, 5, "cvodeieq");
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 cat /home/daniel/repos/pyDentate/mechs_7-6/x86_64/tca.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
 static double FARADAY = 96520.0;
 static double R = 8.3134;
 static double KTOMV = .0853;
 static double *_t_hinf;
 static double *_t_minf;
 static double *_t_m_tau;
 static double *_t_h_tau;
static int _reset;
static char *modelname = "T-calcium channel From Migliore CA3";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
 
static int _ode_spec1(_threadargsproto_);
/*static int _ode_matsol1(_threadargsproto_);*/
 static double _n_h_tau(double);
 static double _n_m_tau(double);
 static double _n_minf(double);
 static double _n_hinf(double);
 static int _slist1[2], _dlist1[2];
 static int states(_threadargsproto_);
 
/*CVODE*/
 static int _ode_spec1 () {_reset=0;
 {
   Dm = ( minf ( _threadargscomma_ v ) - m ) / m_tau ( _threadargscomma_ v ) ;
   Dh = ( hinf ( _threadargscomma_ v ) - h ) / h_tau ( _threadargscomma_ v ) ;
   }
 return _reset;
}
 static int _ode_matsol1 () {
 Dm = Dm  / (1. - dt*( ( ( ( - 1.0 ) ) ) / m_tau ( _threadargscomma_ v ) )) ;
 Dh = Dh  / (1. - dt*( ( ( ( - 1.0 ) ) ) / h_tau ( _threadargscomma_ v ) )) ;
  return 0;
}
 /*END CVODE*/
 static int states () {_reset=0;
 {
    m = m + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / m_tau ( _threadargscomma_ v ))))*(- ( ( ( minf ( _threadargscomma_ v ) ) ) / m_tau ( _threadargscomma_ v ) ) / ( ( ( ( - 1.0 ) ) ) / m_tau ( _threadargscomma_ v ) ) - m) ;
    h = h + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / h_tau ( _threadargscomma_ v ))))*(- ( ( ( hinf ( _threadargscomma_ v ) ) ) / h_tau ( _threadargscomma_ v ) ) / ( ( ( ( - 1.0 ) ) ) / h_tau ( _threadargscomma_ v ) ) - h) ;
   }
  return 0;
}
 
double ghk (  double _lv , double _lci , double _lco ) {
   double _lghk;
 double _lnu , _lf ;
 _lf = KTF ( _threadargscomma_ celsius ) / 2.0 ;
   _lnu = _lv / _lf ;
   _lghk = - _lf * ( 1. - ( _lci / _lco ) * exp ( _lnu ) ) * efun ( _threadargscomma_ _lnu ) ;
   
return _lghk;
 }
 
static void _hoc_ghk(void) {
  double _r;
   _r =  ghk (  *getarg(1) , *getarg(2) , *getarg(3) );
 hoc_retpushx(_r);
}
 
double KTF (  double _lcelsius ) {
   double _lKTF;
 _lKTF = ( ( 25. / 293.15 ) * ( _lcelsius + 273.15 ) ) ;
   
return _lKTF;
 }
 
static void _hoc_KTF(void) {
  double _r;
   _r =  KTF (  *getarg(1) );
 hoc_retpushx(_r);
}
 
double efun (  double _lz ) {
   double _lefun;
 if ( fabs ( _lz ) < 1e-4 ) {
     _lefun = 1.0 - _lz / 2.0 ;
     }
   else {
     _lefun = _lz / ( exp ( _lz ) - 1.0 ) ;
     }
   
return _lefun;
 }
 
static void _hoc_efun(void) {
  double _r;
   _r =  efun (  *getarg(1) );
 hoc_retpushx(_r);
}
 static double _mfac_hinf, _tmin_hinf;
 static void _check_hinf();
 static void _check_hinf() {
  static int _maktable=1; int _i, _j, _ix = 0;
  double _xi, _tmax;
  if (!usetable) {return;}
  if (_maktable) { double _x, _dx; _maktable=0;
   _tmin_hinf =  - 150.0 ;
   _tmax =  150.0 ;
   _dx = (_tmax - _tmin_hinf)/200.; _mfac_hinf = 1./_dx;
   for (_i=0, _x=_tmin_hinf; _i < 201; _x += _dx, _i++) {
    _t_hinf[_i] = _f_hinf(_x);
   }
  }
 }

 double hinf(double _lv){ _check_hinf();
 return _n_hinf(_lv);
 }

 static double _n_hinf(double _lv){ int _i, _j;
 double _xi, _theta;
 if (!usetable) {
 return _f_hinf(_lv); 
}
 _xi = _mfac_hinf * (_lv - _tmin_hinf);
 if (isnan(_xi)) {
  return _xi; }
 if (_xi <= 0.) {
 return _t_hinf[0];
 }
 if (_xi >= 200.) {
 return _t_hinf[200];
 }
 _i = (int) _xi;
 return _t_hinf[_i] + (_xi - (double)_i)*(_t_hinf[_i+1] - _t_hinf[_i]);
 }

 
double _f_hinf (  double _lv ) {
   double _lhinf;
 double _la , _lb ;
 _la = 1.e-6 * exp ( - _lv / 16.26 ) ;
   _lb = 1.0 / ( exp ( ( - _lv + 29.79 ) / 10.0 ) + 1.0 ) ;
   _lhinf = _la / ( _la + _lb ) ;
   
return _lhinf;
 }
 
static void _hoc_hinf(void) {
  double _r;
    _r =  hinf (  *getarg(1) );
 hoc_retpushx(_r);
}
 static double _mfac_minf, _tmin_minf;
 static void _check_minf();
 static void _check_minf() {
  static int _maktable=1; int _i, _j, _ix = 0;
  double _xi, _tmax;
  if (!usetable) {return;}
  if (_maktable) { double _x, _dx; _maktable=0;
   _tmin_minf =  - 150.0 ;
   _tmax =  150.0 ;
   _dx = (_tmax - _tmin_minf)/200.; _mfac_minf = 1./_dx;
   for (_i=0, _x=_tmin_minf; _i < 201; _x += _dx, _i++) {
    _t_minf[_i] = _f_minf(_x);
   }
  }
 }

 double minf(double _lv){ _check_minf();
 return _n_minf(_lv);
 }

 static double _n_minf(double _lv){ int _i, _j;
 double _xi, _theta;
 if (!usetable) {
 return _f_minf(_lv); 
}
 _xi = _mfac_minf * (_lv - _tmin_minf);
 if (isnan(_xi)) {
  return _xi; }
 if (_xi <= 0.) {
 return _t_minf[0];
 }
 if (_xi >= 200.) {
 return _t_minf[200];
 }
 _i = (int) _xi;
 return _t_minf[_i] + (_xi - (double)_i)*(_t_minf[_i+1] - _t_minf[_i]);
 }

 
double _f_minf (  double _lv ) {
   double _lminf;
 double _la , _lb ;
 _la = 0.2 * ( - 1.0 * _lv + 19.26 ) / ( exp ( ( - 1.0 * _lv + 19.26 ) / 10.0 ) - 1.0 ) ;
   _lb = 0.009 * exp ( - _lv / 22.03 ) ;
   _lminf = _la / ( _la + _lb ) ;
   
return _lminf;
 }
 
static void _hoc_minf(void) {
  double _r;
    _r =  minf (  *getarg(1) );
 hoc_retpushx(_r);
}
 static double _mfac_m_tau, _tmin_m_tau;
 static void _check_m_tau();
 static void _check_m_tau() {
  static int _maktable=1; int _i, _j, _ix = 0;
  double _xi, _tmax;
  if (!usetable) {return;}
  if (_maktable) { double _x, _dx; _maktable=0;
   _tmin_m_tau =  - 150.0 ;
   _tmax =  150.0 ;
   _dx = (_tmax - _tmin_m_tau)/200.; _mfac_m_tau = 1./_dx;
   for (_i=0, _x=_tmin_m_tau; _i < 201; _x += _dx, _i++) {
    _t_m_tau[_i] = _f_m_tau(_x);
   }
  }
 }

 double m_tau(double _lv){ _check_m_tau();
 return _n_m_tau(_lv);
 }

 static double _n_m_tau(double _lv){ int _i, _j;
 double _xi, _theta;
 if (!usetable) {
 return _f_m_tau(_lv); 
}
 _xi = _mfac_m_tau * (_lv - _tmin_m_tau);
 if (isnan(_xi)) {
  return _xi; }
 if (_xi <= 0.) {
 return _t_m_tau[0];
 }
 if (_xi >= 200.) {
 return _t_m_tau[200];
 }
 _i = (int) _xi;
 return _t_m_tau[_i] + (_xi - (double)_i)*(_t_m_tau[_i+1] - _t_m_tau[_i]);
 }

 
double _f_m_tau (  double _lv ) {
   double _lm_tau;
 double _la , _lb ;
 _la = 0.2 * ( - 1.0 * _lv + 19.26 ) / ( exp ( ( - 1.0 * _lv + 19.26 ) / 10.0 ) - 1.0 ) ;
   _lb = 0.009 * exp ( - _lv / 22.03 ) ;
   _lm_tau = 1.0 / ( _la + _lb ) ;
   
return _lm_tau;
 }
 
static void _hoc_m_tau(void) {
  double _r;
    _r =  m_tau (  *getarg(1) );
 hoc_retpushx(_r);
}
 static double _mfac_h_tau, _tmin_h_tau;
 static void _check_h_tau();
 static void _check_h_tau() {
  static int _maktable=1; int _i, _j, _ix = 0;
  double _xi, _tmax;
  if (!usetable) {return;}
  if (_maktable) { double _x, _dx; _maktable=0;
   _tmin_h_tau =  - 150.0 ;
   _tmax =  150.0 ;
   _dx = (_tmax - _tmin_h_tau)/200.; _mfac_h_tau = 1./_dx;
   for (_i=0, _x=_tmin_h_tau; _i < 201; _x += _dx, _i++) {
    _t_h_tau[_i] = _f_h_tau(_x);
   }
  }
 }

 double h_tau(double _lv){ _check_h_tau();
 return _n_h_tau(_lv);
 }

 static double _n_h_tau(double _lv){ int _i, _j;
 double _xi, _theta;
 if (!usetable) {
 return _f_h_tau(_lv); 
}
 _xi = _mfac_h_tau * (_lv - _tmin_h_tau);
 if (isnan(_xi)) {
  return _xi; }
 if (_xi <= 0.) {
 return _t_h_tau[0];
 }
 if (_xi >= 200.) {
 return _t_h_tau[200];
 }
 _i = (int) _xi;
 return _t_h_tau[_i] + (_xi - (double)_i)*(_t_h_tau[_i+1] - _t_h_tau[_i]);
 }

 
double _f_h_tau (  double _lv ) {
   double _lh_tau;
 double _la , _lb ;
 _la = 1.e-6 * exp ( - _lv / 16.26 ) ;
   _lb = 1.0 / ( exp ( ( - _lv + 29.79 ) / 10. ) + 1. ) ;
   _lh_tau = 1.0 / ( _la + _lb ) ;
   
return _lh_tau;
 }
 
static void _hoc_h_tau(void) {
  double _r;
    _r =  h_tau (  *getarg(1) );
 hoc_retpushx(_r);
}
 
static int _ode_count(int _type){ return 2;}
 
static void _ode_spec(_NrnThread* _nt, _Memb_list* _ml, int _type) {
   Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
  etca = _ion_etca;
  cai = _ion_cai;
  cao = _ion_cao;
     _ode_spec1 ();
  }}
 
static void _ode_map(int _ieq, double** _pv, double** _pvdot, double* _pp, Datum* _ppd, double* _atol, int _type) { 
 	int _i; _p = _pp; _ppvar = _ppd;
	_cvode_ieq = _ieq;
	for (_i=0; _i < 2; ++_i) {
		_pv[_i] = _pp + _slist1[_i];  _pvdot[_i] = _pp + _dlist1[_i];
		_cvode_abstol(_atollist, _atol, _i);
	}
 }
 
static void _ode_matsol_instance1(_threadargsproto_) {
 _ode_matsol1 ();
 }
 
static void _ode_matsol(_NrnThread* _nt, _Memb_list* _ml, int _type) {
   Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
  etca = _ion_etca;
  cai = _ion_cai;
  cao = _ion_cao;
 _ode_matsol_instance1(_threadargs_);
 }}
 extern void nrn_update_ion_pointer(Symbol*, Datum*, int, int);
 static void _update_ion_pointer(Datum* _ppvar) {
   nrn_update_ion_pointer(_tca_sym, _ppvar, 0, 0);
   nrn_update_ion_pointer(_tca_sym, _ppvar, 1, 3);
   nrn_update_ion_pointer(_tca_sym, _ppvar, 2, 4);
   nrn_update_ion_pointer(_ca_sym, _ppvar, 3, 1);
   nrn_update_ion_pointer(_ca_sym, _ppvar, 4, 2);
 }

static void initmodel() {
  int _i; double _save;_ninits++;
 _save = t;
 t = 0.0;
{
  h = h0;
  m = m0;
 {
   m = minf ( _threadargscomma_ v ) ;
   h = hinf ( _threadargscomma_ v ) ;
   
/*VERBATIM*/
	cai=_ion_cai;
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
  etca = _ion_etca;
  cai = _ion_cai;
  cao = _ion_cao;
 initmodel();
 }}

static double _nrn_current(double _v){double _current=0.;v=_v;{ {
   gcat = gcatbar * m * m * h ;
   itca = gcat * ghk ( _threadargscomma_ v , cai , cao ) ;
   }
 _current += itca;

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
  etca = _ion_etca;
  cai = _ion_cai;
  cao = _ion_cao;
 _g = _nrn_current(_v + .001);
 	{ double _ditca;
  _ditca = itca;
 _rhs = _nrn_current(_v);
  _ion_ditcadv += (_ditca - itca)/.001 ;
 	}
 _g = (_g - _rhs)/.001;
  _ion_itca += itca ;
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
  etca = _ion_etca;
  cai = _ion_cai;
  cao = _ion_cao;
 { error =  states();
 if(error){fprintf(stderr,"at line 49 in file tca.mod:\n	SOLVE states METHOD cnexp\n"); nrn_complain(_p); abort_run(error);}
 } }}

}

static void terminal(){}

static void _initlists() {
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = &(m) - _p;  _dlist1[0] = &(Dm) - _p;
 _slist1[1] = &(h) - _p;  _dlist1[1] = &(Dh) - _p;
   _t_hinf = makevector(201*sizeof(double));
   _t_minf = makevector(201*sizeof(double));
   _t_m_tau = makevector(201*sizeof(double));
   _t_h_tau = makevector(201*sizeof(double));
_first = 0;
}
