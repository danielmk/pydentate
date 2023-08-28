/* Created by Language version: 7.7.0 */
/* NOT VECTORIZED */
#define NRN_VECTORIZED 0
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
 
#define nrn_init _nrn_init__cagk
#define _nrn_initial _nrn_initial__cagk
#define nrn_cur _nrn_cur__cagk
#define _nrn_current _nrn_current__cagk
#define nrn_jacob _nrn_jacob__cagk
#define nrn_state _nrn_state__cagk
#define _net_receive _net_receive__cagk 
#define rate rate__cagk 
#define state state__cagk 
 
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
#define gkbar _p[0]
#define gkbar_columnindex 0
#define ik _p[1]
#define ik_columnindex 1
#define gkca _p[2]
#define gkca_columnindex 2
#define o _p[3]
#define o_columnindex 3
#define ek _p[4]
#define ek_columnindex 4
#define lcai _p[5]
#define lcai_columnindex 5
#define ncai _p[6]
#define ncai_columnindex 6
#define tcai _p[7]
#define tcai_columnindex 7
#define Do _p[8]
#define Do_columnindex 8
#define _g _p[9]
#define _g_columnindex 9
#define _ion_ncai	*_ppvar[0]._pval
#define _ion_lcai	*_ppvar[1]._pval
#define _ion_tcai	*_ppvar[2]._pval
#define _ion_ek	*_ppvar[3]._pval
#define _ion_ik	*_ppvar[4]._pval
#define _ion_dikdv	*_ppvar[5]._pval
 
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
 static void _hoc_alp(void);
 static void _hoc_bet(void);
 static void _hoc_exp1(void);
 static void _hoc_rate(void);
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
 "setdata_cagk", _hoc_setdata,
 "alp_cagk", _hoc_alp,
 "bet_cagk", _hoc_bet,
 "exp1_cagk", _hoc_exp1,
 "rate_cagk", _hoc_rate,
 0, 0
};
#define alp alp_cagk
#define bet bet_cagk
#define exp1 exp1_cagk
 extern double alp( double , double );
 extern double bet( double , double );
 extern double exp1( double , double , double );
 /* declare global and static user variables */
#define abar abar_cagk
 double abar = 0.28;
#define bbar bbar_cagk
 double bbar = 0.48;
#define cai cai_cagk
 double cai = 5e-05;
#define d2 d2_cagk
 double d2 = 1;
#define d1 d1_cagk
 double d1 = 0.84;
#define k2 k2_cagk
 double k2 = 1.3e-07;
#define k1 k1_cagk
 double k1 = 0.00048;
#define otau otau_cagk
 double otau = 0;
#define oinf oinf_cagk
 double oinf = 0;
#define st st_cagk
 double st = 1;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "cai_cagk", "mM",
 "k1_cagk", "mM",
 "k2_cagk", "mM",
 "abar_cagk", "/ms",
 "bbar_cagk", "/ms",
 "st_cagk", "1",
 "otau_cagk", "ms",
 "gkbar_cagk", "mho/cm2",
 "ik_cagk", "mA/cm2",
 "gkca_cagk", "mho/cm2",
 0,0
};
 static double delta_t = 0.01;
 static double o0 = 0;
 static double v = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "cai_cagk", &cai_cagk,
 "d1_cagk", &d1_cagk,
 "d2_cagk", &d2_cagk,
 "k1_cagk", &k1_cagk,
 "k2_cagk", &k2_cagk,
 "abar_cagk", &abar_cagk,
 "bbar_cagk", &bbar_cagk,
 "st_cagk", &st_cagk,
 "oinf_cagk", &oinf_cagk,
 "otau_cagk", &otau_cagk,
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
static void _ode_map(int, double**, double**, double*, Datum*, double*, int);
static void _ode_spec(NrnThread*, _Memb_list*, int);
static void _ode_matsol(NrnThread*, _Memb_list*, int);
 
#define _cvode_ieq _ppvar[6]._i
 static void _ode_matsol_instance1(_threadargsproto_);
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"cagk",
 "gkbar_cagk",
 0,
 "ik_cagk",
 "gkca_cagk",
 0,
 "o_cagk",
 0,
 0};
 static Symbol* _nca_sym;
 static Symbol* _lca_sym;
 static Symbol* _tca_sym;
 static Symbol* _k_sym;
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 10, _prop);
 	/*initialize range parameters*/
 	gkbar = 0.01;
 	_prop->param = _p;
 	_prop->param_size = 10;
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 7, _prop);
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 prop_ion = need_memb(_nca_sym);
 nrn_promote(prop_ion, 1, 0);
 	_ppvar[0]._pval = &prop_ion->param[1]; /* ncai */
 prop_ion = need_memb(_lca_sym);
 nrn_promote(prop_ion, 1, 0);
 	_ppvar[1]._pval = &prop_ion->param[1]; /* lcai */
 prop_ion = need_memb(_tca_sym);
 nrn_promote(prop_ion, 1, 0);
 	_ppvar[2]._pval = &prop_ion->param[1]; /* tcai */
 prop_ion = need_memb(_k_sym);
 nrn_promote(prop_ion, 0, 1);
 	_ppvar[3]._pval = &prop_ion->param[0]; /* ek */
 	_ppvar[4]._pval = &prop_ion->param[3]; /* ik */
 	_ppvar[5]._pval = &prop_ion->param[4]; /* _ion_dikdv */
 
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
extern void _nrn_thread_table_reg(int, void(*)(double*, Datum*, Datum*, NrnThread*, int));
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 void _CaBK_reg() {
	int _vectorized = 0;
  _initlists();
 	ion_reg("nca", 2.0);
 	ion_reg("lca", 2.0);
 	ion_reg("tca", 2.0);
 	ion_reg("k", -10000.);
 	_nca_sym = hoc_lookup("nca_ion");
 	_lca_sym = hoc_lookup("lca_ion");
 	_tca_sym = hoc_lookup("tca_ion");
 	_k_sym = hoc_lookup("k_ion");
 	register_mech(_mechanism, nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init, hoc_nrnpointerindex, 0);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
     _nrn_thread_reg(_mechtype, 2, _update_ion_pointer);
 #if NMODL_TEXT
  hoc_reg_nmodl_text(_mechtype, nmodl_file_text);
  hoc_reg_nmodl_filename(_mechtype, nmodl_filename);
#endif
  hoc_register_prop_size(_mechtype, 10, 7);
  hoc_register_dparam_semantics(_mechtype, 0, "nca_ion");
  hoc_register_dparam_semantics(_mechtype, 1, "lca_ion");
  hoc_register_dparam_semantics(_mechtype, 2, "tca_ion");
  hoc_register_dparam_semantics(_mechtype, 3, "k_ion");
  hoc_register_dparam_semantics(_mechtype, 4, "k_ion");
  hoc_register_dparam_semantics(_mechtype, 5, "k_ion");
  hoc_register_dparam_semantics(_mechtype, 6, "cvodeieq");
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 cagk /Users/temma/ghq/pydentate/mechs/CaBK.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
 
#define FARADAY _nrnunit_FARADAY[_nrnunit_use_legacy_]
static double _nrnunit_FARADAY[2] = {0x1.81f0fae775425p+6, 96.4853}; /* 96.4853321233100161 */
 static double R = 8.313424;
static int _reset;
static char *modelname = "CaGk";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static int rate(double, double);
 
static int _ode_spec1(_threadargsproto_);
/*static int _ode_matsol1(_threadargsproto_);*/
 static int _slist1[1], _dlist1[1];
 static int state(_threadargsproto_);
 
/*CVODE*/
 static int _ode_spec1 () {_reset=0;
 {
   rate ( _threadargscomma_ v , cai ) ;
   Do = ( oinf - o ) / otau ;
   }
 return _reset;
}
 static int _ode_matsol1 () {
 rate ( _threadargscomma_ v , cai ) ;
 Do = Do  / (1. - dt*( ( ( ( - 1.0 ) ) ) / otau )) ;
  return 0;
}
 /*END CVODE*/
 static int state () {_reset=0;
 {
   rate ( _threadargscomma_ v , cai ) ;
    o = o + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / otau)))*(- ( ( ( oinf ) ) / otau ) / ( ( ( ( - 1.0 ) ) ) / otau ) - o) ;
   }
  return 0;
}
 
double alp (  double _lv , double _lc ) {
   double _lalp;
 _lalp = _lc * abar / ( _lc + exp1 ( _threadargscomma_ k1 , d1 , _lv ) ) ;
   
return _lalp;
 }
 
static void _hoc_alp(void) {
  double _r;
   _r =  alp (  *getarg(1) , *getarg(2) );
 hoc_retpushx(_r);
}
 
double bet (  double _lv , double _lc ) {
   double _lbet;
 _lbet = bbar / ( 1.0 + _lc / exp1 ( _threadargscomma_ k2 , d2 , _lv ) ) ;
   
return _lbet;
 }
 
static void _hoc_bet(void) {
  double _r;
   _r =  bet (  *getarg(1) , *getarg(2) );
 hoc_retpushx(_r);
}
 
double exp1 (  double _lk , double _ld , double _lv ) {
   double _lexp1;
 _lexp1 = _lk * exp ( - 2.0 * _ld * FARADAY * _lv / R / ( 273.15 + celsius ) ) ;
   
return _lexp1;
 }
 
static void _hoc_exp1(void) {
  double _r;
   _r =  exp1 (  *getarg(1) , *getarg(2) , *getarg(3) );
 hoc_retpushx(_r);
}
 
static int  rate (  double _lv , double _lc ) {
   double _la ;
 _la = alp ( _threadargscomma_ _lv , _lc ) ;
   otau = 1.0 / ( _la + bet ( _threadargscomma_ _lv , _lc ) ) ;
   oinf = _la * otau ;
    return 0; }
 
static void _hoc_rate(void) {
  double _r;
   _r = 1.;
 rate (  *getarg(1) , *getarg(2) );
 hoc_retpushx(_r);
}
 
static int _ode_count(int _type){ return 1;}
 
static void _ode_spec(NrnThread* _nt, _Memb_list* _ml, int _type) {
   Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
  ncai = _ion_ncai;
  lcai = _ion_lcai;
  tcai = _ion_tcai;
  ek = _ion_ek;
     _ode_spec1 ();
  }}
 
static void _ode_map(int _ieq, double** _pv, double** _pvdot, double* _pp, Datum* _ppd, double* _atol, int _type) { 
 	int _i; _p = _pp; _ppvar = _ppd;
	_cvode_ieq = _ieq;
	for (_i=0; _i < 1; ++_i) {
		_pv[_i] = _pp + _slist1[_i];  _pvdot[_i] = _pp + _dlist1[_i];
		_cvode_abstol(_atollist, _atol, _i);
	}
 }
 
static void _ode_matsol_instance1(_threadargsproto_) {
 _ode_matsol1 ();
 }
 
static void _ode_matsol(NrnThread* _nt, _Memb_list* _ml, int _type) {
   Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
  ncai = _ion_ncai;
  lcai = _ion_lcai;
  tcai = _ion_tcai;
  ek = _ion_ek;
 _ode_matsol_instance1(_threadargs_);
 }}
 extern void nrn_update_ion_pointer(Symbol*, Datum*, int, int);
 static void _update_ion_pointer(Datum* _ppvar) {
   nrn_update_ion_pointer(_nca_sym, _ppvar, 0, 1);
   nrn_update_ion_pointer(_lca_sym, _ppvar, 1, 1);
   nrn_update_ion_pointer(_tca_sym, _ppvar, 2, 1);
   nrn_update_ion_pointer(_k_sym, _ppvar, 3, 0);
   nrn_update_ion_pointer(_k_sym, _ppvar, 4, 3);
   nrn_update_ion_pointer(_k_sym, _ppvar, 5, 4);
 }

static void initmodel() {
  int _i; double _save;_ninits++;
 _save = t;
 t = 0.0;
{
  o = o0;
 {
   cai = ncai + lcai + tcai ;
   rate ( _threadargscomma_ v , cai ) ;
   o = oinf ;
   }
  _sav_indep = t; t = _save;

}
}

static void nrn_init(NrnThread* _nt, _Memb_list* _ml, int _type){
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
  ncai = _ion_ncai;
  lcai = _ion_lcai;
  tcai = _ion_tcai;
  ek = _ion_ek;
 initmodel();
 }}

static double _nrn_current(double _v){double _current=0.;v=_v;{ {
   gkca = gkbar * pow( o , st ) ;
   ik = gkca * ( v - ek ) ;
   }
 _current += ik;

} return _current;
}

static void nrn_cur(NrnThread* _nt, _Memb_list* _ml, int _type){
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
  ncai = _ion_ncai;
  lcai = _ion_lcai;
  tcai = _ion_tcai;
  ek = _ion_ek;
 _g = _nrn_current(_v + .001);
 	{ double _dik;
  _dik = ik;
 _rhs = _nrn_current(_v);
  _ion_dikdv += (_dik - ik)/.001 ;
 	}
 _g = (_g - _rhs)/.001;
  _ion_ik += ik ;
#if CACHEVEC
  if (use_cachevec) {
	VEC_RHS(_ni[_iml]) -= _rhs;
  }else
#endif
  {
	NODERHS(_nd) -= _rhs;
  }
 
}}

static void nrn_jacob(NrnThread* _nt, _Memb_list* _ml, int _type){
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

static void nrn_state(NrnThread* _nt, _Memb_list* _ml, int _type){
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
  ncai = _ion_ncai;
  lcai = _ion_lcai;
  tcai = _ion_tcai;
  ek = _ion_ek;
 { error =  state();
 if(error){fprintf(stderr,"at line 66 in file CaBK.mod:\n	SOLVE state METHOD cnexp\n"); nrn_complain(_p); abort_run(error);}
 } }}

}

static void terminal(){}

static void _initlists() {
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = o_columnindex;  _dlist1[0] = Do_columnindex;
_first = 0;
}

#if NMODL_TEXT
static const char* nmodl_filename = "/Users/temma/ghq/pydentate/mechs/CaBK.mod";
static const char* nmodl_file_text = 
  "TITLE CaGk\n"
  ": Calcium activated K channel.\n"
  ": Modified from Moczydlowski and Latorre (1983) J. Gen. Physiol. 82\n"
  "\n"
  "UNITS {\n"
  "	(molar) = (1/liter)\n"
  "}\n"
  "\n"
  "UNITS {\n"
  "	(mV) =	(millivolt)\n"
  "	(mA) =	(milliamp)\n"
  "	(mM) =	(millimolar)\n"
  "}\n"
  "\n"
  "\n"
  "NEURON {\n"
  "	SUFFIX cagk\n"
  "	USEION nca READ ncai VALENCE 2\n"
  "	USEION lca READ lcai VALENCE 2\n"
  "	USEION tca READ tcai VALENCE 2\n"
  "	USEION k READ ek WRITE ik\n"
  "	RANGE gkbar,gkca, ik\n"
  "	GLOBAL oinf, otau\n"
  "}\n"
  "\n"
  "UNITS {\n"
  "	FARADAY = (faraday)  (kilocoulombs)\n"
  "	R = 8.313424 (joule/degC)\n"
  "}\n"
  "\n"
  "PARAMETER {\n"
  "	celsius		(degC)\n"
  "	v		(mV)\n"
  "	gkbar=.01	(mho/cm2)	: Maximum Permeability\n"
  "	cai = 5.e-5	(mM)\n"
  "	ek		(mV)\n"
  "\n"
  "	d1 = .84\n"
  "	d2 = 1.\n"
  "	k1 = .48e-3	(mM)\n"
  "	k2 = .13e-6	(mM)\n"
  "	abar = .28	(/ms)\n"
  "	bbar = .48	(/ms)\n"
  "        st=1            (1)\n"
  "	lcai		(mV)\n"
  "	ncai		(mV)\n"
  "	tcai		(mV)\n"
  "}\n"
  "\n"
  "ASSIGNED {\n"
  "	ik		(mA/cm2)\n"
  "	oinf\n"
  "	otau		(ms)\n"
  "        gkca          (mho/cm2)\n"
  "}\n"
  "\n"
  "INITIAL {\n"
  "	cai= ncai + lcai + tcai\n"
  "        rate(v,cai)\n"
  "        o=oinf\n"
  "}\n"
  "\n"
  "STATE {	o }		: fraction of open channels\n"
  "\n"
  "BREAKPOINT {\n"
  "	SOLVE state METHOD cnexp\n"
  "	gkca = gkbar*o^st\n"
  "	ik = gkca*(v - ek)\n"
  "}\n"
  "\n"
  "DERIVATIVE state {	: exact when v held constant; integrates over dt step\n"
  "	rate(v, cai)\n"
  "	o' = (oinf - o)/otau\n"
  "}\n"
  "\n"
  "FUNCTION alp(v (mV), c (mM)) (1/ms) { :callable from hoc\n"
  "	alp = c*abar/(c + exp1(k1,d1,v))\n"
  "}\n"
  "\n"
  "FUNCTION bet(v (mV), c (mM)) (1/ms) { :callable from hoc\n"
  "	bet = bbar/(1 + c/exp1(k2,d2,v))\n"
  "}\n"
  "\n"
  "FUNCTION exp1(k (mM), d, v (mV)) (mM) { :callable from hoc\n"
  "	exp1 = k*exp(-2*d*FARADAY*v/R/(273.15 + celsius))\n"
  "}\n"
  "\n"
  "PROCEDURE rate(v (mV), c (mM)) { :callable from hoc\n"
  "	LOCAL a\n"
  "	a = alp(v,c)\n"
  "	otau = 1/(a + bet(v, c))\n"
  "	oinf = a*otau\n"
  "}\n"
  "\n"
  ;
#endif
