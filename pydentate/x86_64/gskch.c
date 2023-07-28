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
 
#define nrn_init _nrn_init__gskch
#define _nrn_initial _nrn_initial__gskch
#define nrn_cur _nrn_cur__gskch
#define _nrn_current _nrn_current__gskch
#define nrn_jacob _nrn_jacob__gskch
#define nrn_state _nrn_state__gskch
#define _net_receive _net_receive__gskch 
#define rate rate__gskch 
#define state state__gskch 
 
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
#define gskbar _p[0]
#define gskbar_columnindex 0
#define isk _p[1]
#define isk_columnindex 1
#define gsk _p[2]
#define gsk_columnindex 2
#define qinf _p[3]
#define qinf_columnindex 3
#define qtau _p[4]
#define qtau_columnindex 4
#define q _p[5]
#define q_columnindex 5
#define esk _p[6]
#define esk_columnindex 6
#define ncai _p[7]
#define ncai_columnindex 7
#define lcai _p[8]
#define lcai_columnindex 8
#define tcai _p[9]
#define tcai_columnindex 9
#define Dq _p[10]
#define Dq_columnindex 10
#define qexp _p[11]
#define qexp_columnindex 11
#define _g _p[12]
#define _g_columnindex 12
#define _ion_esk	*_ppvar[0]._pval
#define _ion_isk	*_ppvar[1]._pval
#define _ion_diskdv	*_ppvar[2]._pval
#define _ion_ncai	*_ppvar[3]._pval
#define _ion_lcai	*_ppvar[4]._pval
#define _ion_tcai	*_ppvar[5]._pval
 
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
 static void _hoc_rate(void);
 static void _hoc_state(void);
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
 "setdata_gskch", _hoc_setdata,
 "rate_gskch", _hoc_rate,
 "state_gskch", _hoc_state,
 0, 0
};
 /* declare global and static user variables */
#define cai cai_gskch
 double cai = 0;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "cai_gskch", "mM",
 "gskbar_gskch", "mho/cm2",
 "isk_gskch", "mA/cm2",
 "gsk_gskch", "mho/cm2",
 "qtau_gskch", "ms",
 0,0
};
 static double delta_t = 1;
 static double q0 = 0;
 static double v = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "cai_gskch", &cai_gskch,
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
"gskch",
 "gskbar_gskch",
 0,
 "isk_gskch",
 "gsk_gskch",
 "qinf_gskch",
 "qtau_gskch",
 0,
 "q_gskch",
 0,
 0};
 static Symbol* _sk_sym;
 static Symbol* _nca_sym;
 static Symbol* _lca_sym;
 static Symbol* _tca_sym;
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 13, _prop);
 	/*initialize range parameters*/
 	gskbar = 0;
 	_prop->param = _p;
 	_prop->param_size = 13;
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 6, _prop);
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 prop_ion = need_memb(_sk_sym);
 nrn_promote(prop_ion, 0, 1);
 	_ppvar[0]._pval = &prop_ion->param[0]; /* esk */
 	_ppvar[1]._pval = &prop_ion->param[3]; /* isk */
 	_ppvar[2]._pval = &prop_ion->param[4]; /* _ion_diskdv */
 prop_ion = need_memb(_nca_sym);
 nrn_promote(prop_ion, 1, 0);
 	_ppvar[3]._pval = &prop_ion->param[1]; /* ncai */
 prop_ion = need_memb(_lca_sym);
 nrn_promote(prop_ion, 1, 0);
 	_ppvar[4]._pval = &prop_ion->param[1]; /* lcai */
 prop_ion = need_memb(_tca_sym);
 nrn_promote(prop_ion, 1, 0);
 	_ppvar[5]._pval = &prop_ion->param[1]; /* tcai */
 
}
 static void _initlists();
 static void _update_ion_pointer(Datum*);
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
extern void _nrn_thread_table_reg(int, void(*)(double*, Datum*, Datum*, NrnThread*, int));
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 void _gskch_reg() {
	int _vectorized = 0;
  _initlists();
 	ion_reg("sk", 1.0);
 	ion_reg("nca", 2.0);
 	ion_reg("lca", 2.0);
 	ion_reg("tca", 2.0);
 	_sk_sym = hoc_lookup("sk_ion");
 	_nca_sym = hoc_lookup("nca_ion");
 	_lca_sym = hoc_lookup("lca_ion");
 	_tca_sym = hoc_lookup("tca_ion");
 	register_mech(_mechanism, nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init, hoc_nrnpointerindex, 0);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
     _nrn_thread_reg(_mechtype, 2, _update_ion_pointer);
 #if NMODL_TEXT
  hoc_reg_nmodl_text(_mechtype, nmodl_file_text);
  hoc_reg_nmodl_filename(_mechtype, nmodl_filename);
#endif
  hoc_register_prop_size(_mechtype, 13, 6);
  hoc_register_dparam_semantics(_mechtype, 0, "sk_ion");
  hoc_register_dparam_semantics(_mechtype, 1, "sk_ion");
  hoc_register_dparam_semantics(_mechtype, 2, "sk_ion");
  hoc_register_dparam_semantics(_mechtype, 3, "nca_ion");
  hoc_register_dparam_semantics(_mechtype, 4, "lca_ion");
  hoc_register_dparam_semantics(_mechtype, 5, "tca_ion");
 	hoc_register_cvode(_mechtype, _ode_count, 0, 0, 0);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 gskch /Users/temma/ghq/pydentate/mechs/gskch.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
 static double _zq10 ;
static int _reset;
static char *modelname = "gskch.mod  calcium-activated potassium channel (non-voltage-dependent)";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static int rate(double);
static int state();
 
static int  state (  ) {
   cai = ncai + lcai + tcai ;
   rate ( _threadargscomma_ cai ) ;
   q = q + ( qinf - q ) * qexp ;
    return 0; }
 
static void _hoc_state(void) {
  double _r;
   _r = 1.;
 state (  );
 hoc_retpushx(_r);
}
 
static int  rate (  double _lcai ) {
   double _lalpha , _lbeta , _ltinc ;
 _zq10 = pow( 3.0 , ( ( celsius - 6.3 ) / 10.0 ) ) ;
   _lalpha = 1.25e1 * _lcai * _lcai ;
   _lbeta = 0.00025 ;
   qtau = 1.0 / ( _lalpha + _lbeta ) ;
   qinf = _lalpha * qtau ;
   _ltinc = - dt * _zq10 ;
   qexp = 1.0 - exp ( _ltinc / qtau ) * _zq10 ;
    return 0; }
 
static void _hoc_rate(void) {
  double _r;
   _r = 1.;
 rate (  *getarg(1) );
 hoc_retpushx(_r);
}
 
static int _ode_count(int _type){ hoc_execerror("gskch", "cannot be used with CVODE"); return 0;}
 extern void nrn_update_ion_pointer(Symbol*, Datum*, int, int);
 static void _update_ion_pointer(Datum* _ppvar) {
   nrn_update_ion_pointer(_sk_sym, _ppvar, 0, 0);
   nrn_update_ion_pointer(_sk_sym, _ppvar, 1, 3);
   nrn_update_ion_pointer(_sk_sym, _ppvar, 2, 4);
   nrn_update_ion_pointer(_nca_sym, _ppvar, 3, 1);
   nrn_update_ion_pointer(_lca_sym, _ppvar, 4, 1);
   nrn_update_ion_pointer(_tca_sym, _ppvar, 5, 1);
 }

static void initmodel() {
  int _i; double _save;_ninits++;
 _save = t;
 t = 0.0;
{
  q = q0;
 {
   cai = ncai + lcai + tcai ;
   rate ( _threadargscomma_ cai ) ;
   q = qinf ;
   
/*VERBATIM*/
	ncai = _ion_ncai;
	lcai = _ion_lcai;
	tcai = _ion_tcai;
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
  esk = _ion_esk;
  ncai = _ion_ncai;
  lcai = _ion_lcai;
  tcai = _ion_tcai;
 initmodel();
 }}

static double _nrn_current(double _v){double _current=0.;v=_v;{ {
   gsk = gskbar * q * q ;
   isk = gsk * ( v - esk ) ;
   }
 _current += isk;

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
  esk = _ion_esk;
  ncai = _ion_ncai;
  lcai = _ion_lcai;
  tcai = _ion_tcai;
 _g = _nrn_current(_v + .001);
 	{ double _disk;
  _disk = isk;
 _rhs = _nrn_current(_v);
  _ion_diskdv += (_disk - isk)/.001 ;
 	}
 _g = (_g - _rhs)/.001;
  _ion_isk += isk ;
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
  esk = _ion_esk;
  ncai = _ion_ncai;
  lcai = _ion_lcai;
  tcai = _ion_tcai;
 { error =  state();
 if(error){fprintf(stderr,"at line 47 in file gskch.mod:\n	SOLVE state\n"); nrn_complain(_p); abort_run(error);}
 } }}

}

static void terminal(){}

static void _initlists() {
 int _i; static int _first = 1;
  if (!_first) return;
_first = 0;
}

#if NMODL_TEXT
static const char* nmodl_filename = "/Users/temma/ghq/pydentate/mechs/gskch.mod";
static const char* nmodl_file_text = 
  "TITLE gskch.mod  calcium-activated potassium channel (non-voltage-dependent)\n"
  "\n"
  "COMMENT\n"
  "\n"
  "gsk granule\n"
  "\n"
  "ENDCOMMENT\n"
  "\n"
  "UNITS {\n"
  "        (molar) = (1/liter)\n"
  "        (mM)    = (millimolar)\n"
  "	(mA)	= (milliamp)\n"
  "	(mV)	= (millivolt)\n"
  "}\n"
  "\n"
  "NEURON {\n"
  "	SUFFIX gskch\n"
  "	USEION sk READ esk WRITE isk VALENCE 1\n"
  "	USEION nca READ ncai VALENCE 2\n"
  "	USEION lca READ lcai VALENCE 2\n"
  "	USEION tca READ tcai VALENCE 2\n"
  "	RANGE gsk, gskbar, qinf, qtau, isk\n"
  "}\n"
  "\n"
  "INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}\n"
  "\n"
  "PARAMETER {\n"
  "	celsius=6.3 (degC)\n"
  "	v		(mV)\n"
  "	dt		(ms)\n"
  "	gskbar  (mho/cm2)\n"
  "	esk	(mV)\n"
  "	cai (mM)\n"
  "	ncai (mM)\n"
  "	lcai (mM)\n"
  "	tcai (mM)\n"
  "}\n"
  "\n"
  "STATE { q }\n"
  "\n"
  "ASSIGNED {\n"
  "	isk (mA/cm2) gsk (mho/cm2) qinf qtau (ms) qexp\n"
  "}\n"
  "\n"
  "\n"
  "BREAKPOINT {          :Computes i=g*q^2*(v-esk)\n"
  "	SOLVE state\n"
  "        gsk = gskbar * q*q\n"
  "	isk = gsk * (v-esk)\n"
  "}\n"
  "\n"
  "UNITSOFF\n"
  "\n"
  "INITIAL {\n"
  "	cai = ncai + lcai + tcai	\n"
  "	rate(cai)\n"
  "	q=qinf\n"
  "	VERBATIM\n"
  "	ncai = _ion_ncai;\n"
  "	lcai = _ion_lcai;\n"
  "	tcai = _ion_tcai;\n"
  "	ENDVERBATIM\n"
  "}\n"
  "\n"
  "\n"
  "PROCEDURE state() {  :Computes state variable q at current v and dt.\n"
  "	cai = ncai + lcai + tcai\n"
  "	rate(cai)\n"
  "	q = q + (qinf-q) * qexp\n"
  "\n"
  "}\n"
  "\n"
  "LOCAL q10\n"
  "PROCEDURE rate(cai) {  :Computes rate and other constants at current v.\n"
  "	LOCAL alpha, beta, tinc\n"
  "	q10 = 3^((celsius - 6.3)/10)\n"
  "		:\"q\" activation system\n"
  "alpha = 1.25e1 * cai * cai\n"
  "beta = 0.00025 \n"
  "\n"
  ":	alpha = 0.00246/exp((12*log10(cai)+28.48)/-4.5)\n"
  ":	beta = 0.006/exp((12*log10(cai)+60.4)/35)\n"
  ": alpha = 0.00246/fctrap(cai)\n"
  ": beta = 0.006/fctrap(cai)\n"
  "	qtau = 1 / (alpha + beta)\n"
  "	qinf = alpha * qtau\n"
  "	tinc = -dt*q10\n"
  "	qexp = 1 - exp(tinc/qtau)*q10\n"
  "}\n"
  "\n"
  "UNITSON\n"
  ;
#endif
