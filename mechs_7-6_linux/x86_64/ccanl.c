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
 
#define nrn_init _nrn_init__ccanl
#define _nrn_initial _nrn_initial__ccanl
#define nrn_cur _nrn_cur__ccanl
#define _nrn_current _nrn_current__ccanl
#define nrn_jacob _nrn_jacob__ccanl
#define nrn_state _nrn_state__ccanl
#define _net_receive _net_receive__ccanl 
#define integrate integrate__ccanl 
 
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
#define catau _p[0]
#define caiinf _p[1]
#define cai _p[2]
#define eca _p[3]
#define inca _p[4]
#define ilca _p[5]
#define itca _p[6]
#define enca _p[7]
#define elca _p[8]
#define etca _p[9]
#define ncai _p[10]
#define Dncai _p[11]
#define lcai _p[12]
#define Dlcai _p[13]
#define tcai _p[14]
#define Dtcai _p[15]
#define _g _p[16]
#define _ion_ncai	*_ppvar[0]._pval
#define _ion_inca	*_ppvar[1]._pval
#define _ion_enca	*_ppvar[2]._pval
#define _style_nca	*((int*)_ppvar[3]._pvoid)
#define _ion_lcai	*_ppvar[4]._pval
#define _ion_ilca	*_ppvar[5]._pval
#define _ion_elca	*_ppvar[6]._pval
#define _style_lca	*((int*)_ppvar[7]._pvoid)
#define _ion_tcai	*_ppvar[8]._pval
#define _ion_itca	*_ppvar[9]._pval
#define _ion_etca	*_ppvar[10]._pval
#define _style_tca	*((int*)_ppvar[11]._pvoid)
 
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
 static void _hoc_ktf(void);
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
 "setdata_ccanl", _hoc_setdata,
 "ktf_ccanl", _hoc_ktf,
 0, 0
};
#define ktf ktf_ccanl
 extern double ktf( );
 /* declare global and static user variables */
#define cao cao_ccanl
 double cao = 2;
#define depth depth_ccanl
 double depth = 200;
#define ica ica_ccanl
 double ica = 0;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "depth_ccanl", "nm",
 "cao_ccanl", "mM",
 "ica_ccanl", "mA/cm2",
 "catau_ccanl", "ms",
 "caiinf_ccanl", "mM",
 "cai_ccanl", "mM",
 "eca_ccanl", "mV",
 0,0
};
 static double delta_t = 1;
 static double lcai0 = 0;
 static double ncai0 = 0;
 static double tcai0 = 0;
 static double v = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "depth_ccanl", &depth_ccanl,
 "cao_ccanl", &cao_ccanl,
 "ica_ccanl", &ica_ccanl,
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
 
#define _cvode_ieq _ppvar[12]._i
 static void _ode_matsol_instance1(_threadargsproto_);
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.5.0",
"ccanl",
 "catau_ccanl",
 "caiinf_ccanl",
 "cai_ccanl",
 0,
 "eca_ccanl",
 0,
 0,
 0};
 static Symbol* _nca_sym;
 static Symbol* _lca_sym;
 static Symbol* _tca_sym;
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 17, _prop);
 	/*initialize range parameters*/
 	catau = 9;
 	caiinf = 5e-05;
 	cai = 5e-05;
 	_prop->param = _p;
 	_prop->param_size = 17;
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 13, _prop);
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 prop_ion = need_memb(_nca_sym);
 nrn_check_conc_write(_prop, prop_ion, 1);
 nrn_promote(prop_ion, 3, 3);
 	_ppvar[0]._pval = &prop_ion->param[1]; /* ncai */
 	_ppvar[1]._pval = &prop_ion->param[3]; /* inca */
 	_ppvar[2]._pval = &prop_ion->param[0]; /* enca */
 	_ppvar[3]._pvoid = (void*)(&(prop_ion->dparam[0]._i)); /* iontype for nca */
 prop_ion = need_memb(_lca_sym);
 nrn_check_conc_write(_prop, prop_ion, 1);
 nrn_promote(prop_ion, 3, 3);
 	_ppvar[4]._pval = &prop_ion->param[1]; /* lcai */
 	_ppvar[5]._pval = &prop_ion->param[3]; /* ilca */
 	_ppvar[6]._pval = &prop_ion->param[0]; /* elca */
 	_ppvar[7]._pvoid = (void*)(&(prop_ion->dparam[0]._i)); /* iontype for lca */
 prop_ion = need_memb(_tca_sym);
 nrn_check_conc_write(_prop, prop_ion, 1);
 nrn_promote(prop_ion, 3, 3);
 	_ppvar[8]._pval = &prop_ion->param[1]; /* tcai */
 	_ppvar[9]._pval = &prop_ion->param[3]; /* itca */
 	_ppvar[10]._pval = &prop_ion->param[0]; /* etca */
 	_ppvar[11]._pvoid = (void*)(&(prop_ion->dparam[0]._i)); /* iontype for tca */
 
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

 void _ccanl_reg() {
	int _vectorized = 0;
  _initlists();
 	ion_reg("nca", 2.0);
 	ion_reg("lca", 2.0);
 	ion_reg("tca", 2.0);
 	_nca_sym = hoc_lookup("nca_ion");
 	_lca_sym = hoc_lookup("lca_ion");
 	_tca_sym = hoc_lookup("tca_ion");
 	register_mech(_mechanism, nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init, hoc_nrnpointerindex, 0);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
     _nrn_thread_reg(_mechtype, 2, _update_ion_pointer);
  hoc_register_prop_size(_mechtype, 17, 13);
  hoc_register_dparam_semantics(_mechtype, 0, "nca_ion");
  hoc_register_dparam_semantics(_mechtype, 1, "nca_ion");
  hoc_register_dparam_semantics(_mechtype, 2, "nca_ion");
  hoc_register_dparam_semantics(_mechtype, 3, "#nca_ion");
  hoc_register_dparam_semantics(_mechtype, 4, "lca_ion");
  hoc_register_dparam_semantics(_mechtype, 5, "lca_ion");
  hoc_register_dparam_semantics(_mechtype, 6, "lca_ion");
  hoc_register_dparam_semantics(_mechtype, 7, "#lca_ion");
  hoc_register_dparam_semantics(_mechtype, 8, "tca_ion");
  hoc_register_dparam_semantics(_mechtype, 9, "tca_ion");
  hoc_register_dparam_semantics(_mechtype, 10, "tca_ion");
  hoc_register_dparam_semantics(_mechtype, 11, "#tca_ion");
  hoc_register_dparam_semantics(_mechtype, 12, "cvodeieq");
 	nrn_writes_conc(_mechtype, 0);
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 ccanl /home/daniel/repos/pyDentate/mechs_7-6/x86_64/ccanl.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
 static double FARADAY = 96520.0;
 static double R = 8.3134;
static int _reset;
static char *modelname = "";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
 static int _deriv1_advance = 0;
 
static int _ode_spec1(_threadargsproto_);
/*static int _ode_matsol1(_threadargsproto_);*/
 static int _slist2[3]; static double _dlist2[3];
 static double _savstate1[3], *_temp1 = _savstate1;
 static int _slist1[3], _dlist1[3];
 static int integrate(_threadargsproto_);
 
/*CVODE*/
 static int _ode_spec1 () {_reset=0;
 {
   Dncai = - ( inca ) / depth / FARADAY * ( 1e7 ) + ( caiinf / 3.0 - ncai ) / catau ;
   Dlcai = - ( ilca ) / depth / FARADAY * ( 1e7 ) + ( caiinf / 3.0 - lcai ) / catau ;
   Dtcai = - ( itca ) / depth / FARADAY * ( 1e7 ) + ( caiinf / 3.0 - tcai ) / catau ;
   }
 return _reset;
}
 static int _ode_matsol1 () {
 Dncai = Dncai  / (1. - dt*( ( ( ( - 1.0 ) ) ) / catau )) ;
 Dlcai = Dlcai  / (1. - dt*( ( ( ( - 1.0 ) ) ) / catau )) ;
 Dtcai = Dtcai  / (1. - dt*( ( ( ( - 1.0 ) ) ) / catau )) ;
  return 0;
}
 /*END CVODE*/
 
static int integrate () {_reset=0;
 { static int _recurse = 0;
 int _counte = -1;
 if (!_recurse) {
 _recurse = 1;
 {int _id; for(_id=0; _id < 3; _id++) { _savstate1[_id] = _p[_slist1[_id]];}}
 error = newton(3,_slist2, _p, integrate, _dlist2);
 _recurse = 0; if(error) {abort_run(error);}}
 {
   Dncai = - ( inca ) / depth / FARADAY * ( 1e7 ) + ( caiinf / 3.0 - ncai ) / catau ;
   Dlcai = - ( ilca ) / depth / FARADAY * ( 1e7 ) + ( caiinf / 3.0 - lcai ) / catau ;
   Dtcai = - ( itca ) / depth / FARADAY * ( 1e7 ) + ( caiinf / 3.0 - tcai ) / catau ;
   {int _id; for(_id=0; _id < 3; _id++) {
if (_deriv1_advance) {
 _dlist2[++_counte] = _p[_dlist1[_id]] - (_p[_slist1[_id]] - _savstate1[_id])/dt;
 }else{
_dlist2[++_counte] = _p[_slist1[_id]] - _savstate1[_id];}}}
 } }
 return _reset;}
 
double ktf (  ) {
   double _lktf;
 _lktf = ( 1000.0 ) * R * ( celsius + 273.15 ) / ( 2.0 * FARADAY ) ;
   
return _lktf;
 }
 
static void _hoc_ktf(void) {
  double _r;
   _r =  ktf (  );
 hoc_retpushx(_r);
}
 
static int _ode_count(int _type){ return 3;}
 
static void _ode_spec(_NrnThread* _nt, _Memb_list* _ml, int _type) {
   Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
  ncai = _ion_ncai;
  inca = _ion_inca;
  enca = _ion_enca;
  ncai = _ion_ncai;
  lcai = _ion_lcai;
  ilca = _ion_ilca;
  elca = _ion_elca;
  lcai = _ion_lcai;
  tcai = _ion_tcai;
  itca = _ion_itca;
  etca = _ion_etca;
  tcai = _ion_tcai;
     _ode_spec1 ();
  _ion_enca = enca;
  _ion_ncai = ncai;
  _ion_elca = elca;
  _ion_lcai = lcai;
  _ion_etca = etca;
  _ion_tcai = tcai;
 }}
 
static void _ode_map(int _ieq, double** _pv, double** _pvdot, double* _pp, Datum* _ppd, double* _atol, int _type) { 
 	int _i; _p = _pp; _ppvar = _ppd;
	_cvode_ieq = _ieq;
	for (_i=0; _i < 3; ++_i) {
		_pv[_i] = _pp + _slist1[_i];  _pvdot[_i] = _pp + _dlist1[_i];
		_cvode_abstol(_atollist, _atol, _i);
	}
 	_pv[0] = &(_ion_ncai);
 	_pv[1] = &(_ion_lcai);
 	_pv[2] = &(_ion_tcai);
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
  ncai = _ion_ncai;
  inca = _ion_inca;
  enca = _ion_enca;
  ncai = _ion_ncai;
  lcai = _ion_lcai;
  ilca = _ion_ilca;
  elca = _ion_elca;
  lcai = _ion_lcai;
  tcai = _ion_tcai;
  itca = _ion_itca;
  etca = _ion_etca;
  tcai = _ion_tcai;
 _ode_matsol_instance1(_threadargs_);
 }}
 extern void nrn_update_ion_pointer(Symbol*, Datum*, int, int);
 static void _update_ion_pointer(Datum* _ppvar) {
   nrn_update_ion_pointer(_nca_sym, _ppvar, 0, 1);
   nrn_update_ion_pointer(_nca_sym, _ppvar, 1, 3);
   nrn_update_ion_pointer(_nca_sym, _ppvar, 2, 0);
   nrn_update_ion_pointer(_lca_sym, _ppvar, 4, 1);
   nrn_update_ion_pointer(_lca_sym, _ppvar, 5, 3);
   nrn_update_ion_pointer(_lca_sym, _ppvar, 6, 0);
   nrn_update_ion_pointer(_tca_sym, _ppvar, 8, 1);
   nrn_update_ion_pointer(_tca_sym, _ppvar, 9, 3);
   nrn_update_ion_pointer(_tca_sym, _ppvar, 10, 0);
 }

static void initmodel() {
  int _i; double _save;_ninits++;
 _save = t;
 t = 0.0;
{
 {
   
/*VERBATIM*/
	ncai = _ion_ncai;
	lcai = _ion_lcai;
	tcai = _ion_tcai; 
 ncai = caiinf / 3.0 ;
   lcai = caiinf / 3.0 ;
   tcai = caiinf / 3.0 ;
   cai = caiinf ;
   eca = 130.0 ;
   enca = eca ;
   elca = eca ;
   etca = eca ;
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
  ncai = _ion_ncai;
  inca = _ion_inca;
  enca = _ion_enca;
  ncai = _ion_ncai;
  lcai = _ion_lcai;
  ilca = _ion_ilca;
  elca = _ion_elca;
  lcai = _ion_lcai;
  tcai = _ion_tcai;
  itca = _ion_itca;
  etca = _ion_etca;
  tcai = _ion_tcai;
 initmodel();
  _ion_enca = enca;
  _ion_ncai = ncai;
  nrn_wrote_conc(_nca_sym, (&(_ion_ncai)) - 1, _style_nca);
  _ion_elca = elca;
  _ion_lcai = lcai;
  nrn_wrote_conc(_lca_sym, (&(_ion_lcai)) - 1, _style_lca);
  _ion_etca = etca;
  _ion_tcai = tcai;
  nrn_wrote_conc(_tca_sym, (&(_ion_tcai)) - 1, _style_tca);
}}

static double _nrn_current(double _v){double _current=0.;v=_v;{
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
double _dtsav = dt;
if (secondorder) { dt *= 0.5; }
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
  inca = _ion_inca;
  enca = _ion_enca;
  ncai = _ion_ncai;
  lcai = _ion_lcai;
  ilca = _ion_ilca;
  elca = _ion_elca;
  lcai = _ion_lcai;
  tcai = _ion_tcai;
  itca = _ion_itca;
  etca = _ion_etca;
  tcai = _ion_tcai;
 { error = _deriv1_advance = 1;
 derivimplicit(_ninits, 3, _slist1, _dlist1, _p, &t, dt, integrate, &_temp1);
_deriv1_advance = 0;
 if(error){fprintf(stderr,"at line 72 in file ccanl.mod:\n	SOLVE integrate METHOD derivimplicit\n"); nrn_complain(_p); abort_run(error);}
    if (secondorder) {
    int _i;
    for (_i = 0; _i < 3; ++_i) {
      _p[_slist1[_i]] += dt*_p[_dlist1[_i]];
    }}
 } {
   cai = ncai + lcai + tcai ;
   eca = ktf ( _threadargs_ ) * log ( cao / cai ) ;
   enca = eca ;
   elca = eca ;
   etca = eca ;
   }
  _ion_enca = enca;
  _ion_ncai = ncai;
  _ion_elca = elca;
  _ion_lcai = lcai;
  _ion_etca = etca;
  _ion_tcai = tcai;
}}
 dt = _dtsav;
}

static void terminal(){}

static void _initlists() {
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = &(ncai) - _p;  _dlist1[0] = &(Dncai) - _p;
 _slist1[1] = &(lcai) - _p;  _dlist1[1] = &(Dlcai) - _p;
 _slist1[2] = &(tcai) - _p;  _dlist1[2] = &(Dtcai) - _p;
 _slist2[0] = &(lcai) - _p;
 _slist2[1] = &(ncai) - _p;
 _slist2[2] = &(tcai) - _p;
_first = 0;
}
