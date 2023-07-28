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
 
#define nrn_init _nrn_init__tmgexp2syn
#define _nrn_initial _nrn_initial__tmgexp2syn
#define nrn_cur _nrn_cur__tmgexp2syn
#define _nrn_current _nrn_current__tmgexp2syn
#define nrn_jacob _nrn_jacob__tmgexp2syn
#define nrn_state _nrn_state__tmgexp2syn
#define _net_receive _net_receive__tmgexp2syn 
#define state state__tmgexp2syn 
 
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
#define e _p[0]
#define e_columnindex 0
#define tau_1 _p[1]
#define tau_1_columnindex 1
#define tau_2 _p[2]
#define tau_2_columnindex 2
#define tau_rec _p[3]
#define tau_rec_columnindex 3
#define tau_facil _p[4]
#define tau_facil_columnindex 4
#define U _p[5]
#define U_columnindex 5
#define u0 _p[6]
#define u0_columnindex 6
#define i _p[7]
#define i_columnindex 7
#define g _p[8]
#define g_columnindex 8
#define A _p[9]
#define A_columnindex 9
#define B _p[10]
#define B_columnindex 10
#define x _p[11]
#define x_columnindex 11
#define factor _p[12]
#define factor_columnindex 12
#define DA _p[13]
#define DA_columnindex 13
#define DB _p[14]
#define DB_columnindex 14
#define v _p[15]
#define v_columnindex 15
#define _g _p[16]
#define _g_columnindex 16
#define _tsav _p[17]
#define _tsav_columnindex 17
#define _nd_area  *_ppvar[0]._pval
 
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
 /* declaration of user functions */
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

 extern Prop* nrn_point_prop_;
 static int _pointtype;
 static void* _hoc_create_pnt(Object* _ho) { void* create_point_process(int, Object*);
 return create_point_process(_pointtype, _ho);
}
 static void _hoc_destroy_pnt(void*);
 static double _hoc_loc_pnt(void* _vptr) {double loc_point_process(int, void*);
 return loc_point_process(_pointtype, _vptr);
}
 static double _hoc_has_loc(void* _vptr) {double has_loc_point(void*);
 return has_loc_point(_vptr);
}
 static double _hoc_get_loc_pnt(void* _vptr) {
 double get_loc_point_process(void*); return (get_loc_point_process(_vptr));
}
 extern void _nrn_setdata_reg(int, void(*)(Prop*));
 static void _setdata(Prop* _prop) {
 _extcall_prop = _prop;
 }
 static void _hoc_setdata(void* _vptr) { Prop* _prop;
 _prop = ((Point_process*)_vptr)->_prop;
   _setdata(_prop);
 }
 /* connect user functions to hoc names */
 static VoidFunc hoc_intfunc[] = {
 0,0
};
 static Member_func _member_func[] = {
 "loc", _hoc_loc_pnt,
 "has_loc", _hoc_has_loc,
 "get_loc", _hoc_get_loc_pnt,
 0, 0
};
 /* declare global and static user variables */
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 "U", 0, 1,
 "tau_facil", 0, 1e+09,
 "tau_rec", 1e-09, 1e+09,
 "tau_1", 1e-09, 1e+09,
 "u0", 0, 1,
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "e", "mV",
 "tau_1", "ms",
 "tau_rec", "ms",
 "tau_facil", "ms",
 "U", "1",
 "u0", "1",
 "A", "uS",
 "B", "uS",
 "i", "nA",
 "g", "umho",
 0,0
};
 static double A0 = 0;
 static double B0 = 0;
 static double delta_t = 0.01;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
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
 static void _hoc_destroy_pnt(void* _vptr) {
   destroy_point_process(_vptr);
}
 
static int _ode_count(int);
static void _ode_map(int, double**, double**, double*, Datum*, double*, int);
static void _ode_spec(NrnThread*, _Memb_list*, int);
static void _ode_matsol(NrnThread*, _Memb_list*, int);
 
#define _cvode_ieq _ppvar[2]._i
 static void _ode_matsol_instance1(_threadargsproto_);
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"tmgexp2syn",
 "e",
 "tau_1",
 "tau_2",
 "tau_rec",
 "tau_facil",
 "U",
 "u0",
 0,
 "i",
 "g",
 0,
 "A",
 "B",
 0,
 0};
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
  if (nrn_point_prop_) {
	_prop->_alloc_seq = nrn_point_prop_->_alloc_seq;
	_p = nrn_point_prop_->param;
	_ppvar = nrn_point_prop_->dparam;
 }else{
 	_p = nrn_prop_data_alloc(_mechtype, 18, _prop);
 	/*initialize range parameters*/
 	e = -90;
 	tau_1 = 2;
 	tau_2 = 20;
 	tau_rec = 100;
 	tau_facil = 1000;
 	U = 0.04;
 	u0 = 0;
  }
 	_prop->param = _p;
 	_prop->param_size = 18;
  if (!nrn_point_prop_) {
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 3, _prop);
  }
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 
}
 static void _initlists();
  /* some states have an absolute tolerance */
 static Symbol** _atollist;
 static HocStateTolerance _hoc_state_tol[] = {
 0,0
};
 static void _net_receive(Point_process*, double*, double);
 static void _net_init(Point_process*, double*, double);
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
extern void _nrn_thread_table_reg(int, void(*)(double*, Datum*, Datum*, NrnThread*, int));
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 void _tmgexp2syn_reg() {
	int _vectorized = 1;
  _initlists();
 	_pointtype = point_register_mech(_mechanism,
	 nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init,
	 hoc_nrnpointerindex, 1,
	 _hoc_create_pnt, _hoc_destroy_pnt, _member_func);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
 #if NMODL_TEXT
  hoc_reg_nmodl_text(_mechtype, nmodl_file_text);
  hoc_reg_nmodl_filename(_mechtype, nmodl_filename);
#endif
  hoc_register_prop_size(_mechtype, 18, 3);
  hoc_register_dparam_semantics(_mechtype, 0, "area");
  hoc_register_dparam_semantics(_mechtype, 1, "pntproc");
  hoc_register_dparam_semantics(_mechtype, 2, "cvodeieq");
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 pnt_receive[_mechtype] = _net_receive;
 pnt_receive_init[_mechtype] = _net_init;
 pnt_receive_size[_mechtype] = 5;
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 tmgexp2syn /Users/temma/ghq/pydentate/mechs/tmgexp2syn.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
static int _reset;
static char *modelname = "";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
 
static int _ode_spec1(_threadargsproto_);
/*static int _ode_matsol1(_threadargsproto_);*/
 static int _slist1[2], _dlist1[2];
 static int state(_threadargsproto_);
 
/*CVODE*/
 static int _ode_spec1 (double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt) {int _reset = 0; {
   DA = - A / tau_1 ;
   DB = - B / tau_2 ;
   }
 return _reset;
}
 static int _ode_matsol1 (double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt) {
 DA = DA  / (1. - dt*( ( - 1.0 ) / tau_1 )) ;
 DB = DB  / (1. - dt*( ( - 1.0 ) / tau_2 )) ;
  return 0;
}
 /*END CVODE*/
 static int state (double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt) { {
    A = A + (1. - exp(dt*(( - 1.0 ) / tau_1)))*(- ( 0.0 ) / ( ( - 1.0 ) / tau_1 ) - A) ;
    B = B + (1. - exp(dt*(( - 1.0 ) / tau_2)))*(- ( 0.0 ) / ( ( - 1.0 ) / tau_2 ) - B) ;
   }
  return 0;
}
 
static void _net_receive (Point_process* _pnt, double* _args, double _lflag) 
{  double* _p; Datum* _ppvar; Datum* _thread; NrnThread* _nt;
   _thread = (Datum*)0; _nt = (NrnThread*)_pnt->_vnt;   _p = _pnt->_prop->param; _ppvar = _pnt->_prop->dparam;
  if (_tsav > t){ extern char* hoc_object_name(); hoc_execerror(hoc_object_name(_pnt->ob), ":Event arrived out of order. Must call ParallelContext.set_maxstep AFTER assigning minimum NetCon.delay");}
 _tsav = t; {
   _args[2] = _args[2] * exp ( - ( t - _args[4] ) / tau_rec ) ;
   _args[2] = _args[2] + ( _args[1] * ( exp ( - ( t - _args[4] ) / tau_1 ) - exp ( - ( t - _args[4] ) / tau_rec ) ) / ( ( tau_1 / tau_rec ) - 1.0 ) ) ;
   _args[1] = _args[1] * exp ( - ( t - _args[4] ) / tau_1 ) ;
   x = 1.0 - _args[1] - _args[2] ;
   if ( tau_facil > 0.0 ) {
     _args[3] = _args[3] * exp ( - ( t - _args[4] ) / tau_facil ) ;
     }
   else {
     _args[3] = U ;
     }
   if ( tau_facil > 0.0 ) {
       _args[3] = _args[3] + U * ( 1.0 - _args[3] )  ;
     }
     g = g + _args[0] * factor * x * _args[3]  ;
     _args[1] = _args[1] + x * _args[3]  ;
   _args[4] = t ;
     if (nrn_netrec_state_adjust && !cvode_active_){
    /* discon state adjustment for cnexp case (rate uses no local variable) */
    double __state = A;
    double __primary = (A + g) - __state;
     __primary += ( 1. - exp( 0.5*dt*( ( - 1.0 ) / tau_1 ) ) )*( - ( 0.0 ) / ( ( - 1.0 ) / tau_1 ) - __primary );
    A += __primary;
  } else {
 A = A + g ;
     }
   if (nrn_netrec_state_adjust && !cvode_active_){
    /* discon state adjustment for cnexp case (rate uses no local variable) */
    double __state = B;
    double __primary = (B + g) - __state;
     __primary += ( 1. - exp( 0.5*dt*( ( - 1.0 ) / tau_2 ) ) )*( - ( 0.0 ) / ( ( - 1.0 ) / tau_2 ) - __primary );
    B += __primary;
  } else {
 B = B + g ;
     }
 } }
 
static void _net_init(Point_process* _pnt, double* _args, double _lflag) {
       double* _p = _pnt->_prop->param;
    Datum* _ppvar = _pnt->_prop->dparam;
    Datum* _thread = (Datum*)0;
    NrnThread* _nt = (NrnThread*)_pnt->_vnt;
 _args[1] = 0.0 ;
   _args[2] = 0.0 ;
   _args[3] = u0 ;
   _args[4] = t ;
   }
 
static int _ode_count(int _type){ return 2;}
 
static void _ode_spec(NrnThread* _nt, _Memb_list* _ml, int _type) {
   double* _p; Datum* _ppvar; Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
     _ode_spec1 (_p, _ppvar, _thread, _nt);
 }}
 
static void _ode_map(int _ieq, double** _pv, double** _pvdot, double* _pp, Datum* _ppd, double* _atol, int _type) { 
	double* _p; Datum* _ppvar;
 	int _i; _p = _pp; _ppvar = _ppd;
	_cvode_ieq = _ieq;
	for (_i=0; _i < 2; ++_i) {
		_pv[_i] = _pp + _slist1[_i];  _pvdot[_i] = _pp + _dlist1[_i];
		_cvode_abstol(_atollist, _atol, _i);
	}
 }
 
static void _ode_matsol_instance1(_threadargsproto_) {
 _ode_matsol1 (_p, _ppvar, _thread, _nt);
 }
 
static void _ode_matsol(NrnThread* _nt, _Memb_list* _ml, int _type) {
   double* _p; Datum* _ppvar; Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
 _ode_matsol_instance1(_threadargs_);
 }}

static void initmodel(double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt) {
  int _i; double _save;{
  A = A0;
  B = B0;
 {
   double _ltp ;
 g = 0.0 ;
   if ( tau_1 / tau_2 > 0.9999 ) {
     tau_1 = 0.9999 * tau_2 ;
     }
   if ( tau_1 / tau_2 < 1e-9 ) {
     tau_1 = tau_2 * 1e-9 ;
     }
   A = 0.0 ;
   B = 0.0 ;
   _ltp = ( tau_1 * tau_2 ) / ( tau_2 - tau_1 ) * log ( tau_2 / tau_1 ) ;
   factor = - exp ( - _ltp / tau_1 ) + exp ( - _ltp / tau_2 ) ;
   factor = 1.0 / factor ;
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
 _tsav = -1e20;
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
 initmodel(_p, _ppvar, _thread, _nt);
}
}

static double _nrn_current(double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt, double _v){double _current=0.;v=_v;{ {
   g = B - A ;
   i = g * ( v - e ) ;
   }
 _current += i;

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
 _g = _nrn_current(_p, _ppvar, _thread, _nt, _v + .001);
 	{ _rhs = _nrn_current(_p, _ppvar, _thread, _nt, _v);
 	}
 _g = (_g - _rhs)/.001;
 _g *=  1.e2/(_nd_area);
 _rhs *= 1.e2/(_nd_area);
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
 {   state(_p, _ppvar, _thread, _nt);
  }}}

}

static void terminal(){}

static void _initlists(){
 double _x; double* _p = &_x;
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = A_columnindex;  _dlist1[0] = DA_columnindex;
 _slist1[1] = B_columnindex;  _dlist1[1] = DB_columnindex;
_first = 0;
}

#if defined(__cplusplus)
} /* extern "C" */
#endif

#if NMODL_TEXT
static const char* nmodl_filename = "/Users/temma/ghq/pydentate/mechs/tmgexp2syn.mod";
static const char* nmodl_file_text = 
  "COMMENT\n"
  "Revised 12/15/2000 in light of a personal communication \n"
  "from Misha Tsodyks that u is incremented _before_ x is \n"
  "converted to y--a point that was not clear in the paper.\n"
  "If u is incremented _after_ x is converted to y, then \n"
  "the first synaptic activation after a long interval of \n"
  "silence will produce smaller and smaller postsynaptic \n"
  "effect as the length of the silent interval increases, \n"
  "eventually becoming vanishingly small.\n"
  "\n"
  "Implementation of a model of short-term facilitation and depression \n"
  "based on the kinetics described in\n"
  "  Tsodyks et al.\n"
  "  Synchrony generation in recurrent networks \n"
  "  with frequency-dependent synapses\n"
  "  Journal of Neuroscience 20:RC50:1-5, 2000.\n"
  "Their mechanism represented synapses as current sources.\n"
  "The mechanism implemented here uses a conductance change instead.\n"
  "\n"
  "The basic scheme is\n"
  "\n"
  "x -------> y    Instantaneous, spike triggered.\n"
  "                Increment is u*x (see discussion of u below).\n"
  "                x == fraction of \"synaptic resources\" that have \n"
  "                     \"recovered\" (fraction of xmtr pool that is \n"
  "                     ready for release, or fraction of postsynaptic \n"
  "                     channels that are ready to be opened, or some \n"
  "                     joint function of these two factors)\n"
  "                y == fraction of \"synaptic resources\" that are in the \n"
  "                     \"active state.\"  This is proportional to the \n"
  "                     number of channels that are open, or the \n"
  "                     fraction of max synaptic current that is \n"
  "                     being delivered. \n"
  "  tau_1\n"
  "y -------> z    z == fraction of \"synaptic resources\" that are \n"
  "                     in the \"inactive state\"\n"
  "\n"
  "  tau_rec\n"
  "z -------> x\n"
  "\n"
  "where x + y + z = 1\n"
  "\n"
  "The active state y is multiplied by a synaptic weight to compute\n"
  "the actual synaptic conductance (or current, in the original form \n"
  "of the model).\n"
  "\n"
  "In addition, there is a \"facilition\" term u that \n"
  "governs the fraction of x that is converted to y \n"
  "on each synaptic activation.\n"
  "\n"
  "  -------> u    Instantaneous, spike triggered, \n"
  "                happens _BEFORE_ x is converted to y.\n"
  "                Increment is U*(1-u) where U and u both \n"
  "                lie in the range 0 - 1.\n"
  "  tau_facil\n"
  "u ------->      decay of facilitation\n"
  "\n"
  "This implementation for NEURON offers the user a parameter \n"
  "u0 that has a default value of 0 but can be used to specify \n"
  "a nonzero initial value for u.\n"
  "\n"
  "When tau_facil = 0, u is supposed to equal U.\n"
  "\n"
  "Note that the synaptic conductance in this mechanism \n"
  "has the same kinetics as y, i.e. decays with time \n"
  "constant tau_1.\n"
  "\n"
  "This mechanism can receive multiple streams of \n"
  "synaptic input via NetCon objects.  \n"
  "Each stream keeps track of its own \n"
  "weight and activation history.\n"
  "\n"
  "The printf() statements are for testing purposes only.\n"
  "ENDCOMMENT\n"
  "\n"
  ": THE NEURON BLOCK DEFINES HOW THE MECHANISM LOOKS LIKE FROM THE OUTSIDE\n"
  ": IT TELLS US WHAT PARTS OF THE MECHANISM ARE EXPOSED TO OTHER PARTS OF THE NEURON LANGUAGE\n"
  "NEURON {\n"
  "	POINT_PROCESS tmgexp2syn\n"
  "	: THE RANGE STATEMENT ASSERTS THAT THESE VARIABLES CAN BE ACCESSED IN HOC THROUGH RANGE VARIABLE SYNTAX\n"
  "	RANGE e, i, g\n"
  "	RANGE tau_1, tau_2, tau_rec, tau_facil, U, u0\n"
  "	: THE NONSPECIFIC CURRENT DEFINITION HAS TWO CONSEQUENCES\n"
  "	: FIRST, i HAS TO BE ACCOUNTED FOR IN CHARGE BALANCE EQUATIONS\n"
  "	: SECOND, i DOES NOT CONTRIBUTE DIRECTLY TO MASS BALANCS OF ANY IONS\n"
  "	NONSPECIFIC_CURRENT i\n"
  "	: RANGE g\n"
  "}\n"
  ": UNITS BLOCK JUST DEFINES NEW NAMES FOR UNITS IN TERMS OF EXISTING UNITS IN THE UNIX UNITS DATABASE\n"
  "UNITS {\n"
  "	(nA) = (nanoamp)\n"
  "	(mV) = (millivolt)\n"
  "	(umho) = (micromho)\n"
  "}\n"
  "\n"
  ": PARAMETERS ARE NORMALLY SPECIFIED BY THE USER\n"
  ": PARAMETERS GENERALLY REMAIN CONSTANT BUT THEY CAN BE CHANGED IN MID-RUN IF NECESSARY\n"
  ": THE PARAMETERS BLOCK SETS DEFAULT VALUES\n"
  ": IF A PARAMETER DOES NOT APPEAR IN A NEURON BLOCKS RANGE STATEMENT, IT WILL HAVE GLOBAL SCOPE\n"
  ": GLOBAL SCOPE MEANS CHANGING ITS VALUE WILL AFFECT EVERY INSTANCE OF THAT MECHANISM\n"
  "PARAMETER {\n"
  "	: e = -90 mV for inhibitory synapses,\n"
  "	:     0 mV for excitatory\n"
  "	e = -90	(mV)\n"
  "	: tau_1 was the same for inhibitory and excitatory synapses\n"
  "	: in the models used by T et al.\n"
  "	tau_1 = 2 (ms) < 1e-9, 1e9 >\n"
  "	tau_2 = 20\n"
  "	: tau_rec = 100 ms for inhibitory synapses,\n"
  "	:           800 ms for excitatory\n"
  "	tau_rec = 100 (ms) < 1e-9, 1e9 >\n"
  "	: tau_facil = 1000 ms for inhibitory synapses,\n"
  "	:             0 ms for excitatory\n"
  "	tau_facil = 1000 (ms) < 0, 1e9 >\n"
  "	: U = 0.04 for inhibitory synapses, \n"
  "	:     0.5 for excitatory\n"
  "	: the (1) is needed for the < 0, 1 > to be effective\n"
  "	:   in limiting the values of U and u0\n"
  "	U = 0.04 (1) < 0, 1 >\n"
  "	: initial value for the \"facilitation variable\"\n"
  "	u0 = 0 (1) < 0, 1 >\n"
  "}\n"
  "\n"
  ": THE ASSIGNED BLOCK IS USEDTO DECLARE TWO KINDS OF VARIABLES\n"
  ": ONE TYPE THAT IS GIVEN VALUES OUTSIDE THE MOD FILE\n"
  ": THE OTHER TYPE THAT APPEARS ON THE LEFT HAND SIDE OF AN ASSIGNMENT IN THE MOD FILE\n"
  ": HERE v IS ASSIGNED OUTSIDE THE MOD FILE. IT COMES FROM THE SEGMENT THE MECHANISM IS ATTACHED TO\n"
  ": i IS ASSIGNED IN THE MOD FILE. IT IS THE CURRENT THAT FLOWS AS A CONSEQUENCE OF THE CONDUCTANCE CHANGE\n"
  ": AND THE DIFFERENCE BETWEEN v AND THE REVERSAL POTENTIAL\n"
  ": x IS ALSO ASSIGNED IN THE MOD FILE\n"
  "ASSIGNED {\n"
  "	v (mV)\n"
  "	i (nA)\n"
  "	g (umho)\n"
  "	x\n"
  "	factor\n"
  "}\n"
  "\n"
  ": IF A MODEL INVOLVES DIFFERENTIAL EQUATIONS THEIR DEPENDENT VARIABLES OR UNKNOWNS ARE TO BE LISTED\n"
  ": IN THE STATE BLOCK\n"
  "STATE {\n"
  "	A (uS)\n"
  "	B (uS)\n"
  "\n"
  "}\n"
  "\n"
  ": THE INITIAL BLOCK IS CALLED BE finitialize() TO SET GOOD INITIAL CONDITIONS OF THE MECHANISM\n"
  "INITIAL {\n"
  "	LOCAL tp\n"
  "	g=0\n"
  "	if (tau_1/tau_2 > 0.9999) {\n"
  "		tau_1 = 0.9999*tau_2\n"
  "	}\n"
  "	if (tau_1/tau_2 < 1e-9) {\n"
  "		tau_1 = tau_2*1e-9\n"
  "	}\n"
  "	A = 0\n"
  "	B = 0\n"
  "	tp = (tau_1*tau_2)/(tau_2 - tau_1) * log(tau_2/tau_1)\n"
  "	factor = -exp(-tp/tau_1) + exp(-tp/tau_2)\n"
  "	factor = 1/factor\n"
  "}\n"
  "\n"
  ": BREAKPOINT IS THE MAIN COMPUTATION BLOCK OF THE MECHANISM\n"
  ": BY THE END OF THE BLOCK ALLE VARIABLES ARE CONSISTENT WITH THE NEW TIME\n"
  ": SOLVE SPECIFIES THE MAGIC BY WHICH THE STATE IS SOLVED\n"
  "BREAKPOINT {\n"
  "	SOLVE state METHOD cnexp\n"
  "	g = B - A\n"
  "	i = g*(v - e)\n"
  "}\n"
  "\n"
  ": THE DERIVATIVE STATE BLOCK ASSIGNES VALUES TO THE DERIVATIVES\n"
  "DERIVATIVE state {\n"
  "	A' = -A/tau_1\n"
  "	B' = -B/tau_2\n"
  "	: g' = -g/tau_2\n"
  "}\n"
  "\n"
  ": THE NET_RECEIVE BLOCK SPECIFIES WHAT HAPPENS TO A NET_RECIEVE EVENT FROM A NETCON OBJECT\n"
  "NET_RECEIVE(weight (umho), y, z, u, tsyn (ms)) {\n"
  "INITIAL {\n"
  ": these are in NET_RECEIVE to be per-stream\n"
  "\n"
  "	y = 0\n"
  "	z = 0\n"
  ":	u = 0\n"
  "	u = u0\n"
  "	tsyn = t\n"
  ": this header will appear once per stream\n"
  ": printf(\"t\\t t-tsyn\\t y\\t z\\t u\\t newu\\t g\\t dg\\t newg\\t newy\\n\")\n"
  "}\n"
  "\n"
  "	: first calculate z at event-\n"
  "	:   based on prior y and z\n"
  "	z = z*exp(-(t - tsyn)/tau_rec)\n"
  "	z = z + ( y*(exp(-(t - tsyn)/tau_1) - exp(-(t - tsyn)/tau_rec)) / ((tau_1/tau_rec)-1) )\n"
  "	: now calc y at event-\n"
  "	y = y*exp(-(t - tsyn)/tau_1)\n"
  "\n"
  "	x = 1-y-z\n"
  "\n"
  "	: calc u at event--\n"
  "	if (tau_facil > 0) {\n"
  "		u = u*exp(-(t - tsyn)/tau_facil)\n"
  "	} else {\n"
  "		u = U\n"
  "	}\n"
  "\n"
  ": printf(\"%g\\t%g\\t%g\\t%g\\t%g\", t, t-tsyn, y, z, u)\n"
  "\n"
  "	if (tau_facil > 0) {\n"
  "		state_discontinuity(u, u + U*(1-u))\n"
  "	}\n"
  "\n"
  ": printf(\"\\t%g\\t%g\\t%g\", u, g, weight*factor*x*u)\n"
  "\n"
  "	state_discontinuity(g, g + weight*factor*x*u)\n"
  "	state_discontinuity(y, y + x*u)\n"
  "\n"
  "	tsyn = t\n"
  "	A = A + g\n"
  "	B = B + g\n"
  "\n"
  ": printf(\"\\t%g\\t%g\\n\", g, y)\n"
  "}\n"
  ;
#endif
