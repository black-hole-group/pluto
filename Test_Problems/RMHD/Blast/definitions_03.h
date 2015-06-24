#define  PHYSICS                 RMHD
#define  DIMENSIONS              3
#define  COMPONENTS              3
#define  GEOMETRY                CARTESIAN
#define  BODY_FORCE              NO
#define  COOLING                 NO
#define  INTERPOLATION           LINEAR
#define  TIME_STEPPING           RK2
#define  DIMENSIONAL_SPLITTING   NO
#define  NTRACER                 0
#define  USER_DEF_PARAMETERS     8
#define  USER_DEF_CONSTANTS      0

/* -- physics dependent declarations -- */

#define  EOS                     IDEAL
#define  ENTROPY_SWITCH          NO
#define  MHD_FORMULATION         CONSTRAINED_TRANSPORT

/* -- user-defined parameters (labels) -- */

#define  PRS_IN                  0
#define  PRS_OUT                 1
#define  RHO_IN                  2
#define  RHO_OUT                 3
#define  BMAG                    4
#define  THETA                   5
#define  PHI                     6
#define  RADIUS                  7

/* -- user-defined symbolic constants -- */


/* -- supplementary constants (user editable) -- */ 

#define  INITIAL_SMOOTHING         NO
#define  WARNING_MESSAGES          YES
#define  PRINT_TO_FILE             YES
#define  INTERNAL_BOUNDARY         NO
#define  SHOCK_FLATTENING          NO
#define  ARTIFICIAL_VISCOSITY      NO
#define  CHAR_LIMITING             NO
#define  LIMITER                   MINMOD_LIM
#define  CT_EMF_AVERAGE            UCT_HLL
#define  CT_EN_CORRECTION          YES
#define  ASSIGN_VECTOR_POTENTIAL   YES
#define  UPDATE_VECTOR_POTENTIAL   NO
