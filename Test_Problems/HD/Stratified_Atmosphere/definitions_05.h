#define  PHYSICS                 HD
#define  DIMENSIONS              3
#define  COMPONENTS              3
#define  GEOMETRY                CARTESIAN
#define  BODY_FORCE              POTENTIAL
#define  COOLING                 NO
#define  INTERPOLATION           PARABOLIC
#define  TIME_STEPPING           CHARACTERISTIC_TRACING
#define  DIMENSIONAL_SPLITTING   NO
#define  NTRACER                 0
#define  USER_DEF_PARAMETERS     1
#define  USER_DEF_CONSTANTS      0

/* -- physics dependent declarations -- */

#define  EOS                     IDEAL
#define  ENTROPY_SWITCH          NO
#define  THERMAL_CONDUCTION      NO
#define  VISCOSITY               NO
#define  ROTATING_FRAME          NO

/* -- user-defined parameters (labels) -- */

#define  ALPHA                   0

/* -- user-defined symbolic constants -- */


/* -- supplementary constants (user editable) -- */ 

#define  INITIAL_SMOOTHING      NO
#define  WARNING_MESSAGES       YES
#define  PRINT_TO_FILE          YES
#define  INTERNAL_BOUNDARY      NO
#define  SHOCK_FLATTENING       NO
#define  ARTIFICIAL_VISCOSITY   NO
#define  CHAR_LIMITING          NO
#define  LIMITER                VANLEER_LIM
