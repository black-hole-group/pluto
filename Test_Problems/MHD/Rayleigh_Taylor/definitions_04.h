#define  PHYSICS                 HD
#define  DIMENSIONS              2
#define  COMPONENTS              2
#define  GEOMETRY                CARTESIAN
#define  BODY_FORCE              POTENTIAL
#define  COOLING                 NO
#define  INTERPOLATION           LINEAR
#define  TIME_STEPPING           RK2
#define  DIMENSIONAL_SPLITTING   NO
#define  NTRACER                 0
#define  USER_DEF_PARAMETERS     3
#define  USER_DEF_CONSTANTS      1

/* -- physics dependent declarations -- */

#define  EOS                     IDEAL
#define  ENTROPY_SWITCH          NO
#define  THERMAL_CONDUCTION      NO
#define  VISCOSITY               NO
#define  ROTATING_FRAME          NO

/* -- user-defined parameters (labels) -- */

#define  ETA                     0
#define  GRAV                    1
#define  CHI                     2

/* -- user-defined symbolic constants -- */

#define  USE_RANDOM_PERTURBATION   NO

/* -- supplementary constants (user editable) -- */ 

#define  INITIAL_SMOOTHING      NO
#define  WARNING_MESSAGES       YES
#define  PRINT_TO_FILE          YES
#define  INTERNAL_BOUNDARY      NO
#define  SHOCK_FLATTENING       NO
#define  CHOMBO_EN_SWITCH       NO
#define  CHOMBO_REF_VAR         RHO
#define  CHOMBO_LOGR            NO
#define  ARTIFICIAL_VISCOSITY   NO
#define  CHAR_LIMITING          YES
#define  LIMITER                MC_LIM
