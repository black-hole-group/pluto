#define  PHYSICS                 RHD
#define  DIMENSIONS              2
#define  COMPONENTS              2
#define  GEOMETRY                CARTESIAN
#define  BODY_FORCE              NO
#define  COOLING                 NO
#define  INTERPOLATION           LINEAR
#define  TIME_STEPPING           RK2
#define  DIMENSIONAL_SPLITTING   NO
#define  NTRACER                 0
#define  USER_DEF_PARAMETERS     1
#define  USER_DEF_CONSTANTS      0

/* -- physics dependent declarations -- */

#define  EOS                     IDEAL
#define  ENTROPY_SWITCH          NO
#define  USE_FOUR_VELOCITY       NO

/* -- user-defined parameters (labels) -- */

#define  SCRH                    0

/* -- user-defined symbolic constants -- */


/* -- supplementary constants (user editable) -- */ 

#define  INITIAL_SMOOTHING      NO
#define  WARNING_MESSAGES       YES
#define  PRINT_TO_FILE          NO
#define  INTERNAL_BOUNDARY      NO
#define  SHOCK_FLATTENING       NO
#define  CHOMBO_EN_SWITCH       NO
#define  CHOMBO_REF_VAR         ENG
#define  CHOMBO_LOGR            NO
#define  ARTIFICIAL_VISCOSITY   NO
#define  CHAR_LIMITING          NO
#define  LIMITER                VANALBADA_LIM
