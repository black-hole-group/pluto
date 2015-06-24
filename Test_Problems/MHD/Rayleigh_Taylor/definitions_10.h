#define  PHYSICS                 MHD
#define  DIMENSIONS              3
#define  COMPONENTS              3
#define  GEOMETRY                CARTESIAN
#define  BODY_FORCE              VECTOR
#define  COOLING                 NO
#define  INTERPOLATION           LINEAR
#define  TIME_STEPPING           HANCOCK
#define  DIMENSIONAL_SPLITTING   NO
#define  NTRACER                 0
#define  USER_DEF_PARAMETERS     3
#define  USER_DEF_CONSTANTS      1

/* -- physics dependent declarations -- */

#define  EOS                     IDEAL
#define  ENTROPY_SWITCH          NO
#define  MHD_FORMULATION         DIV_CLEANING
#define  BACKGROUND_FIELD        NO
#define  RESISTIVE_MHD           NO
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

#define  INITIAL_SMOOTHING         NO
#define  WARNING_MESSAGES          YES
#define  PRINT_TO_FILE             YES
#define  INTERNAL_BOUNDARY         NO
#define  SHOCK_FLATTENING          NO
#define  CHOMBO_EN_SWITCH          NO
#define  CHOMBO_REF_VAR            RHO
#define  CHOMBO_LOGR               NO
#define  ARTIFICIAL_VISCOSITY      NO
#define  CHAR_LIMITING             YES
#define  LIMITER                   MC_LIM
#define  ASSIGN_VECTOR_POTENTIAL   NO
#define  PRIMITIVE_HANCOCK         YES
