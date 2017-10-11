#define  PHYSICS                 HD
#define  DIMENSIONS              2
#define  COMPONENTS              3
#define  GEOMETRY                SPHERICAL
#define  BODY_FORCE              VECTOR
#define  COOLING                 NO
#define  RECONSTRUCTION          PARABOLIC
#define  TIME_STEPPING           RK3
#define  DIMENSIONAL_SPLITTING   NO
#define  NTRACER                 1
#define  USER_DEF_PARAMETERS     6

/* -- physics dependent declarations -- */

#define  EOS                     IDEAL
#define  ENTROPY_SWITCH          NO
#define  THERMAL_CONDUCTION      NO
#define  VISCOSITY               EXPLICIT
#define  ROTATING_FRAME          NO

/* -- user-defined parameters (labels) -- */

#define  RMIN                    0
#define  RMAX                    1
#define  RHO_CUT                 2
#define  BETA                    3
#define  ETA                     4
#define  SCALE_HEIGHT            5

/* [Beg] user-defined constants (do not change this line) */

#define  USE_DIPOLE              NO

/* [End] user-defined constants (do not change this line) */

/* -- supplementary constants (user editable) -- */ 

#define  INITIAL_SMOOTHING   NO
#define  WARNING_MESSAGES    YES
#define  PRINT_TO_FILE       YES
#define  INTERNAL_BOUNDARY   NO
#define  SHOCK_FLATTENING    NO
#define  CHAR_LIMITING       NO
#define  LIMITER             MC_LIM
