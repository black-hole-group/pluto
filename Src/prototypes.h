/* ---------------------------------------------------------------------
                    General function prototypes 
   --------------------------------------------------------------------- */

void   AdvectFlux (const State_1D *, int, int, Grid *);
void   Analysis (const Data *, Grid *);
#if EOS == BAROTROPIC
 double BAROTROPIC_PR(double);
#endif

void BodyForceVectorGet (double **v, double **g,
                         double *, double *, double *, int beg, int end, 
                         Grid *grid);
void BodyForcePotentialGet (double **v, double **gphi,
                            double *phi_p, double *phi_c,
                            int beg, int end, Grid *grid);


double BodyForcePotential(double, double, double);
void   BodyForceVector(double *, double *, double, double, double);

void  ChangeDumpVar ();
void  CheckPrimStates (double **, double **, double **, int, int);
int   CheckNaN (double **, int, int, int);
unsigned char CheckZone (int z, int bit);
void  ComputeUserVar (const Data *, Grid *);
float ***Convert_dbl2flt (double ***, double, int);

void ConsToPrim3D(Data_Arr, Data_Arr, Grid *);
void PrimToCons3D(Data_Arr, Data_Arr, Grid *);
void CreateImage (char *);

void ComputeEntropy      (const Data *, Grid *);
void EntropySwitch       (const Data *, Grid *);
void EntropyOhmicHeating (const Data *, Data_Arr, double, Grid *);

#ifdef FINITE_DIFFERENCE  
 Riemann_Solver FD_Flux;
 Reconstruct MP5_Reconstruct, PPM_Reconstruct, LIMO3_Reconstruct,
             WENOZ_Reconstruct, WENO3_Reconstruct;
 void FD_GetMaxEigenvalues (const Data *d, State_1D *state, Grid *grid);
#endif

void FindShock (const Data *, Grid *);
void FlagShock (const Data *, Grid *);
void FlagReset (const Data *d);
void Flatten (const State_1D *, int, int, Grid *);
void FreeGrid (Grid *);

void   GetCGSUnits (double *u);
Image  *GetImage (char *);
double *GetInverse_dl (const Grid *);
int    GetNghost (Input *);
char   *GetOutputDir();
double ***GetUserVar (char *);

int LocateIndex(double *, int, int, double);

void Init (double *, double, double, double);
void Initialize(int argc, char *argv[], Data *, Input *, Grid *, Cmd_Line *);

void InternalBoundaryReset (const State_1D *, Time_Step *, int, int, Grid *);

void InputDataFree (void);
void InputDataInterpolate (double *, double, double, double);
void InputDataRead (char *, char *);
void InputDataSet (char *, int *);

int  IsLittleEndian (void);

void MakeState (State_1D *);
void MakeGeometry (Grid *);
double Median (double a, double b, double c);

void   ParabolicFlux(Data_Arr, Data_Arr J, double ***, const State_1D *,
                     double **, int, int, Grid *);
double ParabolicRHS (const Data *, Data_Arr, double, Grid *);
void   ParseCmdLineArgs (int, char *argv[], char *, Cmd_Line *);

int    ParamFileRead    (char *);
char  *ParamFileGet     (const char *, int );
int    ParamExist       (const char *);
int    ParamFileHasBoth (const char *, const char *);

void   PrimToChar (double **, double *, double *); 

void ResetState (const Data *, State_1D *, Grid *);
void RightHandSide (const State_1D *, Time_Step *, int, int, double, Grid *);
void RKC (const Data *d, Time_Step *, Grid *);

void SetColorMap (unsigned char *, unsigned char *, unsigned char *, char *);
void SetDefaultVarNames(Output *);
int  SetDumpVar (char *, int, int);
void SetIndexes (Index *indx, Grid *grid);
int  SetLogFile(char *, Cmd_Line *);
void SetOutput (Data *d, Input *input);
void SetRBox(RBox *, RBox *, RBox *, RBox *);
Riemann_Solver *SetSolver (const char *);
int  Setup (Input *, Cmd_Line *, char *);
void SetGrid (struct INPUT *INI, Grid *);
void SetJetDomain   (const Data *, int, int, Grid *);
void SetOutputDir(char *);
void SplitSource (const Data *, double, Time_Step *, Grid *);
void Startup (Data *, Grid *);
void States (const State_1D *, int, int, Grid *);


void UnsetJetDomain (const Data *, int, Grid *);

void VectorPotentialDiff (double *, int, int, int, Grid *);
void VISC_FLUX (const State_1D *, int beg, int end, Grid *);

/* ---------------------------------------------------------------------
            Prototyping for memory allocation functions
   --------------------------------------------------------------------- */

char    *Array1D (int, size_t);
char   **Array2D (int, int, size_t);
char  ***Array3D (int, int, int, size_t);
char ****Array4D (int, int, int, int, size_t);

double ***ArrayBox(long int, long int, long int, long int, long int, long int);

double ***ArrayBoxMap (int, int, int, int, int, int, double *);
double ***ArrayMap (int, int, int, double *);

unsigned char ***ArrayCharMap (int, int, int, unsigned char *);

void FreeArray1D (void *);
void FreeArray2D (void **);
void FreeArray3D (void ***);
void FreeArray4D (void ****);

void FreeArrayBox(double ***, long, long, long);

void FreeArrayBoxMap (double ***, int, int, int, int, int, int);
void FreeArrayMap (double ***);

void FreeArrayCharMap(unsigned char ***);

#define ARRAY_1D(nx,type)          (type    *)Array1D(nx,sizeof(type))
#define ARRAY_2D(nx,ny,type)       (type   **)Array2D(nx,ny,sizeof(type))
#define ARRAY_3D(nx,ny,nz,type)    (type  ***)Array3D(nx,ny,nz,sizeof(type))
#define ARRAY_4D(nx,ny,nz,nv,type) (type ****)Array4D(nx,ny,nz,nv,sizeof(type))

/* ---------------------------------------------------------------------
            Prototyping for time-stepping functions
   --------------------------------------------------------------------- */

void CharTracingStep(const State_1D *, int, int, Grid *);
void HancockStep    (const State_1D *, int, int, Grid *);
void STS (const Data *d, Time_Step *, Grid *);
int  UpdateSolution (const Data *, Riemann_Solver *, Time_Step *, Grid *);
void UpdateStage(const Data *, Data_Arr, double **, Riemann_Solver *,
                 double, Time_Step *, Grid *);

/* ---------------------------------------------------------------------
            Prototyping for standard output/debugging
   --------------------------------------------------------------------- */

void print  (const char *fmt, ...);
void print1 (const char *fmt, ...);
void Show (double **, int);
void ShowMatrix(double **, int n, double);
void ShowVector (double *, int n);
void ShowConfig(int, char *a[], char *);
void ShowDomainDecomposition (int, Grid *);
void ShowUnits ();
void Trace (double);
void Where (int, Grid *);

/* ---------------------------------------------------------------------
            Prototyping for Boundary condition functions
   --------------------------------------------------------------------- */

void Boundary    (const Data *, int, Grid *);
void FlipSign       (int, int, int *);
void OutflowBound   (double ***, RBox *, int, Grid *);
void PeriodicBound  (double ***, RBox *, int);
void ReflectiveBound(double ***, RBox *, int, int);
void UserDefBoundary (const Data *, RBox *, int,  Grid *); 

/* ---------------------------------------------------------------------
              Prototyping for I/O functions
   --------------------------------------------------------------------- */

#ifdef USE_ASYNC_IO
void Async_EndWriteData (Input *ini);
void Async_BegWriteData (const Data *d, Output *output, Grid *grid);
#endif
int   CloseBinaryFile (FILE *, int);
FILE *OpenBinaryFile  (char *, int, char *);
void ReadBinaryArray (void *, size_t, int, FILE *, int, int);
void ReadHDF5 (Output *output, Grid *grid);

void Restart (Input *, int, int, Grid *);
void RestartDump (Input *);
void RestartGet  (Input *, int, int, int);

void SwapEndian (void *, const int); 

void WriteData (const Data *, Output *, Grid *);
void WriteBinaryArray (void *, size_t, int, FILE *, int);
void WriteHDF5        (Output *output, Grid *grid);
void WriteVTK_Header (FILE *, Grid *);
void WriteVTK_Vector (FILE *, Data_Arr, double, char *, Grid *);
void WriteVTK_Scalar (FILE *, double ***, double, char *, Grid *);
void WriteTabArray (Output *, char *, Grid *);
void WritePPM (double ***, char *, char *, Grid *);
void WritePNG (double ***, char *, char *, Grid *);

/* --------------------------------------------------------------------- 
          Prototype for cooling functions
   --------------------------------------------------------------------- */
   
#if COOLING != NO
 void Numerical_Jacobian (double *v, double **J);
 void Jacobian (double *v, double *rhs, double **dfdy);
 void CoolingSource (const Data *, double, Time_Step *, Grid *);
 #if COOLING == POWER_LAW
  void  PowerLawCooling (Data_Arr, double, Time_Step *, Grid *);
 #endif
 /* move the following elsewhere ?  */
/*
double SolveODE_CK45  (double *, double *, double *, double, double);
double SolveODE_RKF23 (double *, double *, double *, double, double);
double SolveODE_RKF12 (double *, double *, double *, double, double);
*/
 double SolveODE_CK45  (double *, double *, double *, double, double, intList *);
 double SolveODE_RKF23 (double *, double *, double *, double, double, intList *);
 double SolveODE_RKF12 (double *, double *, double *, double, double, intList *);

 double SolveODE_ROS34 (double *, double *, double *, double, double);
 double SolveODE_RK4   (double *, double *, double *, double, intList *);
 double SolveODE_RK2   (double *, double *, double *, double);
 double SolveODE_Euler (double *, double *, double *, double);
#endif

/* ----------------------------------------------
           functions in tools.c 
   ---------------------------------------------- */

void PlutoError (int, char *);


double Length_1 (int i, int j, int k, Grid *);
double Length_2 (int i, int j, int k, Grid *);
double Length_3 (int i, int j, int k, Grid *);

#if UPDATE_VECTOR_POTENTIAL == YES
 void VectorPotentialUpdate (const Data *d, const void *vp, 
                             const State_1D *state, const Grid *grid);
#endif

/* --------------- EOS ------------------------- */

void Enthalpy    (double **, double *, int, int);
void Entropy     (double **, double *, int, int);
void SoundSpeed2 (double **, double *, double *, int, int,  int, Grid *);

double GetEntropy (double x);


/* --- New stuffs -- */

void   WriteAsciiFile (char *fname, double *q, int nvar);
