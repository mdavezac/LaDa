#ifndef FAKE_MPI
#define IAGA_MPI
#endif

#ifndef UTIL_INCLUDED
#include "util.h"
#endif

#ifndef ALLOYGEOM_INCLUDED
#include "alloygeom.h"
#endif

#ifndef SOLVER_INCLUDED
#include "solver.h"
#endif

#ifndef TESTFUNCT_INCLUDED
#include "testfunct.h"
#endif

#ifndef ALLOYFUNCTIONS_INCLUDED
#include "alloy_functions.h"
#endif

#ifndef SOLUTION_INCLUDED
#include "solution.h"
#endif

#ifndef PROBLEMS_INCLUDED
#include "problems.h"
#endif

#define IAGA_FITNESS_FUNCTION 4 /* Use IAGA fitness function */

#define IAGA_HDF_CREATE 1  /* Create HDF5 File */
#define IAGA_HDF_OPENED 2  /* Use Already Created and Opened HDF5 File*/
#define IAGA_HDF_OPEN 3    /* Open Existing (Created) HDF5 File */

#define IAGA_RUN_PHYS_FUNCT 1
#define IAGA_RUN_TEST_FUNCT 2
#define IAGA_FINALIZE 3

#define IAGA_MAXPATHSIZE 512


/* separation has been introduced:  between (1) thing we evolve, an array of integers, and
   (2) atomic configuration, a description/prescription of/for the exact positions of every
   atom in the material. Also between Alloy / Ion (ie the C structures) and IAGA itself. see alloy_functions.h.
   see alloygeom.h for position/spatial stuff. */

#define GROUP_PGA 1
#define GROUP_FUNCTIONAL_ONLY 2
typedef struct 
{
#ifdef IAGA_MPI
  MPI_Comm comm_pga;
  MPI_Comm comm_funct;
  MPI_Comm comm_iaga;
  MPI_Comm comm;
/* Should Disable: Shouldn't be needed */
  MPI_Group group;
  int group_pga;
/*Number of processor to use for each functional evaluation */
  int size_funct;
#endif
/* Not Used? */
  int enabled;
} MPIContext;

#if 0
 /* some function types */
typedef int (*FnInitStringType)(PGAContext *ctx, int p, int pop);
typedef int (*FnMutateType)(PGAContext *ctx, int p, int pop, double mr);
typedef int (*FnCrosseverType)(PGAContext *ctx, int p1, int p2, int p_pop, int c1, int c2, int c_pop);
typedef int (*FnGenAtomicConfig)(SOLUTION *string, int len, AtomicConfig *patoms);
typedef int (*FnProblemPreSetup)();
typedef int (*FnProblemSetup)();
typedef int (*FnProblemCleanup)();


typedef struct
{
    int ProbNum;
	int StringStorageSize;  /* default to alloy->total */
    FnInitStringType fnInitString;
    FnMutateType fnMutate;
    FnCrosseverType fnCrossover;
    FnGenAtomicConfig fnGenAtomicConfig;
    FnProblemPreSetup fnPreSetup;
    FnProblemCleanup fnCleanup;
    FnProblemSetup fnSetup;
    void *ProbData; /* will be cast to appropriate type for each problem */
} ProblemContext;
#endif

/* the four main functions for every problem to implement a version of */
/* the above types are to match these. These are in iaga.c and are just stubs
which call the above problem specific (but unkown to iaga) functions */
/* these are also callback fns installed in PGA. see below for some more */
int AlloyMutation(PGAContext *ctx, int p, int pop, double mr);
void AlloyCrossover(PGAContext *ctx, int p1, int p2, int p_pop, int c1, int c2, int c_pop);
void AlloyInitString(PGAContext *ctx, int p, int pop);
int GenAtomicConfig(SOLUTION *string, int len, AtomicConfig *patoms);

typedef struct ConfigList
{
	AtomicConfig *pConfig;
	double FunctionalValue;
	struct ConfigList *pNext;
  #ifdef _IAGA_LAMARCK_
	double concentration;
  #endif // _IAGA_LAMARCK_
} ConfigList;

typedef struct 
{
/*    AlloyInfo *alloy; */
    MPIContext *MPIctx;
    PGAContext *PGActx;
    ProblemContext *Probctx; 
    SolverContext *Solverctx;
	
	ConfigList *pAllConfigs;
	int NumConfigsTried;
 
    int IAGA_Phys;
    int IAGA_PESCAN_WritePhysicsFiles;
    int TotalCountFinal;
    int TotalCount;
    int IAGA_Target;
    float IAGA_TargetValue;
    int IAGA_PGARun;
    int CSA_Flag;
    int CSA_Frequency;
    float CSA_MutationReduction;
    float CSA_UniformCrossoverReduction;
    float CSA_CrossoverReduction;
    int IAGA_MinimizeDistanceFlag;
    int IAGA_MinimizeFinalDistanceFlag;
    int IAGA_MinimizeDistanceInCrossoverFlag;
    int ConserveSigmaFlag;
    float ConserveSigmaValue;
    int IAGA_PrintNumIndividualsValue; /* _when_ we print (which is kept track of by PGA), how _many_ to print */
    int IAGA_RestartFrequency; /* how often to write _everything_, thus enabling restart from this generation */
    char Infilename[IAGA_MAXPATHSIZE];  /* input file name (e.g. iaga.in) */
    char Outfilename[IAGA_MAXPATHSIZE]; /* output file name (e.g. iaga.h5) */
    char Restartfilename[IAGA_MAXPATHSIZE]; /* restart file name (e.g. last created iaga.h5 file, but doesn't need to be the same as the output file) */
    int IAGA_Test;    /* for test cases, whether to do 1d or 3d Ising model */
    float IAGA_Test_IsingCC;  /* Ising coupling constant */
    float IAGA_Test_IsingCC2;  /* Ising coupling constant for 2nd neighbor*/
    float IAGA_Test_IsingRadius;  /* Ising cutoff radius for interaction */
    int IAGA_ThresholdStop;
    float IAGA_ThresholdStopVal;
} IAGAContext;

/* initialization/shutdown: */
void IAGA_SetContext( IAGAContext *ctx );
IAGAContext *IAGA_GetContext( );
PGAContext *IAGA_GetPGAContext();
int IAGA_MPI_Setup( int *argc, char ***argv, MPIContext *info);
int IAGA_MPI_Init( int *argc, char ***argv, MPIContext *info );
int IAGA_PGASetup(int *argc, char **argv, IAGAContext *IAGActx);
void IAGA_PGADestroy (PGAContext *ctx);
/*int IAGA_AlloySetup(int *argc, char **argv, IAGAContext *IAGActx);*/
int IAGA_ParamSetup(int *argc, char **argv, IAGAContext *IAGActx);
int IAGA_ProblemPreSetup(IAGAContext *iagactx);
int IAGA_ProblemSetup(IAGAContext *iagactx);
int IAGA_ProblemCleanup(IAGAContext *iagactx);
int IAGA_Finalize( int status );

/* more callbacks from PGA,  but not re-implemented by each problem */
int AlloyCheckDuplicate(PGAContext *ctx, int p1, int pop1, int p2, int pop2);
int IAGACheckStopping(PGAContext *ctx);
void AlloyCreateString (PGAContext *ctx, int p, int pop, int InitFlag);
void AlloyCopyString (PGAContext *ctx, int p1, int pop1, int p2, int pop2);
void AlloyPrintString ( PGAContext *ctx, FILE *fp, int p, int pop);
MPI_Datatype AlloyBuildDatatype(PGAContext *ctx, int p, int pop);

/*void IAGA_WritePhysicsInputFiles(VEC3 cell);*/

/* this is the callback to evaluate an individual. it calls problem dependent GenAtomicConfig,
 then it calls iaga_functional */
double AlloyFunctional(PGAContext *ctx, int p, int pop);
double AlloyFunctionalNoPGA(AtomicConfig *pAtoms, char *tag);

void  IAGAOutputHDF(IAGAContext *IAGActx, int pop, int IAGA_HDF_Create_Open);
int hdf5_read_restart( IAGAContext *IAGActx, int *Iter, int restart_from_end);

herr_t IAGA_Write_Attribute_Float(hid_t group_id, const char* dataset_name, int dataset_num, const char* attribute_name, float *attribute_value);
herr_t IAGA_Write_Dataset(hid_t group_id, const char* dataset_name, int dataset_num, int *dataset, int dataset_size);

int IAGA_CSA(IAGAContext *IAGActx, PGAContext *PGActx);
int IAGA_MinimizeDistance(IAGAContext *IAGActx, PGAContext *PGActx, int p0, int p1, int pop);
int IAGA_CheckDuplicateViaDistance(IAGAContext *IAGActx, PGAContext *PGActx, int p1, int pop1, int p2, int pop2);
void IAGA_Fitness ( IAGAContext *IAGActx, PGAContext *ctx, int popindex );
void IAGA_FitnessFunction(IAGAContext *IAGActx, PGAContext *ctx, PGAIndividual *pop);
void IAGA_exit(int status);
int IAGA_CheckComp(IonInfo *Ion, PGAContext *ctx, SOLUTION *p, const char *string);
/* in AlloyFuncional.c */
int IAGA_RunFunctional( MPIContext *info );
/* in iaga_std.c */
int iaga_restdio( int *id);
/* in AlloySum.c*/
double AlloySum(PGAContext *ctx, int p, int pop);
 /* in AlloyInitString.c */
void AlloyInitIon( IonInfo *Ion, PGAContext *ctx, int p, int pop);
void AlloyInitIonString(int *c, IonInfo *Ion);
/* declarations of other alloy fns, those not involving PGA, moved to alloy_functions.h*/
/* in iaga_PGARun.c*/
int IAGA_PGARun(IAGAContext *IAGActx, int restart);
/* in AlloyMutation.c */
int IonMutation( IonInfo *Ion, PGAContext *ctx, int p, int pop, double mr);
/* in AlloyCrossover.c */
int IonCrossoverConserve( IonInfo *Ion, PGAContext *ctx, int *parent1, int *parent2, int *child);
/* in iaga_Distance.c */
int IAGA_GetMinimizeDistanceInCrossoverFlag(IAGAContext *ctx);
int IAGA_GetMinimizeDistanceFlag(IAGAContext *ctx);
int IAGA_MinimizePopulation(IAGAContext *IAGActx,int pop, int unconditional);
/* in AlloyConserve.c */
int IonNotConserveValue( IonInfo *Ion, PGAContext *ctx, SOLUTION *AlloyArray);

/***************** hiding PGA ************/
void int2PGAcpy(PGAInteger *dst, int *src, int len);
int *PGA2int(PGAInteger *src, int len);


