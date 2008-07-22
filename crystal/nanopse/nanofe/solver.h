#ifndef SOLVER_INCLUDED
#define SOLVER_INCLUDED

#include "mpi.h"
#include "alloygeom.h"

/* we begin the task of hiding the forward solver from the IAGA mechanics */
enum
{
	IAGA_SOLVER_ISING3D = 0,
        IAGA_SOLVER_BANDGAP = 1,
        IAGA_SOLVER_EFFECTIVE_MASS = 2,
        IAGA_SOLVER_SURFACE_TO_VOL = 3,
        IAGA_SOLVER_STRAIN = 4,
	IAGA_SOLVER_CLUSTER_EXPANSION = 5,
	IAGA_SOLVER_SQS = 6,
	IAGA_SOLVER_ISING3D_FREECELL = 7
        #ifdef _IAGA_LAMARCK_
        , IAGA_SOLVER_VA = 8
        #endif
};

/* for passing to iaga functional */
/* this structure is paralleled by a Fortran common block in iaga_funct/iaga_funct_api.f90 */
typedef struct 
{
  int vffmode;   /*0 = novff, only genpot/escan, 1 = vff + gp/escan, 2 = only vff */
  int ierror;
  double functional_value;
  MPI_Comm comm_funct;   /* handle of MPI communicator  */
} IagaFunctData;

typedef struct 
{
  int vffmode;   /*0 = novff, only genpot/escan, 1 = vff + gp/escan, 2 = only vff */
  int ierror;
  double emass[6];
  MPI_Comm comm_funct;   /* handle of MPI communicator  */
} EmassFunctData;


typedef int (*FnSolverSetup)();
typedef int (*FnSolverCleanup)();
typedef int (*FnSolverEvaluate)(void *, void *, void *);
typedef int (*FnSolverEvaluateAtoms)(void *, void *, void *, void *, void *, void *);

#define INFILECHARS 1000
typedef struct
{
    int SolverNum;
	int EPMtype; /* eg algaas etc. */
    int WritePhysicsFiles;  /* these three should be in a solver specific subpart !*/
	int DoWave;  /* whether to compute real space wave fn. for viewing, etc. eventually majority rep...*/
    char infilename[INFILECHARS];
    FnSolverSetup fnSolverSetup;
    FnSolverCleanup fnSolverCleanup;
    FnSolverEvaluate fnSolverEvaluate;
    FnSolverEvaluateAtoms fnSolverEvaluateAtoms;
} SolverContext;

/* these are fortran functions. the iaga_functional interface now uses IagaFunctData (ie the void * is IagaFunctData *) */
/* the emass one will too the next time we use it ! */
void bandgap_functional(void *);
void bandgap_functional_(void *);
void emass_functional(void *, void *, void *);
void emass_functional_(void *, void *, void *);
void raw_emass_functional(void *);
void raw_emass_functional_(void *);

/* These are C functions which get things organized and call the fortran */
void bandgap_c_functional(SolverContext *solver, AtomicConfig *pAtoms, double *functional_p, MPI_Comm *comm_funct, int *error, char *tag);
void emass_c_functional(SolverContext *solver, AtomicConfig *pAtoms, double *functional_p, MPI_Comm *comm_funct, int *error, char *tag);
void strain_c_functional(SolverContext *solver, AtomicConfig *pAtoms, double *functional_p, MPI_Comm *comm_funct, int *error, char *tag);
void surface_to_vol_functional(SolverContext *solver, AtomicConfig *pAtoms, double *functional_p, MPI_Comm *comm_funct, int *error, char *tag);
/*void surface_to_vol_functional(void *, void *, void *, void *);*/


int SolverInit(char *infilename, SolverContext **ppSolver);
int SolverSetup( SolverContext *solver);
int SolverCleanup( SolverContext **ppSolver);
void SolverCallFunctional(SolverContext *solver, AtomicConfig *patoms, double *egap, MPI_Comm *comm_funct, int *ierror, char *tag);




#endif
