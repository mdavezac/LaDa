#include <config.h>

#include <stdlib.h>
#include <string.h>

#include "paramio.h"
#include "util.h"
#include "testfunct.h"
#include "ce.h"
#include "writephys.h"
#include "solver.h"

/**********************************
for now, solver set up is also in here. they are related
************************************/

void FillSolverCtxBandgap(SolverContext *ctx);
void FillSolverCtxEffectiveMass(SolverContext *ctx);
void FillSolverCtxStrain(SolverContext *ctx);
void FillSolverCtxSurfaceToVol(SolverContext *ctx);
void FillSolverCtxIsing3d(SolverContext *ctx);
void FillSolverCtxIsing3dFreeCell(SolverContext *ctx);
void FillSolverCtxClusterExpansion(SolverContext *ctx);
void FillSolverCtxSQS(SolverContext *ctx);

#ifdef _IAGA_LAMARCK_
  extern void FillSolverCtxVA(SolverContext *ctx);
#endif // _IAGA_LAMARCK_

int SolverCleanup( SolverContext **ppSolver)
{
	if ((*ppSolver)->fnSolverCleanup)
		 ((*ppSolver)->fnSolverCleanup)(ppSolver);
	free(*ppSolver);
	*ppSolver = NULL;
	return 0;
}

int SolverInit(char *infilename, SolverContext **ppSolver)
{
    *ppSolver = (SolverContext *)malloc (sizeof(SolverContext));
    if (!ppSolver)
		MDOAbort("failed to allocated solver\n");
	strcpy((*ppSolver)->infilename, infilename);
    (*ppSolver)->SolverNum = IAGA_SOLVER_BANDGAP;
    PARAMIO_ReadFileInt(infilename,"IAGA_SolverNumber", &(*ppSolver)->SolverNum);
    (*ppSolver)->WritePhysicsFiles = FALSE;  
    (*ppSolver)->EPMtype = -1;  /*uninitialized */  	
    (*ppSolver)->DoWave = FALSE;  
	PARAMIO_ReadFileLogical(infilename,"IAGA_PESCAN_WritePhysicsFiles", &(*ppSolver)->WritePhysicsFiles);
	PARAMIO_ReadFileLogical(infilename,"IAGA_PESCAN_DoWave", &(*ppSolver)->DoWave);
	
	(*ppSolver)->fnSolverSetup = NULL;
    (*ppSolver)->fnSolverCleanup = NULL;
    (*ppSolver)->fnSolverEvaluate = NULL;
	(*ppSolver)->fnSolverEvaluateAtoms = NULL;
	
    switch ((*ppSolver)->SolverNum)
    {
	case IAGA_SOLVER_BANDGAP:
	    FillSolverCtxBandgap(*ppSolver);
	    printf ("IAGA: initializing solver: bandgap \n");
	    break;

	case IAGA_SOLVER_EFFECTIVE_MASS:
	    FillSolverCtxEffectiveMass(*ppSolver);
	    printf ("IAGA: initializing solver: effective mass \n");
	    break;
	    	    
       case IAGA_SOLVER_STRAIN:
	    FillSolverCtxStrain(*ppSolver);
	    printf ("IAGA: initializing solver: vff strain energy \n");
	    break;
	    	    
	case IAGA_SOLVER_SURFACE_TO_VOL:
	    FillSolverCtxSurfaceToVol(*ppSolver);
	    printf ("IAGA: initializing solver: surface to volume \n");
	    break;
	    	    
	case IAGA_SOLVER_ISING3D:
	    FillSolverCtxIsing3d(*ppSolver);
	    printf ("IAGA: initializing solver: ising 3d \n");
	    break;
	case IAGA_SOLVER_ISING3D_FREECELL:
	    FillSolverCtxIsing3dFreeCell(*ppSolver);
	    printf ("IAGA: initializing solver: ising 3d for non-fixed supercell\n");
	    break;
	case IAGA_SOLVER_CLUSTER_EXPANSION:
	    FillSolverCtxClusterExpansion(*ppSolver);
	    printf ("IAGA: initializing solver: cluster expansion\n");
	    break;
	case IAGA_SOLVER_SQS:
	    FillSolverCtxSQS(*ppSolver);
	    printf ("IAGA: initializing solver: special quasirandom structure (SQS) generation\n");
	    break;
        #ifdef _IAGA_LAMARCK_ 
          case IAGA_SOLVER_VA:
              FillSolverCtxVA(*ppSolver);
              printf ("IAGA: initializing solver: VA\n");
              break;
        #endif // _IAGA_LAMARCK_
	    	    
	default:
		MDOAbort("unknown solver specified in input file");
	    return -1;

    }
    return 0;
}

int SolverSetup(SolverContext *solver)
{
	if (solver->fnSolverSetup)
		return (*solver->fnSolverSetup)(solver);
	else
		return 0;
}

void check_write_phys(SolverContext *solver)
{
    int rank;
    /* perhaps check and rewrite physics input files */
    /* this needs the global lattice set up, which is direct from
       iaga.in and should not be part of the problem setup 
       unless it is in the problemPREsetup */
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    if (rank == 0)
    {
		if (solver->WritePhysicsFiles) 
		{
			printf ("rewriting physics files \n");
			/* this function is not declared, a la iaga_functional.
			we are going to change its name someday */
			PESCAN_WritePhysicsInputFilesEpm(GetGlobalLattice(), solver->EPMtype);
		}
		PESCAN_ConsistencyCheck();
    }
}


int BandgapSolverSetup(SolverContext *solver)
{
	check_write_phys(solver);
    return 0;
}
	
int EffectiveMassSolverSetup(SolverContext *solver)
{
	check_write_phys(solver);
    return 0;
}
	
int StrainSolverSetup(SolverContext *solver)
{
	check_write_phys(solver);
    return 0;
}
	
int SurfaceToVolSolverSetup(SolverContext *solver)
{
	check_write_phys(solver);
    return 0;
}
	
int Ising3dSolverSetup(SolverContext *solver)
{
    return 0;
}
	
int ClusterExpansionSolverSetup(SolverContext *solver)
{
	/* read in clusters */
	char *figfilename = PARAMIO_ReadFileString(solver->infilename, "CE_FigureFile");
	ce_readfigures(figfilename, TRUE);
	free(figfilename);	
	return 0;
}

int SQSSolverSetup(SolverContext *solver)
{
	/* read in clusters */
	char *figfilename = PARAMIO_ReadFileString(solver->infilename, "CE_FigureFile");
	ce_readfigures(figfilename, FALSE);
	free(figfilename);	
	return 0;
}

int ClusterExpansionSolverCleanup(SolverContext *solver)
{
	ce_freeallfigs();
	return 0;
}

void FillSolverCtxClusterExpansion(SolverContext *ctx)
{
    ctx->fnSolverSetup = ClusterExpansionSolverSetup;
    ctx->fnSolverCleanup = ClusterExpansionSolverCleanup;
    ctx->fnSolverEvaluateAtoms = (FnSolverEvaluateAtoms)cluster_expansion_functional;
    ctx->fnSolverEvaluate = NULL;
}

void FillSolverCtxSQS(SolverContext *ctx)
{
    ctx->fnSolverSetup = SQSSolverSetup;
    ctx->fnSolverCleanup = ClusterExpansionSolverCleanup;
    ctx->fnSolverEvaluateAtoms = (FnSolverEvaluateAtoms)sqs_functional;
    ctx->fnSolverEvaluate = NULL;
}

void FillSolverCtxBandgap(SolverContext *ctx)
{
    ctx->fnSolverSetup = BandgapSolverSetup;
    ctx->fnSolverEvaluateAtoms = (FnSolverEvaluateAtoms)bandgap_c_functional;
    ctx->fnSolverEvaluate = NULL;
} 

void FillSolverCtxEffectiveMass(SolverContext *ctx)	
{
    ctx->fnSolverSetup = EffectiveMassSolverSetup;
    ctx->fnSolverEvaluateAtoms = (FnSolverEvaluateAtoms)emass_c_functional;
    ctx->fnSolverEvaluate = NULL;
}

void FillSolverCtxStrain(SolverContext *ctx)	
{
    ctx->fnSolverSetup = StrainSolverSetup;
    ctx->fnSolverEvaluateAtoms = (FnSolverEvaluateAtoms)strain_c_functional;
    ctx->fnSolverEvaluate = NULL;
}

void FillSolverCtxSurfaceToVol(SolverContext *ctx)
{
    ctx->fnSolverSetup = SurfaceToVolSolverSetup;
    /*    ctx->fnSolverEvaluateAtoms = (FnSolverEvaluateAtoms)surface_to_vol_functional;*/
    MDOAbort("surface to volume functional currently disabled.");
    ctx->fnSolverEvaluateAtoms = NULL;
    ctx->fnSolverEvaluate = NULL;
} 

void FillSolverCtxIsing3d(SolverContext *ctx)
{
    ctx->fnSolverSetup = Ising3dSolverSetup;
    ctx->fnSolverEvaluateAtoms = (FnSolverEvaluateAtoms)ising3d_functional;
    ctx->fnSolverEvaluate = NULL;
} 

void FillSolverCtxIsing3dFreeCell(SolverContext *ctx)
{
    ctx->fnSolverSetup = Ising3dSolverSetup;
    ctx->fnSolverEvaluateAtoms = (FnSolverEvaluateAtoms)ising3d_functional_freecell;
    ctx->fnSolverEvaluate = NULL;
} 

void setup_pos_in(AtomicConfig *pAtoms, MPI_Comm *comm_funct)
{
    int mynode_funct, mynode_world;
    MPI_Comm_rank(*comm_funct,&mynode_funct);
    MPI_Comm_rank(MPI_COMM_WORLD,&mynode_world);

    if (mynode_funct == 0)
    {	
        writeposition(pAtoms, mynode_world, FALSE);
    }
}

void bandgap_c_functional(SolverContext *solver, AtomicConfig *pAtoms, double *functional_p, MPI_Comm *comm_funct, int *error, char *tag)
{
	/*  generate position.in file */
	/* then call fortan to do physics */
    IagaFunctData data;
    int mynode_funct, mynode_world;
    MPI_Comm_rank(*comm_funct,&mynode_funct);
    MPI_Comm_rank(MPI_COMM_WORLD,&mynode_world);

    if (mynode_funct == 0)
	setup_pos_in(pAtoms, comm_funct);
    if (mynode_funct == 0 &&  solver->WritePhysicsFiles)
    {
	write_escan_input_c(&pAtoms->Lattice, solver->EPMtype, mynode_world, tag, solver->DoWave);
    }
    data.comm_funct = *comm_funct;
    data.vffmode = 1; /* do whole cycle */

    FC_FUNC_(bandgap_functional, BANDGAP_FUNCTIONAL)(&data);


#if 0  
This (below) has been replace with FC_FUNC macro
#if (defined(WITH_PHYS) || defined (WITH_MPI_PHYS)) 
#ifdef FORTRANUNDERSCORE 
    bandgap_functional_(&data);
#else 
    bandgap_functional(&data);
#endif 
#endif
#endif
    *functional_p = data.functional_value;
    printf("finished bandgap calc\n");
}

void emass_c_functional(SolverContext *solver,AtomicConfig *pAtoms, double *functional_p, MPI_Comm *comm_funct, int *error, char *tag)
{
    setup_pos_in(pAtoms, comm_funct);
    FC_FUNC_(emass_functional, EMASS_FUNCTIONAL)(functional_p, comm_funct, error);

#if 0
#if (defined(WITH_PHYS) || defined (WITH_MPI_PHYS)) 
#ifdef FORTRANUNDERSCORE 
    emass_functional_(functional_p, comm_funct, error);
#else 
	emass_functional(functional_p, comm_funct, error);
#endif 
#endif	
#endif
}

void strain_c_functional(SolverContext *solver,AtomicConfig *pAtoms, double *functional_p, MPI_Comm *comm_funct, int *error, char *tag)
{
  double strain = 0;
    int mynode_funct, mynode_world;
    MPI_Comm_rank(*comm_funct,&mynode_funct);
    MPI_Comm_rank(MPI_COMM_WORLD,&mynode_world);

    if (mynode_funct == 0)
    {
	setup_pos_in(pAtoms, comm_funct);
	FC_FUNC_(relax_vff,RELAX_VFF)(&strain);
#if 0
#if (defined(WITH_PHYS) || defined (WITH_MPI_PHYS)) 
#ifdef FORTRANUNDERSCORE 
	relax_vff_(&strain);
#else 
	relax_vff(&strain);
#endif
#endif
#endif

    }
    *functional_p = strain;
}

/* this function moves the generation of postions and writing of position.in to C instead of fortran.
 The reason we need
to do this is that we want this part of our iaga_functional call to use a structured type which will include
a list of the atoms, and there appears to be no way in fortran to pass such an object */
/* now this is part of an interface such that multiple optimizers can call multiple solvers, none built to
know much (if anything) about eachother */
void SolverCallFunctional(SolverContext *solver, AtomicConfig *patoms, double *egap, MPI_Comm *comm_funct, int *ierror, char *tag)
{
	/* all are done the first (ie through a C fn. that gets AtomicConfig passed to it ) way now */
    if (solver->fnSolverEvaluateAtoms)
    {
		(*solver->fnSolverEvaluateAtoms)(solver, patoms, egap,comm_funct, ierror, tag); 
    }
	else
		MDOAbort("no solver installed");
}


/**********************************************************************************
 end  , rest is scrap
*************************************************************/

#if 0
    else
    {
#if  (defined (WITH_PHYS) || defined (WITH_MPI_PHYS))
		(*solver->fnSolverEvaluate)(egap,comm_funct, ierror); 
/*      iaga_functional_(egap,comm_funct, ierror);  */
#else
    printf ("IAGA: Tried to call iaga_functional, but no physics compiled, calling test functional instead\n");
    iaga_test_funct(patoms, egap, comm_funct, ierror);
#endif
    }


#if  (defined (WITH_PHYS) || defined (WITH_MPI_PHYS))
#ifdef FORTRANUNDERSCORE
    ctx->fnSolverEvaluate = (FnSolverEvaluate)iaga_functional_;
#else
    ctx->fnSolverEvaluate = (FnSolverEvaluate)iaga_functional;
#endif
#endif

#if  (defined (WITH_PHYS) || defined (WITH_MPI_PHYS))
#ifdef FORTRANUNDERSCORE 
    ctx->fnSolverEvaluate = (FnSolverEvaluate)emass_functional_;
#else
    ctx->fnSolverEvaluate = (FnSolverEvaluate)emass_functional;
#endif
#endif

#endif

