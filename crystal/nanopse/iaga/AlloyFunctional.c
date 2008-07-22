/* 
 *  This is code sets up and uses an integer array appropriate for maximizing a
 *  function appropriate for an zincblend alloy.
 */

#ifndef FAKE_MPI
#define IAGA_MPI
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifdef IAGA_MPI
#include "mpi.h"
#endif

#include "hdf5.h"
#include "pgapack.h"
#include "paramio.h"
#include "iaga.h"
#include "testfunct.h"

#ifdef _IAGA_LAMARCK_
  extern double evaluate_convex_hull( double concentration );
  extern double get_config_concentration ( AtomicConfig *pAtoms );
#endif // _IAGA_LAMARCK_

/* this is a C functio. It writes the position.in file, then calls the physics (in bandgap_pe.f),
   iaga_functional.  If there is no physics compiled, it calls the test functional  */
void iaga_test_funct(AtomicConfig *pAtoms, double *functional, MPI_Comm *comm_funct, int *error);

int IAGA_RunFunctional( MPIContext *info )
{
    IAGAContext *GL_IAGActx = IAGA_GetContext( );
    SolverContext *solver = GL_IAGActx->Solverctx;
	int EvalOrFinalize=IAGA_RUN_PHYS_FUNCT;
	double FunctionalEvaluation;
	int error; 
	AtomicConfig atoms;

  atoms.NumAtoms = -1;
  while (EvalOrFinalize != IAGA_FINALIZE ) 
  {
    MPI_Bcast(&EvalOrFinalize, 1, MPI_INT, 0, info->comm_funct);

    if ( EvalOrFinalize == IAGA_RUN_PHYS_FUNCT ) 
    {
		printf ("sub-process running physics:start\n");
		SolverCallFunctional(solver, &atoms, &FunctionalEvaluation, &(info->comm_funct), &error, "");
		printf ("sub-process running physics:end\n");
    }

    if ( EvalOrFinalize == IAGA_RUN_TEST_FUNCT ) 
    {
		iaga_test_funct(NULL, &FunctionalEvaluation, &(info->comm_funct), &error);
    }
  }

  IAGA_Finalize(0);
  return 0;
}

double AlloyFunctionalNoPGA(AtomicConfig *pAtoms, char *tag)
{
    IAGAContext *GL_IAGActx = IAGA_GetContext( );
    SolverContext *solver = GL_IAGActx->Solverctx;
   int EvalOrFinalize;
    int error;
    double fctl_val;
 
 /* the cases should go away as the test and real functionals share a common interface */
 
    if ( GL_IAGActx->IAGA_Phys == IAGA_TRUE ||  GL_IAGActx->Solverctx->SolverNum == 3) 
    {
		EvalOrFinalize = IAGA_RUN_PHYS_FUNCT;
		MPI_Bcast(&EvalOrFinalize, 1, MPI_INT, 0, GL_IAGActx->MPIctx->comm_funct);
/*    iaga_functional_(alloyarray, &alloyarraylen, &BandGap, &(GL_IAGActx->MPIctx->comm_funct), &error);*/
		SolverCallFunctional(solver, pAtoms, &fctl_val, &(GL_IAGActx->MPIctx->comm_funct), &error, tag);
    } 
    else 
    {
		EvalOrFinalize = IAGA_RUN_TEST_FUNCT;
		MPI_Bcast(&EvalOrFinalize, 1, MPI_INT, 0, GL_IAGActx->MPIctx->comm_funct);
		iaga_test_funct(pAtoms, &fctl_val, &(GL_IAGActx->MPIctx->comm_funct), &error);
/*    iaga_test_funct(patoms, functional_p, comm_funct, error);*/
    }
    
    switch (GL_IAGActx->IAGA_Target) 
    {
		case OPTIMIZER_min: fctl_val = fctl_val; break;
		case OPTIMIZER_max: fctl_val = fctl_val; break;
		case OPTIMIZER_target: fctl_val = - fabs(GL_IAGActx->IAGA_TargetValue - fctl_val); break;
    }
    
/*      printf ("returning %f from AlloyFunctional\n", fctl_val);  */
    GL_IAGActx->TotalCount ++;
    return fctl_val;
}

int equal_config(AtomicConfig *a1, AtomicConfig *a2)
{
	int i;
	
	if (a1->NumAtoms !=a2->NumAtoms)
		return FALSE;
	for (i=0; i<min(a1->NumAtoms, a2->NumAtoms); i++)
		if (a1->AtomList[i] != a2->AtomList[i])
			return FALSE;
	return TRUE;
}


#define MAX_CONFIG_SEARCHES 1000
ConfigList *SearchAllConfigs(ConfigList *pList, AtomicConfig *atoms)
{
	ConfigList *pThis = pList;
	int nsearched = 0;
	
	while (pThis && nsearched++ < MAX_CONFIG_SEARCHES)
	{
		if (equal_config(pThis->pConfig, atoms))
			return pThis;
		pThis = pThis->pNext;
	}
	return NULL;
}

void AddToAllConfigs(ConfigList **pList, AtomicConfig *atoms, double fctl_val)
{
	ConfigList *pNew = abort_null((ConfigList *)malloc(sizeof(ConfigList)), "mem failed adding to configs list");

	pNew->pConfig = atoms;
	pNew->FunctionalValue = fctl_val;
	if (!*pList)
	{
		*pList = pNew;
		pNew->pNext = NULL;
	}
	else
	{
		/* insert at head of list */
		pNew->pNext = *pList;
		*pList = pNew;
	}
}

double AlloyFunctional( PGAContext *ctx, int p, int pop) 
{
    IAGAContext *GL_IAGActx = IAGA_GetContext( );
    SOLUTION *soln;
    int alloyarraylen;
    AtomicConfig *atoms = abort_null((AtomicConfig*)malloc(sizeof(AtomicConfig)), "mem fail in AlloyFunctional");
    ConfigList *pList = NULL;
    double fctl_val;
    #ifdef _IAGA_LAMARCK_
      double concentration;
    #endif // _IAGA_LAMARCK_
    char tag[100];
    
    /* this "tag" notion is so we can save every wave function ,etc.  currently not used, so I will disable it because
       it is not working right for case when user DOES NOT want physics files written by us */
    /*    sprintf(tag, "g%di%d", ctx->ga.iter, p); */
    strcpy(tag, ""); 
    soln = (SOLUTION *)PGAGetIndividual(ctx, p, pop)->chrom;
/*     dump_soln(soln); */
    if (soln->atoms[0] == 0)
	 printf("soln not looking good in AlloyFunctional\n");
/*      else */
/*        printf("soln ok\n"); */

    alloyarraylen = (*GL_IAGActx->Probctx->fnGetStorageSize)(GL_IAGActx->Probctx);
/*    MDOAssert(alloyarraylen == soln->natoms, "length mismatch in AlloyFunctional!"); */
    atoms->Lattice.Basis = NULL;
    GenAtomicConfig(soln, alloyarraylen, atoms);

	pList = SearchAllConfigs(GL_IAGActx->pAllConfigs, atoms);
	/* turned off for free basis search */
	pList = NULL;
	if (pList)
	{
		fctl_val = pList->FunctionalValue;
                #ifdef _IAGA_LAMARCK_
                  concentration = pList->concentration;
                #endif // _IAGA_LAMARCK_
		printf("this config already tested, fctl = %f\n", fctl_val);
		free_AtomicConfig(atoms);
		free(atoms);
	}
	else
	{
		fctl_val =  AlloyFunctionalNoPGA(atoms, tag);
		printf("new (non recent) config, fctl = %f\n", fctl_val);
		AddToAllConfigs(&GL_IAGActx->pAllConfigs, atoms, fctl_val);
                #ifdef _IAGA_LAMARCK_
                  concentration = get_config_concentration(atoms);
                  GL_IAGActx->pAllConfigs->concentration = concentration;
                #endif // _IAGA_LAMARCK_
		GL_IAGActx->NumConfigsTried++;
	}
    #ifdef _IAGA_LAMARCK_
      return (fctl_val - evaluate_convex_hull( concentration ));
    #else
      return (fctl_val);
    #endif // _IAGA_LAMARCK_
}


void iaga_test_funct(AtomicConfig *pAtoms, double *functional_p, MPI_Comm *comm_funct, int *error)
{
    IAGAContext *GL_IAGActx = IAGA_GetContext( );

    test_funct_noiaga(pAtoms, functional_p, comm_funct, error, GL_IAGActx->IAGA_Test,
		      GL_IAGActx->IAGA_Test_IsingCC, GL_IAGActx->IAGA_Test_IsingCC2);
}


