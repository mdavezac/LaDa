/* 
 *  This is code sets up and uses an integer array appropriate for maximizing a
 *  function appropriate for an zincblend alloy.
 */


#ifndef FAKE_MPI
#define IAGA_MPI
#endif

#include <stdio.h>
#include <math.h>

#ifdef IAGA_MPI
#include "mpi.h"
#endif

#include "hdf5.h"
#include "pgapack.h"

#include "paramio.h"
#include "iaga.h"

void iaga_FullCheckComp(IAGAContext *IAGActx, int pop);

int iaga_test_evaluptodate(PGAContext *ctx, int pop)
{
    int p;
    PGAInteger PopSize = PGAGetPopSize(ctx);
    PGAIndividual *Indiv;

    for (p = 0; p < PopSize; p++)
    {
	Indiv = PGAGetIndividual(ctx, p, pop);
	if(Indiv->evaluptodate == PGA_TRUE) 
	{
	    if(Indiv->evalfunc != AlloyFunctional( ctx, p, pop)) 
	    {
		printf("IAGA: Individual %d of Pop %d uptodate, but fitness = %f, eval = %f != AlloyFunct = %f\n",p,pop,Indiv->fitness,Indiv->evalfunc,AlloyFunctional( ctx, p, pop) );
		PGAPrintString(ctx,stdout,p,pop);
	    }
	}
    }
    return 0;
}

int iaga_sort_pop(PGAContext *ctx)
{
    int i, j;
    PGAInteger n = PGAGetPopSize(ctx);

    PGASortPop( ctx, PGA_OLDPOP);
    for ( i=0; i < n; i++ )
    {
        j = PGAGetSortedPopIndex( ctx, i );
        PGACopyIndividual ( ctx, j, PGA_OLDPOP, i, PGA_NEWPOP );
    }
    for ( i=0; i < n; i++ )
    {
        PGACopyIndividual ( ctx, i, PGA_NEWPOP, i, PGA_OLDPOP );
    }
	return 0;
}

void copyfile(char *dest, char *source)
{
#define BUFSIZE 1000
    char ptr[BUFSIZE];
    FILE *fd, *fs;
    fs = fopen(source, "r");
    fd = fopen(dest, "w+");
    while (!feof(fs))
    {
	fread(ptr, 1, BUFSIZE, fs);
	fwrite(ptr, 1, BUFSIZE, fd);
    }
    fclose(fs);
    fclose(fd);
}


int IAGA_PGARun(IAGAContext *IAGActx, int restart)
{
    PGAContext *ctx = IAGActx->PGActx;
    int rank;
#ifdef IAGA_MPI
    MPI_Comm  comm  = IAGActx->MPIctx->comm;
#endif
    int Restarted;
    int *TotalCountFinal=&(IAGActx->TotalCountFinal);
    int IAGA_Restart = IAGA_FALSE;
    int IAGA_RestartIter=0;
    char *tempStr;
    
    IAGActx->TotalCount=0;    
/*    printf("IAGA: Initial Mutation Prob =  %f\n", PGAGetMutationProb(ctx));*/    
    rank = PGAGetRank(ctx, comm);
    if( rank == 0 ) 
    {
	PARAMIO_ReadFileLogical(IAGActx->Infilename,"IAGA_Restart", &IAGA_Restart);
	IAGA_Restart = IAGA_Restart || restart;
	if ( IAGA_Restart) 
	{ 
	    IAGA_RestartIter = PGAGetMaxGAIterValue(ctx);
	    PARAMIO_ReadFileInt(IAGActx->Infilename,"IAGA_RestartIter", &IAGA_RestartIter);
	    printf("IAGA: reading IAGA_RestartIter %d\n",IAGA_RestartIter);
	    tempStr = PARAMIO_ReadFileString(IAGActx->Infilename,"IAGA_RestartFilename");
	    if (tempStr)	    
		strcpy(IAGActx->Restartfilename, tempStr);
	    if (strcmp(IAGActx->Outfilename, IAGActx->Restartfilename) != 0)
	    {
		/* copy to outfilename so we get generations before restart iter in our output file */
		copyfile(IAGActx->Outfilename, IAGActx->Restartfilename);
		strcpy(IAGActx->Restartfilename, IAGActx->Outfilename);
	    }
	    hdf5_read_restart(IAGActx,&IAGA_RestartIter, restart); 
	    if (restart)
	    {
		/* default for command line restart is to go MaxGAIter _more_ generations */
		PGASetMaxGAIterValue(ctx, IAGA_RestartIter +  PGAGetMaxGAIterValue(ctx));
	    }
	} 
    }
    *TotalCountFinal=0;
    
    printf("IAGA: MPI Rank %d %d rank %d %d\n",PGAGetRank(ctx,MPI_COMM_WORLD),
	   MPI_COMM_WORLD,PGAGetRank(ctx,PGAGetCommunicator(ctx)),PGAGetCommunicator(ctx));
    
    if( rank == 0 && IAGA_GetMinimizeDistanceFlag(IAGActx) == IAGA_TRUE ) 
    {
	IAGA_MinimizePopulation(IAGActx, PGA_OLDPOP,IAGA_FALSE);
    } 
    PGAEvaluate(ctx, PGA_OLDPOP, AlloyFunctional, comm);
    /* makes sure that evaluations are up-to-date
       with respect to convex hull -- cv may change during evals */
    #ifdef _IAGA_LAMARCK_ 
      if ( convex_hull_has_changed() )
      {
        /* only re-evaluation, no VA, no Convex Hull recomputation */
        set_recalculating(1); 
        PGAEvaluate(ctx, PGA_OLDPOP, AlloyFunctional, comm);
        print_convex_hull();
        /* resets VA and CH computation */
        set_recalculating(0);
      }
    #endif /*  _IAGA_LAMARCK_ */
    if( rank == 0 ) 
    { 
	if ( IAGActx->ConserveSigmaFlag == IAGA_TRUE ) 
	    IAGA_Fitness(IAGActx, ctx, PGA_OLDPOP);
	else 
	    PGAFitness(ctx, PGA_OLDPOP);
    }
    
    if ( rank == 0)
    {
	if (IAGA_Restart == IAGA_FALSE || strcmp(IAGActx->Outfilename, IAGActx->Restartfilename) !=0) 
	    IAGAOutputHDF(IAGActx,  PGA_OLDPOP, IAGA_HDF_CREATE);   
	else 
	    IAGAOutputHDF(IAGActx,  PGA_OLDPOP, IAGA_HDF_OPEN);   
    }
    printf("IAGA: Starting While Loop %d\n",rank);
    
    while(!PGADone(ctx, comm)) 
    {
/*     iaga_test_evaluptodate(ctx,PGA_OLDPOP);  */
	if ( rank == 0 ) 
	{ 
/* 	    iaga_FullCheckComp(IAGActx, PGA_OLDPOP);  */
	    if ( (IAGActx->CSA_Flag  == IAGA_TRUE)
		 && (ctx->ga.ItersOfSame % IAGActx->CSA_Frequency == 0) 
		 && (ctx->ga.ItersOfSame != 0 ) ) 
	    {
		IAGA_CSA( IAGActx, ctx);
	    }
	    
	    Restarted = PGA_FALSE;
	    if ( (ctx->ga.restart == PGA_TRUE) &&
		 (ctx->ga.ItersOfSame % ctx->ga.restartFreq == 0) &&
		 (ctx->ga.ItersOfSame != 0 ) ) 
	    {
		printf("IAGA: Restarting \n");
		ctx->ga.ItersOfSame++;
		Restarted = PGA_TRUE;
		PGARestart(ctx, PGA_OLDPOP, PGA_NEWPOP);
	    } 
	    else 
	    {
		PGASelect(ctx, PGA_OLDPOP);
		if (PGAGetMutationOrCrossoverFlag(ctx)) 
		    PGARunMutationOrCrossover(ctx, PGA_OLDPOP, PGA_NEWPOP);
		else 
		    PGARunMutationAndCrossover(ctx, PGA_OLDPOP, PGA_NEWPOP);
	    }
	}
	MPI_Bcast(&Restarted, 1, MPI_INT, 0, comm);
	
	if( rank == 0 && IAGA_GetMinimizeDistanceFlag(IAGActx) == IAGA_TRUE) 
	{
	    IAGA_MinimizePopulation(IAGActx, PGA_NEWPOP, IAGA_FALSE);
	}

	PGAEvaluate(ctx, PGA_NEWPOP, AlloyFunctional, comm);
        /* makes sure that evaluations are up-to-date
           with respect to convex hull */
        #ifdef _IAGA_LAMARCK_ 
          increment_iterations();
          if ( convex_hull_has_changed() )
          {
            /* only re-evaluation, no VA, no Convex Hull recomputation */
            set_recalculating(1); 
            PGAEvaluate(ctx, PGA_OLDPOP, AlloyFunctional, comm);
            PGAEvaluate(ctx, PGA_NEWPOP, AlloyFunctional, comm);
            print_convex_hull();
            /* resets VA and CH computation */
            set_recalculating(0);
            if ( rank == 0 )
            {
              if ( IAGActx->ConserveSigmaFlag == IAGA_TRUE ) 
		IAGA_Fitness(IAGActx, ctx, PGA_OLDPOP);
              else
	  	PGAFitness(ctx, PGA_OLDPOP);
            }
          }
        #endif /*  _IAGA_LAMARCK_ */
	if( rank == 0 ) 
	{ 
	    if ( IAGActx->ConserveSigmaFlag == IAGA_TRUE ) 
		IAGA_Fitness(IAGActx, ctx, PGA_NEWPOP);
	    else 
		PGAFitness(ctx, PGA_NEWPOP);
	}
	if (!Restarted) 
	{
	    PGAUpdateGeneration(ctx, comm);
	    if ( rank == 0 ) 
		PGAPrintReport (ctx, stdout, PGA_OLDPOP);
	    if ( rank == 0 ) 
	    {
		/*iaga_sort_pop(ctx);*/
		if (ctx->ga.iter > PGAGetMaxGAIterValue(ctx) && IAGActx->IAGA_MinimizeFinalDistanceFlag) 
		{
		    IAGA_MinimizePopulation(IAGActx, PGA_OLDPOP,IAGA_TRUE);
		}
		IAGAOutputHDF(IAGActx,  PGA_OLDPOP, IAGA_HDF_OPEN);   
	    }
	}
	MPI_Reduce(&(IAGActx->TotalCount), TotalCountFinal, 1, MPI_INT, MPI_SUM, 0, comm);
	if (rank == 0)
	{
	    printf("iter = %d, total(unique)configs = %d, totalevals = %d, best = %f\n", ctx->ga.iter, IAGActx->NumConfigsTried,
			*TotalCountFinal, PGAGetEvaluation(ctx, PGAGetBestIndex(ctx, PGA_OLDPOP), PGA_OLDPOP));
	}
    }
    #ifdef _IAGA_LAMARCK_ 
      print_convex_hull(1); /* prints convex hull in xml format */
    #endif /*  _IAGA_LAMARCK_ */
    return 0;
}

void IAGA_exit(int status) {

#ifdef IAGA_MPI
  MPI_Finalize();
#endif
  exit(status);
}

#if 0
 needs to be fixed to be used again
void iaga_FullCheckComp(IAGAContext *IAGActx, int pop)
{
    int p;
    PGAContext *ctx=IAGActx->PGActx;
    PGAInteger PopSize = PGAGetPopSize(ctx);
    SOLUTION *parent;
    int worstdiff = 0;

    for (p = 0; p < PopSize; p++)
    {
	parent = (SOLUTION *)PGAGetIndividual(ctx, p, pop)->chrom;
	worstdiff = IAGA_MAX(IAGA_CheckComp(GET_ALLOY(IAGActx->Probctx)->Cation, ctx, parent,"iaga_PGARun Full check cation"), worstdiff);
    }
    printf("IAGA: Checking composition, largest mismatch = %d\n", worstdiff);
}
#endif
