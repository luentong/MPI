/*
 * Simple parallel program to test for percolation of a cluster
 */
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <time.h>
#include <stdbool.h>
#include <string.h>
#include "arralloc.h"

#include "percolate.h"


int compareFile(FILE * fPtr1, FILE * fPtr2, int * line, int * col)
{
    char ch1, ch2;

    *line = 1;
    *col  = 0;

    do
    {
        // Input character from both files
        ch1 = fgetc(fPtr1);
        ch2 = fgetc(fPtr2);

        // Increment line
        if (ch1 == '\n')
        {
            *line += 1;
            *col = 0;
        }

        // If characters are not same then return -1
        if (ch1 != ch2)
            return -1;

        *col  += 1;

    } while (ch1 != EOF && ch2 != EOF);


    /* If both files have reached end */
    if (ch1 == EOF && ch2 == EOF)
        return 0;
    else
        return -1;
}
/*
 * Simple parallel program to test for percolation of a cluster.
 */

int main(int argc, char *argv[])
{
  /*
   *  Define the main arrays for the simulation
   */
  clock_t t;
  t = clock();
  int L;
  L = atoi(argv[2]);
  int NPROC;
  NPROC = atoi(argv[3]);
  int M = L/NPROC;
  int N = L;
  int **old, **new;
  //int old[M+2][N+2], new[M+2][N+2];

  /*
   *  Additional array WITHOUT halos for initialisation and IO. This
   *  is of size LxL because, even in our parallel program, we do
   *  these two steps in serial
   */
  int **map;
  //int map[L][L];
  int **maptmp;
  //int maptmp[L][L];

  /*
   *  Array to store local part of map
   */
  int **smallmap;
  //int smallmap[M][N];

  /*
   *  Variables that define the simulation
   */

  int seed;
  double rho;

  /*
   *  Local variables
   */

  int i, j, nhole, step, maxstep, oldval, newval;
  int nchangelocal, nchange, printfreq;
  int itop, ibot, perc;
  double r;

  /*
   *  MPI variables
   */

  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Status status;

  int size, rank, prev, next;
  int tag = 1;

  MPI_Init(&argc, &argv);

  MPI_Comm_size(comm, &size);
  MPI_Comm_rank(comm, &rank);

  next = rank + 1;
  prev = rank - 1;

  /*
   * Non-periodic boundary conditions
   *
   * Note that the special rank of MPI_PROC_NULL is a "black hole" for
   * communications. Using this value for processes off the edges of the
   * image means there is no additional logic needed to ensure processes
   * at the edges do not attempt to send to or receive from invalid
   * ranks (i.e. rank = -1 and rank = NPROC).
   *
   * Proper solution would compute neighbours with a Cartesian topology
   * and MPI_Cart_shift, where MPI_PROC_NULL is assigned automatically.
   */

  if (next >= size)
    {
      next = MPI_PROC_NULL;
    }

  if (prev < 0)
    {
      prev = MPI_PROC_NULL;
    }

  if (NPROC != size)
    {
      if (rank == 0)
	{
	  printf("percolate: ERROR, NPROC = %d but running on %d\n",
		 NPROC, size);
	}

      // MPI_Finalize();
      // return 0;
  }

  if (argc != 4)
    {
      if (rank == 0)
	{
	  printf("Usage: percolate <seed> <size> <NPROC>\n");
	}

      MPI_Finalize();
      return 0;
    }
  old = (int **) arralloc(sizeof(int), 2, M+2, N+2);
  new = (int **) arralloc(sizeof(int), 2, M+2, N+2);
  map = (int **) arralloc(sizeof(int), 2, L, L);
  maptmp = (int **) arralloc(sizeof(int), 2, L, L);
  smallmap = (int **) arralloc(sizeof(int), 2, M, N);
  if (NULL == old || NULL == new || NULL == map || NULL == maptmp || NULL == smallmap)
{
  printf("percolate: array allocation failed\n");
  return 1;
}

  /*
   *  Update for a fixed number of steps, periodically report progress
   */

  maxstep = 5*L;
  if (L < 100){
    maxstep = 5 * 100;
  }
  printfreq = 100;


  if (rank == 0)
    {
      printf("percolate: running on %d process(es)\n", NPROC);

      /*
       *  Set most important value: the rock density rho (between 0 and 1)
       */

      rho = 0.4064;

      /*
       *  Set the randum number seed and initialise the generator
       */

      seed = atoi(argv[1]);

      printf("percolate: L = %d, rho = %f, seed = %d, maxstep = %d\n",
	     L, rho, seed, maxstep);

      rinit(seed);

      /*
       *  Initialise map with density rho. Zero indicates rock, a positive
       *  value indicates a hole. For the algorithm to work, all the holes
       *  must be initialised with a unique integer
       */

      nhole = 0;

      for (i=0; i < L; i++)
    	{
    	  for (j=0; j < L; j++)
    	    {
    	      r=uni();

    	      if(r < rho)
        		{
        		  map[i][j] = 0;
        		}
    	      else
        		{
        		  nhole++;
        		  map[i][j] = nhole;
        		}
    	    }
    	}
      printf("n hole is %d\n",nhole);
      printf("percolate: rho = %f, actual density = %f\n",
	      rho, 1.0 - ((double) nhole)/((double) L*L) );
    }

  /*
   * Use broadcast and copy-back to distribute the map. This is not as
   * elegant as using scatter in the 1D decomposition, but generalises
   * to a 2D decomposition (while scatter does not). Use &map[0][0]
   * syntax as this also work for dynamically allocated arrays.
   */

  MPI_Bcast(&map[0][0], L*L, MPI_INT, 0, comm);

  /*
   * Copy the appropriate section back to smallmap. Could probably
   * eliminate use of smallmap in its entirety, but leave it in here
   * for simplicity.
   */

  for (i=0; i < M; i++)
    {
      for (j=0; j < N; j++)
        {
          smallmap[i][j] = map[rank*M+i][j];
        }
    }

  for (i=1; i <= M; i++)
    {
      for (j=1; j <= N; j++)
    	{
    	  old[i][j] = smallmap[i-1][j-1];
    	}
    }

  /*
   * Initialise the old array: copy the LxL array smallmap to the centre of
   * old, and set the halo values to zero.
   */

   for (i=1; i <= M; i++)
    {
      for (j=1; j <= N; j++)
	{
	  old[i][j] = smallmap[i-1][j-1];
	}
    }

   for (i=0; i <= M+1; i++)  // zero the bottom and top halos
    {
      old[i][0]   = 0;
      old[i][N+1] = 0;
    }

   for (j=0; j <= N+1; j++)  // zero the left and right halos
    {
      old[0][j]   = 0;
      old[M+1][j] = 0;
    }

  step = 1;
  nchange = 1;

  while (step <= maxstep)
    {
      /*
       *  Swap halos up and down
       */

      /*
       * Communications is done using the sendrecv routine; a proper
       * solution would use non-blocking communications (e.g. some
       * combination of issend/recv or ssend/irecv)
       */

      MPI_Sendrecv(&old[M][1], N, MPI_INT, next, tag,
		   &old[0][1], N, MPI_INT, prev, tag,
		   comm, &status);

      MPI_Sendrecv(&old[1][1], N, MPI_INT, prev, tag,
		   &old[M+1][1], N, MPI_INT, next, tag,
		   comm, &status);

      nchangelocal = 0;

      for (i=1; i<=M; i++)
	{
	  for (j=1; j<=N; j++)
	    {
	      oldval = old[i][j];
	      newval = oldval;

	      /*
	       * Set new[i][j] to be the maximum value of old[i][j]
	       * and its four nearest neighbours
	       */

	      if (oldval != 0)
    		{
    		  if (old[i][j-1] > newval) newval = old[i][j-1];
    		  if (old[i][j+1] > newval) newval = old[i][j+1];
    		  if (old[i-1][j] > newval) newval = old[i-1][j];
    		  if (old[i+1][j] > newval) newval = old[i+1][j];

    		  if (newval != oldval)
    		    {
    		      ++nchangelocal;
    		    }
    		}

    	      new[i][j] = newval;
    	    }
    	}

      /*
       *  Compute global number of changes on rank 0
       */

      MPI_Reduce(&nchangelocal, &nchange, 1, MPI_INT, MPI_SUM, 0, comm);

      /*
       *  Report progress every now and then
       */

      if (step % printfreq == 0)
    	{
    	  if (rank == 0)
    	    {
                  printf("percolate: changes on step %d is %d\n",
                         step, nchange);
    	    }
    	}

      /*
       *  Copy back in preparation for next step, omitting halos
       */
      int total;
      total = 0;
      for (i=1; i<=M; i++)
    	{
    	  for (j=1; j<=N; j++)
    	    {
            total += new[i][j];
    	      old[i][j] = new[i][j];

    	    }
    	}
      if (step % printfreq == 0)
      {
        if (rank == 0)
          {
                  printf("Average value of map array on step %d is %d\n",
                         step, total/M/N);
          }
      }

      step++;
      // if (rank == 0 && nchange == 0 && step >= maxstep/2){
      //   printf("???");
      //   step = maxstep + 1;
      // }
    }
  if(rank == 0){
    t = clock() - t;
    double time_taken = ((double)t)/CLOCKS_PER_SEC; // in seconds
    printf("fun() took %f seconds to execute. Average time for step is %f \n", time_taken, time_taken/maxstep);
  }

  /*
   *  We set a maximum number of steps to ensure the algorithm always
   *  terminates. However, if we hit this limit before the algorithm
   *  has finished then there must have been a problem (e.g. the value
   *  of maxstep is too small)
   */
  if (rank == 0)
    {
      if (nchange != 0)
    	{
        printf("percolate: WARNING max steps = %d reached but nchange != 0\n",
    	     maxstep);
    	}
    }

  /*
   *  Copy the centre of old, excluding the halos, into smallmap
   */

  for (i=1; i<=M; i++)
    {
      for (j=1; j<=N; j++)
    	{
    	  smallmap[i-1][j-1] = old[i][j];
    	}
    }

  /*
   *  Use copy and reduce to collect the map.  This is not as elegant
   *  as using gather in the 1D decomposition, but generalises to a 2D
   *  decomposition (while gather does not).  Use &map[0][0] syntax in
   *  reduce as this also works for dynamically allocated arrays.
   */

  // Zero maptmp

  for (i=0; i < L; i++)
    {
      for (j=0; j < L; j++)
        {
          maptmp[i][j] = 0;
        }
    }

  /*
   *  Copy smallmap to correct place in maptmp. Could probably
   *  eliminate smallmap entirely, but leave here for simplicity.
   */

  for (i=0; i < M; i++)
    {
      for (j=0; j < N; j++)
        {
          maptmp[rank*M+i][j] = smallmap[i][j];
        }
    }

  MPI_Reduce(&maptmp[0][0], &map[0][0], L*L, MPI_INT, MPI_SUM, 0, comm);

  /*
   *  Test to see if percolation occurred by looking for positive numbers
   *  that appear on both the top and bottom edges
   */

  if (rank == 0)
    {
      perc = 0;

      for (itop=0; itop < L; itop++)
        {
          if (map[itop][L-1] > 0)
            {
              for (ibot=0; ibot < L; ibot++)
                {
                  if (map[ibot][0] == map[itop][L-1])
                    {
                      perc = 1;
                    }
                }
            }
        }

      if (perc != 0)
	{
	  printf("percolate: cluster DOES percolate\n");
	}
      else
	{
	  printf("percolate: cluster DOES NOT percolate\n");
	}

      /*
       *  Write the map to the file "map.pgm", displaying the two
       *  largest clusters. If the last argument here was 3, it would
       *  display the three largest clusters etc. The picture looks
       *  cleanest with only a single cluster, but multiple clusters
       *  are useful for debugging.
       */


      mapwritedynamic("map.pgm", map, L, 2);

      FILE *f = fopen("map_parallel.data", "wb");
      for(int i=0;i<L;++i){
        for(int j=0;j<L;++j){
         fwrite((char*)&(map[i][j]),sizeof(char),sizeof(int), f);
       }
      }
      // fwrite(map, sizeof(char), sizeof(map), f);
      fclose(f);

      /* File pointer to hold reference of input file */
      FILE * fPtr1;
      FILE * fPtr2;
      char path1[100];
      char path2[100];

      int diff;
      int line, col;


      /* Input path of files to compare */
      // printf("Enter path of first file: ");
      // scanf("%s", path1);
      // printf("Enter path of second file: ");
      // scanf("%s", path2);
      strcpy(path1, "map.data");
      strcpy(path2, "map_parallel.data");


      /*  Open all files to compare */
      fPtr1 = fopen(path1, "r");
      fPtr2 = fopen(path2, "r");

      /* fopen() return NULL if unable to open file in given mode. */
      if (fPtr1 == NULL || fPtr2 == NULL)
      {
          /* Unable to open file hence exit */
          perror("Error");
          printf("\nUnable to open file.\n");
          printf("Please check whether file exists and you have read privilege.\n");
          exit(EXIT_FAILURE);
      }


      /* Call function to compare file */
      diff = compareFile(fPtr1, fPtr2, &line, &col);

      if (diff == 0)
      {
          printf("\nBoth files are equal.\n");
      }
      else
      {
          printf("\nFiles are not equal.\n");
          printf("Line: %d, col: %d\n", line, col);
      }


      /* Finally close files to release resources */
      fclose(fPtr1);
      fclose(fPtr2);


        }

  free(old);
  free(new);
  free(map);
  free(maptmp);
  free(smallmap);

  MPI_Finalize();

  return 0;
}
