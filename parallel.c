/*
 * Simple parallel program to test for percolation of a cluster
 */

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "percolate_her.h"
#include "arralloc.h"

// void print_matrix1(int m, int n, int map[m][n], int rank)
// {
// 	if (rank == 2) {
// 		for (int i = 0; i < m; ++i) {
// 			for (int j = 0; j < n; ++j) {
// 				printf("%d ", map[i][j]);
// 			}
// 			printf("\n");
// 		}
// 	}
// }

// void print_matrix(int m, int n, int map[m][n], int rank)
// {
// }

/*
 * Simple parallel program to test for percolation of a cluster.
 */

int main(int argc, char *argv[])
{

	/*
	 * Variables that define the simulation
	 */

	int seed;
	double rho;

	/*
	 * Local variables
	 */

	int i, j, nhole, step, maxstep;
	int oldval, newval;
	int nchangelocal;
	int nchange;
	int printfreq;
	int itop, ibot, perc;
	double r;

	/*
	 * MPI variables
	 */


	int size, rank;
	//int tag = 1;

	//初始化top
	int dims[2] = {0, 0};
	int periods[2] = {0, 0};

	int reorder = 0;
	int coords[2];
	int ndims = 2;
	int disp = 1;
	int left_nbr, right_nbr, up_nbr, down_nbr;


	//初始化MPI
		MPI_Comm comm = MPI_COMM_WORLD;
	//新的communicator
		MPI_Comm new_comm = MPI_COMM_WORLD;

	MPI_Status status_a, status_b, status_c, status_d;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(comm, &size);
	MPI_Comm_rank(comm, &rank);
	MPI_Request requests_a, requests_b, requests_c, requests_d;
	//MPI_Comm_rank(newcommunicator, &rank);


	//初始化MPi_top
		MPI_Dims_create(size, ndims, dims);
	//在由指定维数构成的笛卡尔网络上分解给定数量的进程

	if (rank == 0) {
		printf("dims is %d, %d\n", dims[0], dims[1]);
	}

	int L;
	L = atoi(argv[2]);
	int NPROC;
	NPROC = atoi(argv[3]);
	//(要分解的进程数，维数，分配给每个维度进程的数组)
	int nheight = dims[0];
	int nwidth = dims[1];
	int M = L / nheight;
	int N = L / nwidth;
	MPI_Cart_create(comm, ndims, dims, periods, reorder, &new_comm);

	int my_rank;
	MPI_Comm_rank(new_comm, &my_rank);
	MPI_Cart_coords(new_comm, my_rank, 2, coords);
	if (rank == 0)
		printf("MPI process %d, located at (%d %d).\n", my_rank, coords[0], coords[1]);

	MPI_Cart_shift(new_comm, 0, disp, &up_nbr, &down_nbr);
	MPI_Cart_shift(new_comm, 1, disp, &left_nbr, &right_nbr);

	// (communicator, 位移坐标维数，虚拟移动拓扑的单位数，目的地，不是周期�119 � 的可以用MPI_PROC_NULL)
		if (rank == 0)
		printf("rank %d coords : %d,%d nbrs: up %d, down %d, left %d, right %d\n", my_rank, coords[0], coords[1], up_nbr, down_nbr, left_nbr, right_nbr);


	//初始化左右交换的数据结构
		MPI_Datatype vector_LR;
	MPI_Type_vector(M, 1, N + 2, MPI_INT, &vector_LR);
	MPI_Type_commit(&vector_LR);

  int height, width;
	height = dims[0];
	width = dims[1];
	int count;
	int blocklength;
	count = (int)ceil((double)L / height);
	blocklength = (int)ceil((double)L / width);
	M = count;
	N = blocklength;
	MPI_Datatype vector, halo_j, halo_i;
	MPI_Type_vector(count, blocklength, L + size, MPI_INT, &vector);
	MPI_Type_commit(&vector);

	MPI_Type_vector(1, blocklength, N + 2, MPI_INT, &halo_j);
	MPI_Type_commit(&halo_j);
	MPI_Type_vector(count, 1, N + 2, MPI_INT, &halo_i);
	MPI_Type_commit(&halo_i);

	/*
	 * Define the main arrays for the simulation
	 */

	// int old[M + 2][N + 2];
	// int new[M + 2][N + 2];
	//
	// /*
	//  * Additional array WITHOUT halos for initialisation and IO. This is
	//  * of size LxL because, even in our parallel program, we do these two
	//  * steps in serial
	//  */
	//
	// int map[L][L];
	// int maptmp[L][L];
	//
	// /*
	//  * Array to store local part of map
	//  */
	//
	// int smallmap[M][N];
	MPI_Barrier(comm);
	double time_start = MPI_Wtime();

	int **old;
	int **new;
	int **map;
	int **maptmp;
	int **smallmap;
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
	 * Non-periodic boundary conditions
	 *
	 * Note that the special rank of MPI_PROC_NULL is a "black hole" for
	 * communications. Using this value for processes off the edges of
	 * the image means there is no additional logic needed to ensure
	 * processes at the edges do not attempt to send to or receive from
	 * invalid ranks (i.e. rank = -1 and rank = NPROC).
	 *
	 * Proper solution would compute neighbours with a Cartesian topology
	 * and MPI_Cart_shift, where MPI_PROC_NULL is assigned automatically.
	 */

	/*
	 * if (next >= size) { next = MPI_PROC_NULL; }
	 *
	 * if (prev < 0) { prev = MPI_PROC_NULL; }
	 */

	if (NPROC != size) {
		if (rank == 0) {
			printf("percolate: ERROR, NPROC = %d but running on %d\n",
			       NPROC, size);
		}
		MPI_Finalize();
		return 0;
	}
	if (argc != 4)
    {
      if (rank == 0)
	{
	  printf("Usage: percolate <seed> <size> <NPROC>\n");
	}
}
	/*
	 * Update for a fixed number of steps, periodically report progress
	 */

	maxstep = 5 * L;
	printfreq = 100;

	if (rank == 0) {
		printf("percolate: running on %d process(es)\n", size);

		/*
		 * Set most important value: the rock density rho (between 0
		 * and 1)
		 */

		rho = 0.4064;

		/*
		 * Set the randum number seed and initialise the generator
		 */
		seed = atoi(argv[1]);

		printf("percolate: L = %d, rho = %f, seed = %d, maxstep = %d\n",
		       L, rho, seed, maxstep);

		rinit(seed);

		/*
		 * Initialise map with density rho. Zero indicates rock, a
		 * positive value indicates a hole. For the algorithm to
		 * work, all the holes must be initialised with a unique
		 * integer
		 */

		nhole = 0;

		for (i = 0; i < L; i++) {
			for (j = 0; j < L; j++) {
				r = uni();

				if (r < rho) {
					map[i][j] = 0;
				} else {
					nhole++;
					map[i][j] = nhole;
				}

				//map[i][j] = i * 1000 + j;
			}
		}

		printf("percolate: rho = %f, actual density = %f\n",
		       rho, 1.0 - ((double)nhole) / ((double)L * L));
	}
	/*
	 * Use broadcast and copy-back to distribute the map. This is not as
	 * elegant as using scatter in the 1D decomposition, but generalises
	 * to a 2D decomposition (while scatter does not). Use &map[0][0]
	 * syntax as this also work for dynamically allocated arrays.
	 */

	MPI_Bcast(&map[0][0], L * L, MPI_INT, 0, comm);

	/*
	 * Copy the appropriate section back to smallmap. Could probably
	 * eliminate use of smallmap in its entirety, but leave it in here
	 * for simplicity.
	 */
          // if (rank == 0) {
          // 	printf("map is \n");
          // 	for (i = 0; i < L; i++) {
          // 		for (j = 0; j < L; j++) {
          // 			//smallmap[i][j] = map[coords[0] * M + i][coords[1] * N + j];
          // 			printf("%d ", map[i][j]);
          // 		}
          // 		printf(" \n ");
          // 	}
          // }

	//分发到小地图
	MPI_Cart_coords(new_comm, my_rank, 2, coords);

	//printf("coords=%d  in rank :%d\n", coords[0], my_rank);
	//printf("coords_y=%d  in rank :%d\n", coords[1], my_rank);
 //	printf("small map is \n");
	for (i = 0; i < M; i++) {
		for (j = 0; j < N; j++) {
			smallmap[i][j] = map[coords[0] * M + i][coords[1] * N + j];
			//if (rank == 2)
				//printf("%d ", smallmap[i][j]);
		}
		//printf(" \n ");
	}

	// if (rank == 0)
	// {
	// 	for (int k = 0; k < size; k++)
	// 	{
	// 		MPI_Cart_coords(new_comm, k, 2, coords);
	//
	// 		int x = coords[0] * count;
	// 		int y = coords[1] * blocklength;
	//
	// 		MPI_Bsend(&(map[x][y]), 1, vector, k, 0, comm);
	// 	}
	// }

	//print_matrix(M, N, smallmap, rank);



	//printf("M=%d ------ :\n", M);
	//printf("N=%d ------ :\n", N);

	//printf("coords=%d  in rank :%d ", coords[0], my_rank);


	// for (i = 1; i <= M; i++) {
	// 	for (j = 1; j <= N; j++) {
	// 		//old[i][j] = smallmap[i - 1][j - 1];
	// 	}
	// }

	/*
	 * Initialise the old array: copy the LxL array smallmap to the
	 * centre of old, and set the halo values to zero.
	 */
	//old � ��上halo
		for (i = 1; i <= M; i++) {
		for (j = 1; j <= N; j++) {
			old[i][j] = smallmap[i - 1][j - 1];
		}
	}

	for (i = 0; i <= M + 1; i++)
		//zero the bottom and top halos
	{
		old[i][0] = 0;
		old[i][N + 1] = 0;
	}

	for (j = 0; j <= N + 1; j++)
		//zero the left and right halos
	{
		old[0][j] = 0;
		old[M + 1][j] = 0;
	}

	/*
	 * 打印old 没有问题 printf("Old: \n"); for (i=0; i < M+2; i++)
	 * { for (j=0; j < N+2; j++) { //smallmap[i][j] =
	 * map[coords[0]*M+i][coords[1]*N+j];
	 *
	 * printf("%d ",old[i][j]); } printf(" \n "); }
	 */
	step = 1;
	nchange = 1;

	//if (rank == 0)
	//	printf("old map is \n");
	//print_matrix(M + 2, N + 2, old, rank);

	while (step <= maxstep) {
		/*
		 * Swap halos up and down
		 */

		/*
		 * Communications is done using the sendrecv routine; a
		 * proper solution would use non-blocking communications
		 * (e.g. some combination of issend/recv or ssend/irecv)
		 */


		MPI_Issend(&old[M][1], N, MPI_INT, down_nbr, 1, comm, &requests_a);
		MPI_Issend(&old[1][1], N, MPI_INT, up_nbr, 2, comm, &requests_b);
		MPI_Issend(&old[1][N], 1, vector_LR, right_nbr, 1, comm, &requests_c);
		MPI_Issend(&old[1][1], 1, vector_LR, left_nbr, 2, comm, &requests_d);


		MPI_Recv(&old[0][1], N, MPI_INT, up_nbr, 1, comm, &status_a);
		MPI_Recv(&old[M + 1][1], N, MPI_INT, down_nbr, 2, comm, &status_b);
		MPI_Recv(&old[1][0], 1, vector_LR, left_nbr, 1, comm, &status_c);
		MPI_Recv(&old[1][N + 1], 1, vector_LR, right_nbr, 2, comm, &status_d);

		MPI_Wait(&requests_a, &status_a);
		MPI_Wait(&requests_b, &status_b);
		MPI_Wait(&requests_c, &status_c);
		MPI_Wait(&requests_d, &status_d);

		/*printf每次的更改后的地图
		 * if (rank == 2) { printf("neighbour left %d, right %d, up
		 * %d, down %d\n", left_nbr, right_nbr, up_nbr, down_nbr);
		 * printf("old mao in : --------%d\n ",my_rank);
		 *
		 * for (i=0; i < M+2; i++) { for (j=0; j < N+2; j++) {
		 * //smallmap[i][j] = map[coords[0]*M+i][coords[1]*N+j];
		 *
		 * printf("%d ",old[i][j]); } printf(" \n "); } }
		 */

		nchangelocal = 0;

		for (i = 1; i <= M; i++) {
			for (j = 1; j <= N; j++) {
				oldval = old[i][j];
				newval = oldval;


				//Set new[i][j] to be the maximum value of old[i][j]
				//  and its four nearest neighbours
				if  (oldval != 0) {
					if (old[i][j - 1] > newval)
						newval = old[i][j - 1];
					if (old[i][j + 1] > newval)
						newval = old[i][j + 1];
					if (old[i - 1][j] > newval)
						newval = old[i - 1][j];
					if (old[i + 1][j] > newval)
						newval = old[i + 1][j];

					if (newval != oldval) {
						++nchangelocal;
					}
				}
				    new[i][j] = newval;
			}
		}


		/*
		 * Compute global number of changes on rank 0
		 */

		MPI_Reduce(&nchangelocal, &nchange, 1, MPI_INT, MPI_SUM, 0, comm);

		/*
		 * Report progress every now and then
		 */

		if (step % printfreq == 0) {
			if (rank == 0) {
				printf("percolate: changes on step %d is %d\n",
				       step, nchange);
			}
		}
		/*
		 * Copy back in preparation for next step, omitting halos
		 */

		for (i = 1; i <= M; i++) {
			for (j = 1; j <= N; j++) {
				old[i][j] = new[i][j];
			}
		}

		step++;

	}
		MPI_Barrier(comm);
	if(rank == 0){
		double time_end = MPI_Wtime();
		printf("fun() took %f seconds to execute. Average time for step is %f seconds\n", time_end - time_start, (time_end-time_start)/step);
	}


	/*
	 * We set a maximum number of steps to ensure the algorithm always
	 * terminates. However, if we hit this limit before the algorithm has
	 * finished then there must have been a problem (e.g. the value of
	 * maxstep is too small)
	 */

	if (rank == 0) {
		if (nchange != 0) {
			printf("percolate: WARNING max steps = %d reached but nchange != 0\n",
			       maxstep);
		}
	}
	/*
	 * Copy the centre of old, excluding the halos, into smallmap
	 */

	for (i = 1; i <= M; i++) {
		for (j = 1; j <= N; j++) {
			smallmap[i - 1][j - 1] = old[i][j];
		}
	}

	/*
	 * Use copy and reduce to collect the map.  This is not as elegant as
	 * using gather in the 1D decomposition, but generalises to a 2D
	 * decomposition (while gather does not).  Use &map[0][0] syntax in
	 * reduce as this also works for dynamically allocated arrays.
	 */

	//Zero maptmp

		for (i = 0; i < L; i++) {
		for (j = 0; j < L; j++) {
			maptmp[i][j] = 0;
		}
	}

	/*
	 * Copy smallmap to correct place in maptmp. Could probably eliminate
	 * smallmap entirely, but leave here for simplicity.
	 */

	for (i = 0; i < M; i++) {
		for (j = 0; j < N; j++) {

			maptmp[coords[0] * M + i][coords[1] * N + j] = smallmap[i][j];

		}
	}

	MPI_Reduce(&maptmp[0][0], &map[0][0], L * L, MPI_INT, MPI_SUM, 0, comm);

//打印最终的地图

	// if (rank == 0) {
	// 	for (i = 0; i < L; i++) {
	// 		for (j = 0; j < L; j++) {
	// 			//smallmap[i][j] = map[coords[0] * M + i][coords[0] * N + j];
	// 			printf("%d ", map[i][j]);
	// 		}
	// 		printf(" \n ");
	// 	}
	// }
	/*
	 * Test to see if percolation occurred by looking for positive
	 * numbers that appear on both the top and bottom edges
	 */

	if (rank == 0) {
		perc = 0;

		for (itop = 0; itop < L; itop++) {
			if (map[itop][L - 1] > 0) {
				for (ibot = 0; ibot < L; ibot++) {
					if (map[ibot][0] == map[itop][L - 1]) {
						perc = 1;
					}
				}
			}
		}

		if (perc != 0) {
			printf("percolate: cluster DOES percolate\n");
		} else {
			printf("percolate: cluster DOES NOT percolate\n");
		}

		/*
		 * Write the map to the file "map.pgm", displaying the two
		 * largest clusters. If the last argument here was 3, it
		 * would display the three largest clusters etc. The picture
		 * looks cleanest with only a single cluster, but multiple
		 * clusters are useful for debugging.
		 */

		mapwritedynamic("map.pgm", map, L, 2);
	}
	MPI_Finalize();

	  free(old);
	  free(new);
	  free(map);
	  free(maptmp);
	  free(smallmap);

	return 0;
}
