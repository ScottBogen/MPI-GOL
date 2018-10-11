#include <stdio.h>
#include <mpi.h>
#include <assert.h>
#include <sys/time.h>
#include <stdlib.h>


#define DEAD  0
#define ALIVE 1

int elements_per_row;
int elements_per_subgrid;

typedef struct Cell {
    int status;         // 0 == DEAD, 1 == ALIVE
    int next_status;
}Cell;


int FindNorth(Cell* subgrid, int index, int* top);
int FindNorthEast(Cell* subgrid, int index, int* top);
int FindEast(Cell* subgrid, int index);
int FindSouthEast(Cell* subgrid, int index, int* bottom);
int FindSouth(Cell* subgrid, int index, int* bottom);
int FindSouthWest(Cell* subgrid, int index, int* bottom);
int FindWest(Cell* subgrid, int index);
int FindNorthWest(Cell* subgrid, int index, int* top);

int main(int argc,char *argv[]) {
    int rank,p;
    int n = 16;
    int g = 100;
    int i;

    srand(time(NULL));
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);        // current rank
    MPI_Comm_size(MPI_COMM_WORLD, &p);           // # total processes

    printf("my rank=%d\n",rank);
    printf("Rank=%d: number of processes =%d\n",rank,p);
    printf("grid size n*n = %d\n", n*n);
    printf("subgrid size n*n/p = %d\n", n*n/p);
    assert(p>=2);

    if (!(n%p==0)) {
        printf("Error: n%p!=0");
        return 0;
    }

    Cell* subgrid = (Cell*) malloc(sizeof(subgrid) * n*n/p);    // this is each rank's portion of the grid
    elements_per_row = n;
    elements_per_subgrid = n*n/p;

    printf("elements per row: %d\nelements per subgrid:%d\n\n\n", elements_per_row, elements_per_subgrid);

    GenerateInitialGoL(rank, p, n, subgrid);
    //DisplayGoL(rank, p, n, subgrid);


    // for each generation:
    for (i = 0; i < g; i++) {
    	MPI_Barrier(MPI_COMM_WORLD);
	Simulate(rank, p, n, subgrid);
	Update(subgrid);
	MPI_Barrier(MPI_COMM_WORLD);

	if  (i % 10 == 0) {
	    if (rank == 0) {
	        printf("Grid at generation#%d:\n", i);
	    }
            DisplayGoL(rank, p, n, subgrid);
	    printf("\n");
	}
    }

    MPI_Finalize();
    return 0;
}

void Update(Cell* subgrid) {
	int i;
	for (i = 0; i < elements_per_subgrid; i++) {
		subgrid[i].status = subgrid[i].next_status;
	}
}

void GenerateInitialGoL(int rank, int p, int n, Cell* subgrid) {
    if (rank==0) {
	int i;
        for (i = 1; i < p; i++) {
            // generate random number
            int random = (rand() % 93563) + 1;
            // MPI SEND random to rank#i
            MPI_Send(&random, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
        }

        // just for rank0
        int random = (rand() % 93563) + 1;
        for (i = 0; i < elements_per_subgrid; i++) {
            int tmp = rand() % random;

            if (tmp%2) {
                subgrid[i].status = DEAD;
                subgrid[i].next_status = DEAD;
            }

            else {
                subgrid[i].status = ALIVE;
                subgrid[i].next_status = ALIVE;
	    }

        }
    }

    else {
        // others wait for MPI receive
        int random, i;
	MPI_Status status;
        //	printf("Rank #%d awaiting msg\n", rank);
   	MPI_Recv(&random,1,MPI_INT,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&status);

        for (i = 0; i < elements_per_subgrid; i++) {
            int tmp = (rand() % random) + 1;
            if (tmp%2==0) {
                subgrid[i].status = DEAD;
                subgrid[i].next_status = DEAD;
            }
            else {
                subgrid[i].status = ALIVE;
                subgrid[i].next_status = ALIVE;
            }
        }
    }
}

void Simulate(int rank, int p, int n, Cell* subgrid) {
    int i, j;

    MPI_Status status;
    int* my_top    = (int*) calloc(elements_per_row, sizeof(int));
    int* my_bottom = (int*) calloc(elements_per_row, sizeof(int));
    int* top       = (int*) calloc(elements_per_row, sizeof(int));
    int* bottom    = (int*) calloc(elements_per_row, sizeof(int));

    for (i = 0; i < elements_per_row; i++) {
      my_top[i] = subgrid[i].status;							// first row
      my_bottom[i] = subgrid[elements_per_subgrid - elements_per_row + i].status;	// last row
    }

    MPI_Send(my_top,    elements_per_row, MPI_INT, (rank-1+p)%p, 0, MPI_COMM_WORLD);	// previous rank
    MPI_Send(my_bottom, elements_per_row, MPI_INT, (rank+1)%p,   0, MPI_COMM_WORLD);	// next rank
    //printf("[%d] MPI Sent top and bottom\n", rank);

    // receive my top and bottom neighbor's subgrids
    MPI_Recv(bottom, elements_per_row, MPI_INT,(rank+1)%p,MPI_ANY_TAG,MPI_COMM_WORLD,&status);  // next rank, first row
    MPI_Recv(top, elements_per_row, MPI_INT,(rank-1+p)%p,MPI_ANY_TAG,MPI_COMM_WORLD,&status);	// previous rank, last row

    //printf("[%d] MPI Received top/bottom\n", rank);

    //for (j = 0; j < elements_per_row; j++) {
    //    printf("[%d] Testing top[%d]: %d\n",rank, j, top[j]);
    //}


    //printf("[%d] Barrier up #2\n", rank);
    // block until everybody is ready
    MPI_Barrier(MPI_COMM_WORLD);
    //printf("[%d] Barrier down #2\n", rank);


    for (i = 0; i < n*n/p; i++) {
      int sum = AnalyzeNeighbors(rank, subgrid, &subgrid[i], i, top, bottom);
      if (sum < 3 || sum > 5) {
        subgrid[i].next_status = DEAD;
      }
      else {
        subgrid[i].next_status = ALIVE;
      }
    }

    // free 4 arrs
    free(my_top);
    free(my_bottom);
    free(top);
    free(bottom);

}

int AnalyzeNeighbors(int rank, Cell* subgrid, Cell* cell, int index, int* top, int* bottom) {
    int k = 0;

    k += FindNorth(subgrid, index, top);
    //if (k > 0) { printf("North --- "); sum+=k; k=0;}
    k += FindNorthEast(subgrid, index, top);
    //if (k > 0) { printf("NorthEast --- "); sum+=k; k=0; }
    k += FindEast(subgrid, index);
   // if (k > 0) { printf("East --- "); sum+=k; k=0; }
    k += FindSouthEast(subgrid, index, bottom);
    //if (k > 0) { printf("SouthEast --- "); sum+=k; k=0; }
    k += FindSouth(subgrid, index, bottom);
    //if (k > 0) { printf("South --- "); sum+=k; k=0; }
    k += FindSouthWest(subgrid, index, bottom);
    //if (k > 0) { printf("SouthWest --- "); sum+=k; k=0; }
    k += FindWest(subgrid, index);
    //if (k > 0) { printf("West --- "); sum+=k; k=0; }
    k += FindNorthWest(subgrid, index, top);
    //if (k > 0) { printf("NorthWest --- "); sum+=k; k=0; }

    //printf("---[%d::%d]---Sum = %d---\n", rank,index, sum);

    return k;
}


// FIRST CALL:
// subgrid = [ cell0, cell1, ...]
// i = 0
// top = {1,0,1,0,0,0,0,0}
// return = 1
int FindNorth(Cell* subgrid, int i, int* top) {
    int temp;

    if (i < elements_per_row) {     // at north border
        temp = top[i];
    }
    else {
        temp = subgrid[i - elements_per_row].status;
    }
    return temp;
}


int FindNorthEast(Cell* subgrid, int i, int* top) {
    int temp;
    int temp_index;


    if (i < elements_per_row) {     // at north border
       // printf("A\n");
        if ((i+1) % elements_per_row == 0) { // at east border
            temp = top[i -elements_per_row + 1];
       //     printf("B\n");
        }
        else {
            temp = top[i+1];
       //     printf("C = %d\n", temp);
        }
    }

    else {
        temp_index = i - elements_per_row;  // normal move north

        if ((temp_index+1) % elements_per_row == 0) { // at east border
            temp = subgrid[temp_index-elements_per_row+1].status;
        }
        else {  // normal case
            temp = subgrid[temp_index+1].status;
        }

    }

    return temp;
}

int FindEast(Cell* subgrid, int i) {
    int temp;
    if ((i + 1) % elements_per_row == 0) {  // at east border
        temp = subgrid[i - elements_per_row + 1].status;
    }
    else {      // normal case
        temp = subgrid[i+1].status;
    }
    return temp;
}

int FindSouthEast(Cell* subgrid, int i, int* bottom) {
    int temp;
    int temp_index;

    // SOUTH
    if (i >= elements_per_subgrid - elements_per_row) {     // at south border
        temp_index = i % elements_per_row;

        if ((temp_index+1) % elements_per_row == 0) { // at east border
            temp = bottom[temp_index - elements_per_row + 1];
        }
        else {
            temp = bottom[temp_index + 1];
        }
    }

    else {
        temp_index = i + elements_per_row;  // normal move south

        if ((temp_index+1) % elements_per_row == 0) { // at east border
            temp = subgrid[temp_index - elements_per_row + 1].status;  // wrap around
        }
        else {  // normal case
            temp = subgrid[temp_index + 1].status;
        }
    }
    return temp;
}

int FindSouth(Cell* subgrid, int i, int* bottom) {
    int temp;
    if (i >= elements_per_subgrid - elements_per_row) {     // at north border
        temp = bottom[i % elements_per_row];
    }
    else {
        temp = subgrid[i + elements_per_row].status;
    }
    return temp;
}

int FindSouthWest(Cell* subgrid, int i, int* bottom) {
  int temp;
  int temp_index;

  // SOUTH
  if (i >= elements_per_subgrid - elements_per_row) {     // at south border
      temp_index = i % elements_per_row;

      if ((temp_index) % elements_per_row == 0) { 	  // at west border
          temp = bottom[temp_index + elements_per_row - 1];
      }
      else {
          temp = bottom[temp_index+1];
      }
  }

  else {
      temp_index = i + elements_per_row;  // normal move south

      // WEST
      if (temp_index % elements_per_row == 0) {    // at west border
          temp = subgrid[temp_index + elements_per_row - 1].status;
      }
      else {
          temp = subgrid[temp_index - 1].status;
      }

  }

  return temp;
}

int FindWest(Cell* subgrid, int i) {
    int temp;

    if (i % elements_per_row == 0) {    // at west border
        temp = subgrid[i + elements_per_row - 1].status;
    }
    else {
        temp = subgrid[i - 1].status;
    }
    return temp;
}

int FindNorthWest(Cell* subgrid, int i, int* top) {
    int temp;
    int temp_index;

    if (i < elements_per_row) {     // at north border
        if (i % elements_per_row == 0) {	 // at west border
            temp = top[i + elements_per_row - 1];
        }
        else {
            temp = top[i - 1];
        }
    }

    else {
        temp_index = i - elements_per_row;  // normal move north

        if (temp_index % elements_per_row == 0) {    // at west border
            temp = subgrid[temp_index + elements_per_row - 1].status;
        }
        else {
            temp = subgrid[temp_index - 1].status;
        }
    }
    return temp;
}

void DisplayGoL(int rank, int p, int n, Cell* subgrid) {
	int dummy;
	MPI_Status status;
	int j, i;

	int fullgrid[p][elements_per_subgrid];	// array to hold all subgrids
	int tempgrid[elements_per_subgrid];

	if (rank==0) {

	 	for (i = 0; i < elements_per_subgrid; i++) {
                    tempgrid[i] = subgrid[i].status;
                }

                // add that tempgrid to the fullgrid
                for (j = 0; j < elements_per_subgrid; j++) {
		    fullgrid[0][j] = tempgrid[j]; 
		}

		//fullgrid[0] = tempgrid;
		//memcpy(fullgrid[0], tempgrid, sizeof(int) * elements_per_subgrid);

                // clear tempgrid;
                memset(tempgrid, 0, sizeof(int) * elements_per_subgrid);

                // get remaining processes data into subgrid
                for (i = 1; i < p; i++) {
                    MPI_Recv(tempgrid, elements_per_subgrid, MPI_INT,i,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
                    //fullgrid[i] = tempgrid;

	            for (j = 0; j < elements_per_subgrid; j++) {
                        fullgrid[i][j] = tempgrid[j];
                    }

                    //memcpy(fullgrid[i], tempgrid, sizeof(int) * elements_per_subgrid);
		    memset(tempgrid, 0, sizeof(int) * elements_per_subgrid);
                }

                // fullgrid is now filled up, let's print it out now

		for (j = 0; j < p; j++) {
			for (i = 0; i < elements_per_subgrid; i++) {
                    		if (i%elements_per_row==0) {
                    			printf("\n");
                    		}

            	    		if (fullgrid[j][i] == ALIVE) { printf("X "); }
            	    		else if (fullgrid[j][i] == DEAD) { printf("_ "); }
           	    		else { printf("? "); }	// this shouldn't happen
			}
		}

        }

        // rank x (from 1 to p-1) will package its contents into a handy array and send the array to rank0 
        else {
                int tempgrid[elements_per_subgrid];
                for (i = 0; i < elements_per_subgrid; i++) {
                    tempgrid[i] = subgrid[i].status;
		}
		//						 root
		MPI_Send(tempgrid, elements_per_subgrid, MPI_INT, 0, 0, MPI_COMM_WORLD);
	}

    return;
}
