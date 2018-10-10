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

    struct Cell* N;
    struct Cell* NE;
    struct Cell* E;
    struct Cell* SE;
    struct Cell* S;
    struct Cell* SW;
    struct Cell* W;
    struct Cell* NW;
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

    int n = 8;

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
    DisplayGoL(rank, p, n, subgrid);

    MPI_Barrier(MPI_COMM_WORLD);
    Simulate(rank, p, n, subgrid);
    //DisplayGoL(rank, p, n, subgrid);

    MPI_Finalize();
    return 0;
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
        for (i = 0; i < n*n/p; i++) {
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
            else if (tmp%2==1) {
                subgrid[i].status = ALIVE;
                subgrid[i].next_status = ALIVE;
            }
	    else {
		printf("WHAT WHAT WHAT WHAT WHAT :::: tmp%2==%d ON INT %d~~~~~\n", tmp%2);
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

    //MPI_Barrier(MPI_COMM_WORLD);

    for (i = 0; i < elements_per_row; i++) {
      if (subgrid[i].status == 0 || subgrid[i].status == 1) {
	my_top[i] = subgrid[i].status;
      }
      else {
	printf("ERROR ALERT: subgrid[i].status == %d\n", subgrid[i].status);
      }

      my_top[i] = subgrid[i].status;							// first row
      my_bottom[i] = subgrid[elements_per_subgrid - elements_per_row + i].status;	// last row
    }

    MPI_Send(my_top,    elements_per_row, MPI_INT, (rank-1+p)%p, 0, MPI_COMM_WORLD);	// previous rank
    MPI_Send(my_bottom, elements_per_row, MPI_INT, (rank+1)%p,   0, MPI_COMM_WORLD);	// next rank
    printf("[%d] MPI Sent top and bottom\n", rank);

    // receive my top and bottom neighbor's subgrids
    MPI_Recv(bottom, elements_per_row, MPI_INT,(rank+1)%p,MPI_ANY_TAG,MPI_COMM_WORLD,&status);  // next rank, first row
    MPI_Recv(top, elements_per_row, MPI_INT,(rank-1+p)%p,MPI_ANY_TAG,MPI_COMM_WORLD,&status);	// previous rank, last row

    printf("[%d] MPI Received top/bottom\n", rank);

    //for (j = 0; j < elements_per_row; j++) {
    //    printf("[%d] Testing top[%d]: %d\n",rank, j, top[j]);
    //}


    printf("[%d] Barrier up #2\n", rank);
    // block until everybody is ready
    MPI_Barrier(MPI_COMM_WORLD);
    printf("[%d] Barrier down #2\n", rank);


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
}

int AnalyzeNeighbors(int rank, Cell* subgrid, Cell* cell, int index, int* top, int* bottom) {
    /* 
    if (rank == 0) {
    	printf("[%d] Making call to FindEast(subgrid, index=%d, top={%d,%d,%d,%d,%d,%d,%d,%d}\n",
         	  rank, index, top[0], top[1], top[2], top[3], top[4], top[5], top[6], top[7]);
    	int c  = FindEast(subgrid, index);
    	printf("[%d] FINDEAST returned int = %d\n", rank, c);
    }
    

    cell->N->status  = FindNorth(subgrid, index, top);
    printf("Finished N\n");
    cell->NE->status = FindNorthEast(subgrid, index, top);
    printf("Finished NE\n");
    cell->E->status  = FindEast(subgrid, index);
    printf("Finished E\n");
    cell->SE->status = FindSouthEast(subgrid, index, bottom);
    printf("Finished SE\n");
    cell->S->status  = FindSouth(subgrid, index, bottom);
    cell->SW->status = FindSouthWest(subgrid, index, bottom);
    printf("Finished SW\n");
    cell->W->status  = FindWest(subgrid, index);
    cell->NW->status = FindNorthWest(subgrid, index, top);
    printf("\n---Finished all Finds---\n");
    */ 
    int sum = 0;
    int k;

    printf("\n---[%d::%d] --- ", rank, index);
    k = FindNorth(subgrid, index, top);
    if (k > 0) { printf("North --- "); sum+=k; k=0;}
    k = FindNorthEast(subgrid, index, top);
    if (k > 0) { printf("NorthEast --- "); sum+=k; k=0; }
    k = FindEast(subgrid, index);
    if (k > 0) { printf("East --- "); sum+=k; k=0; }
    k = FindSouthEast(subgrid, index, bottom);
    if (k > 0) { printf("SouthEast --- "); sum+=k; k=0; }
    k = FindSouth(subgrid, index, bottom);
    if (k > 0) { printf("South --- "); sum+=k; k=0; }
    k = FindSouthWest(subgrid, index, bottom);
    if (k > 0) { printf("SouthWest --- "); sum+=k; k=0; }
    k = FindWest(subgrid, index);
    if (k > 0) { printf("West --- "); sum+=k; k=0; }
    k = FindNorthWest(subgrid, index, top);
    if (k > 0) { printf("NorthWest --- "); sum+=k; k=0; }

    printf("---[%d::%d]---Sum = %d---\n", rank,index, sum);

    return sum;
}


// FIRST CALL:
// subgrid = [ cell0, cell1, ...]
// i = 0
// top = {1,0,1,0,0,0,0,0}
// return = 1
int FindNorth(Cell* subgrid, int i, int* top) {
    int temp;
//      0 < 7   --> TRUE
    if (i < elements_per_row) {     // at north border
        //     top[
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

    //printf("INSIDE NE\n");

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

	printf("rank#%d", rank);
	int i;
        for (i = 0; i < elements_per_subgrid; i++) {
            if (i%elements_per_row==0){
                printf("\n");
            }

	    if (subgrid[i].status == ALIVE) { printf("X "); }
	    else if (subgrid[i].status == DEAD) { printf("_ "); }
	    else { printf("? "); }

            //subgrid[i].status == ALIVE ? printf("X ", 254) : printf("_ ");
        }

    printf("\n");
    return;
}

