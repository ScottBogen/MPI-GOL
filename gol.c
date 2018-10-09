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


Cell* FindNorth(subgrid, cell, index, top);
Cell* FindNorthEast(subgrid, cell, index, top);
Cell* FindEast(subgrid, cell, index);
Cell* FindSouthEast(subgrid, cell, index, bottom);
Cell* FindSouth(subgrid, cell, index, bottom);
Cell* FindSouthWest(subgrid, cell, index, bottom);
Cell* FindWest(subgrid, cell, index);
Cell* FindNorthWest(subgrid, cell, index, top);

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
    DisplayGoL(rank, p, n, subgrid);

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

            else {
                subgrid[i].status = ALIVE;
                subgrid[i].next_status = ALIVE;
            }

        }
    }
}

void Simulate(int rank, int p, int n, Cell* subgrid) {
    int i;

    MPI_Status status;
    int* my_top = (int*) malloc(sizeof(my_top) * elements_per_row);
    int* my_bottom = (int*) malloc(sizeof(my_bottom) * elements_per_row);
    int* top = (int*) malloc(sizeof(top) * elements_per_row);
    int* bottom = (int*) malloc(sizeof(bottom) * elements_per_row);

    Cell* temp;
    temp = subgrid;

    memset(my_top, 50, elements_per_row*sizeof(my_top[0]));
    memset(my_bottom, 5, elements_per_row*sizeof(my_bottom[0]));
    memset(top, 10, elements_per_row*sizeof(top[0]));
    memset(bottom, 20, elements_per_row*sizeof(bottom[0]));

    MPI_Barrier(MPI_COMM_WORLD);

    for (i = 0; i < elements_per_row; i++) {
      my_top[i] = subgrid[i].status;
      my_bottom[i] = subgrid[elements_per_subgrid - elements_per_row + i].status;
    }

    int j;
    //for (j = 0; j < elements_per_row; j++) {
    //    printf("[%d] Testing my top[%d].status: %d\n", rank, j, my_top[j]);
    //}

    MPI_Send(my_top,    elements_per_row, MPI_INT, (rank-1+p)%p, 0, MPI_COMM_WORLD);
    MPI_Send(my_bottom, elements_per_row, MPI_INT, (rank+1)%p,   0, MPI_COMM_WORLD);
    printf("[%d] MPI Sent top and bottom\n", rank);

    //memset(top, 10, elements_per_subgrid*sizeof(top[0]));
    //memset(bottom, 20, elements_per_subgrid*sizeof(bottom[0]));

    // receive my top and bottom neighbor's subgrids
    MPI_Recv(top, elements_per_row, MPI_INT,(rank-1+p)%p,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
    MPI_Recv(bottom, elements_per_row, MPI_INT,(rank+1)%p,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
    printf("[%d] MPI Received top/bottom\n", rank);

    for (j = 0; j < elements_per_row; j++) {
        printf("[%d] Testing top[%d]: %d\n",rank, j, top[j]);
    }
    printf("[%d] Barrier up #2\n", rank);
    // block until everybody is ready
    MPI_Barrier(MPI_COMM_WORLD);
    printf("[%d] Barrier down #2\n", rank);

    /*
    for (i = 0; i < n*n/p; i++) {
      int sum = AnalyzeNeighbors(subgrid, subgrid[i], i, top, bottom);
      if (sum < 3 || sum > 5) {
        subgrid[i].next_status = DEAD;
      }
      else {
        subgrid[i].next_status = ALIVE;
      }
    }
    */
}

int AnalyzeNeighbors(Cell* subgrid, Cell* cell, int index, Cell* top, Cell* bottom) {
	// todo: remove cell from each call, it's not needed
    cell->N  = FindNorth(subgrid, cell, index, top);
    cell->NE = FindNorthEast(subgrid, cell, index, top);
    cell->E  = FindEast(subgrid, cell, index);
    cell->SE = FindSouthEast(subgrid, cell, index, bottom);
    cell->S  = FindSouth(subgrid, cell, index, bottom);
    cell->SW = FindSouthWest(subgrid, cell, index, bottom);
    cell->W  = FindWest(subgrid, cell, index);
    cell->NW = FindNorthWest(subgrid, cell, index, top);


    int sum = 0;
    sum += cell->N->status;
    sum += cell->NE->status;
    sum += cell->E->status;
    sum += cell->SE->status;
    sum += cell->S->status;
    sum += cell->SW->status;
    sum += cell->W->status;
    sum += cell->NW->status;

    return sum;
}

Cell* FindNorth(Cell* subgrid, Cell* c, int i, Cell* top) {
    Cell* temp;
    if (i < elements_per_row) {     // at north border
        temp = &top[elements_per_subgrid - (elements_per_row - 1 - i)];
    }
    else {
        temp = &subgrid[i - elements_per_row];
    }
    return temp;
}

Cell* FindNorthEast(Cell* subgrid, Cell* c, int i, Cell* top) {
    Cell* temp;
    int temp_index;

    if (i < elements_per_row) {     // at north border
        temp_index = elements_per_subgrid - (elements_per_row - 1 - i);
        if ((temp_index+1) % elements_per_row == 0) { // at east border
            temp = &top[temp_index-elements_per_row+1];
        }
        else {
            temp = &top[temp_index+1];
        }
    }

    else {
        temp_index = i - elements_per_row + 1;  // normal move north

        if ((temp_index+1) % elements_per_row == 0) { // at east border
            temp = &subgrid[temp_index-elements_per_row+1];
        }
        else {  // normal case
            temp = &subgrid[temp_index+1];
        }

    }

    return temp;
}

Cell* FindEast(Cell* subgrid, Cell* c, int i) {
    Cell* temp;
    if ((i + 1) % elements_per_row == 0) {  // at east border
        temp = &subgrid[i - elements_per_row + 1];
    }
    else {      // normal case
        temp = &subgrid[i+1];
    }
    return temp;
}

Cell* FindSouthEast(Cell* subgrid, Cell* c, int i, Cell* bottom) {
    Cell* temp;
    int temp_index;

    // SOUTH
    if (i >= elements_per_subgrid - elements_per_row) {     // at south border
        temp_index = i % elements_per_row;

        if ((temp_index+1) % elements_per_row == 0) { // at east border
            temp = &bottom[temp_index-elements_per_row+1];
        }
        else {
            temp = &bottom[temp_index+1];
        }
    }

    else {
        temp_index = i + elements_per_row;  // normal move south


        if ((temp_index+1) % elements_per_row == 0) { // at east border
            temp = &subgrid[temp_index-elements_per_row+1];  // wrap around
        }
        else {  // normal case
            temp = &subgrid[temp_index+1];
        }

    }

    return temp;
}

Cell* FindSouth(Cell* subgrid, Cell* c, int i, Cell* bottom) {
    Cell* temp;
    if (i >= elements_per_subgrid - elements_per_row) {     // at north border
        temp = &bottom[i%elements_per_row];
    }
    else {
        temp = &subgrid[i + elements_per_row];
    }
    return temp;
}

Cell* FindSouthWest(Cell* subgrid, Cell* c, int i, Cell* bottom) {
  Cell* temp;
  int temp_index;

  // SOUTH
  if (i >= elements_per_subgrid - elements_per_row) {     // at south border
      temp_index = i % elements_per_row;

      if ((temp_index+1) % elements_per_row == 0) { // at east border
          temp = &bottom[temp_index-elements_per_row+1];
      }
      else {
          temp = &bottom[temp_index+1];
      }
  }

  else {
      temp_index = i + elements_per_row;  // normal move south


      // WEST
      if (temp_index % elements_per_row == 0) {    // at west border
          temp = &subgrid[temp_index + elements_per_row - 1];
      }
      else {
          temp = &subgrid[temp_index - 1];
      }

  }

  return temp;
}

Cell* FindWest(Cell* subgrid, Cell* c, int i) {
    Cell* temp;

    if (i % elements_per_row == 0) {    // at west border
        temp = &subgrid[i + elements_per_row - 1];
    }
    else {
        temp = &subgrid[i - 1];
    }
    return temp;
}

Cell* FindNorthWest(Cell* subgrid, Cell* c, int i, Cell* top) {
    Cell* temp;
    int temp_index;

    if (i < elements_per_row) {     // at north border
        temp_index = elements_per_subgrid - (elements_per_row - 1 - i);
        if ((temp_index+1) % elements_per_row == 0) { // at east border
            temp = &top[temp_index-elements_per_row+1];
        }
        else {
            temp = &top[temp_index+1];
        }
    }

    else {
        temp_index = i - elements_per_row + 1;  // normal move north

        if (temp_index % elements_per_row == 0) {    // at west border
            temp = &subgrid[temp_index + elements_per_row - 1];
        }
        else {
            temp = &subgrid[temp_index - 1];
        }
    }
}

void DisplayGoL(int rank, int p, int n, Cell* subgrid) {

	printf("rank#%d", rank);
	int i;
        for (i = 0; i < n*n/p; i++) {
            if (i%n==0){
                printf("\n");
            }

            subgrid[i].status == ALIVE ? printf("X ", 254) : printf("_ ");
        }

    printf("\n");
    return;
}

