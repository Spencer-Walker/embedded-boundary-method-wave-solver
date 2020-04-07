#include <petscvec.h>
static char help[] = "Solves the wave equation using some embedded boundary method\n\n";
int main(int argc,char **argv)
{
  PetscInt       n = 20;
  PetscMPIInt    rank,size;
  PetscErrorCode ierr;

  ierr = PetscInitialize(&argc,&argv,(char*)0,help);if (ierr) return ierr;
  ierr = PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL);CHKERRQ(ierr);
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
  MPI_Comm_size(PETSC_COMM_WORLD,&size);
  if (rank == 0) 
 		PetscPrintf(PETSC_COMM_SELF,"[%d] n = %D \n",rank,n); 
	

  ierr = PetscFinalize();
  return ierr;
}



