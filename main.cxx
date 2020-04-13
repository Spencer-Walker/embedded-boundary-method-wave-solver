// Include some petsc stuff
#include <petscdmda.h>
// Some extra information
static char help[] = "Solves the wave equation using some embedded boundary method\n\n";
// Main
int main(int argc,char **argv)
{
  //---------------------------------------------------------------------------
  // Number of global gridpoints in the x and y direction
  PetscInt       Nx = 5, Ny = 4;
  // MPI rank and size
  PetscMPIInt    rank,size,message;
  // PETSC error code
  PetscErrorCode ierr;
  // Discritization manager/distributed array
  DM             da;
  // Local and global solution vectors
  Vec            local,global;
  // Some real number 
  PetscScalar    value;
  // PETSC viewer for exporting/printing complicated data structures
  PetscViewer    viewer;
  //-------------------------------------------------------------------------------------
  // Setup PETSC and MPI 
  //-------------------------------------------------------------------------------------
  // Initialize petsc (and MPI)
  ierr = PetscInitialize(&argc,&argv,(char*)0,help);if (ierr) return ierr;

  // MPI rank and size
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
  MPI_Comm_size(PETSC_COMM_WORLD,&size);

  // Print from each rank in order (useful for when we eventually want to look at output)
  if (rank == 0) 
  {
    MPI_Send(&message, 1, MPI_INT, 1, 0, PETSC_COMM_WORLD);
    PetscPrintf(PETSC_COMM_SELF,"HI FROM RANK = %d of SIZE = %d \n", rank, size);
  } 
  else 
  {
    int buffer;
    MPI_Status status;
    MPI_Probe(MPI_ANY_SOURCE, 0, PETSC_COMM_WORLD, &status);
    MPI_Get_count(&status, MPI_INT, &buffer);
    if (buffer == 1) 
    {
      PetscPrintf(PETSC_COMM_SELF,"HI FROM RANK = %d of SIZE = %d \n", rank, size);
      MPI_Recv(&message, buffer, MPI_INT, MPI_ANY_SOURCE, 0, PETSC_COMM_WORLD, &status);
      if (rank + 1 != size) 
      {
        MPI_Send(&message, 1, MPI_INT, ++rank, 0, PETSC_COMM_WORLD);
      }
    }
  }

  // Create a viewer object that prints ascii to terminal
  PetscViewerCreate(PETSC_COMM_WORLD,&viewer);
  PetscViewerSetType(viewer,PETSCVIEWERASCII);	
  
  // Extra command line arguments for the global gridpoints 
  ierr = PetscOptionsGetInt(NULL,NULL,"-Nx",&Nx,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(NULL,NULL,"-Ny",&Ny,NULL);CHKERRQ(ierr);
  
  // Only rank 0 prints things (just gridpoints for now)
  if (rank == 0) 
     PetscPrintf(PETSC_COMM_SELF,"[%d/%d] Nx = %D, Ny = %D \n",rank,size,Nx,Ny); 

  // Create the 2D grid that we will use for this problem
  DMDACreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DMDA_STENCIL_STAR, \
    Nx, Ny, PETSC_DECIDE, PETSC_DECIDE, 1, 1, PETSC_NULL, PETSC_NULL, &da);
  DMSetFromOptions(da);
  DMSetUp(da);
  DMCreateGlobalVector(da,&global);
  DMCreateLocalVector(da,&local);

  value = rank;
  VecSet(local,value);
  DMLocalToGlobalBegin(da,local,INSERT_VALUES,global);
  DMLocalToGlobalEnd(da,local,INSERT_VALUES,global);

  DMView(da,viewer);
  
  VecDestroy(&local);
  VecDestroy(&global);
  DMDestroy(&da);

  // Close PETSC (and MPI)
  ierr = PetscFinalize();
  return ierr;
}



