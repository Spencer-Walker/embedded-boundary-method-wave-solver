// Include some petsc stuff
#include <petscvec.h>
#include <petscksp.h>
#include <petscdm.h>
#include <petscdmda.h>
// Some extra information
static char help[] = "Solves the wave equation using some embedded boundary method\n\n";
// Main
int main(int argc,char **argv)
{
	//---------------------------------------------------------------------------
	// Number of global gridpoints in the x and y direction
  PetscInt       Nx = 10, Ny = 8;
	// MPI rank and size
  PetscMPIInt    rank,size;
  // PETSC error code
	PetscErrorCode ierr;
	// Discritization manager
	DM da;
	// Local and global solution vectors
	Vec local,global;
	// 
	PetscScalar value;
	//
	PetscViewer viewer;
	//----------------------------------------------------------------------------

	// Initialize petsc (and MPI)
  ierr = PetscInitialize(&argc,&argv,(char*)0,help);if (ierr) return ierr;
  PetscViewerCreate(PETSC_COMM_WORLD,&viewer);
	PetscViewerSetType(viewer,PETSCVIEWERASCII);	
	
  // Extra command line arguments for the global gridpoints 
	ierr = PetscOptionsGetInt(NULL,NULL,"-Nx",&Nx,NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetInt(NULL,NULL,"-Ny",&Ny,NULL);CHKERRQ(ierr);
	
	// MPI rank and size
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
  MPI_Comm_size(PETSC_COMM_WORLD,&size);

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

	value = -3;
	VecSet(global,value);
	DMGlobalToLocalBegin(da,global,ADD_VALUES,local);
	DMGlobalToLocalEnd(da,global,ADD_VALUES,local);
	
	value = rank+1;
	VecScale(local,value);
	DMLocalToGlobalBegin(da,local,ADD_VALUES,global);
	DMLocalToGlobalEnd(da,local,ADD_VALUES,global);

	DMView(da,viewer);
	VecView(global,viewer);


	PetscViewerDestroy(&viewer);
	VecDestroy(&local);
	VecDestroy(&global);
	DMDestroy(&da);

  // Close PETSC (and MPI)
	ierr = PetscFinalize();
  return ierr;
}



