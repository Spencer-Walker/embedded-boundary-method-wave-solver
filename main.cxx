// Some extra information
static char help[] = "Solves the wave equation using some embedded boundary method\n\n";
/*
  u_tt - \Delta u = 0

  which we split into two first-order systems

  u_t -     v    = 0 
  v_t - \Delta u = 0 
*/

// Include some petsc stuff
#include <petscdmda.h>
#include <petscts.h>
// Function that makes each processor print in order 
extern void printInOrder(PetscMPIInt rank, PetscMPIInt size, const char* fmt, ...);
extern PetscErrorCode FormFunction(TS,PetscReal,Vec,Vec,void*),FormInitialSolution(DM,Vec);
extern PetscErrorCode MyTSMonitor(TS,PetscInt,PetscReal,Vec,void*);

// Main
int main(int argc,char **argv)
{
  //---------------------------------------------------------------------------
  // Number of global gridpoints in the x and y direction
  PetscInt       Nx = 5, Ny = 4;
  // MPI rank and size
  PetscMPIInt    rank,size;
  // PETSC error code
  PetscErrorCode ierr;
  // Discritization manager/distributed array
  DM             da;
  // Timestepping context
  TS             ts;
  // Solution and residual vectors 
  Vec            x,r;                        
  // PETSC viewer for exporting/printing complicated data structures
  PetscViewer    viewer;
  PetscInt       steps;    
  SNES           ts_snes;
  PetscReal      ftime;   
  //-------------------------------------------------------------------------------------
  // Setup PETSC and MPI 
  //-------------------------------------------------------------------------------------
  // Initialize petsc (and MPI)
  ierr = PetscInitialize(&argc,&argv,(char*)0,help);if (ierr) return ierr;

  // MPI rank and size
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
  MPI_Comm_size(PETSC_COMM_WORLD,&size);

  // Print hi in order
  printInOrder(rank,size,"HI FROM RANK = %d of SIZE = %d \n", rank, size);

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
    Nx, Ny, PETSC_DECIDE, PETSC_DECIDE, 2, 1, PETSC_NULL, PETSC_NULL, &da);
  DMSetFromOptions(da);
  DMSetUp(da);
  DMDASetFieldName(da,0,"u");
  DMDASetFieldName(da,1,"v");


  // Create solution and residual vectors
  DMCreateGlobalVector(da,&x);
  VecDuplicate(x,&r);

  // Create timestepping solver context
  TSCreate(PETSC_COMM_WORLD,&ts);
  TSSetDM(ts,da);
  TSSetProblemType(ts,TS_NONLINEAR);
  TSSetRHSFunction(ts,NULL,FormFunction,da);
  TSSetMaxTime(ts,1.0);
  TSMonitorSet(ts,MyTSMonitor,0,0);
  TSSetExactFinalTime(ts,TS_EXACTFINALTIME_STEPOVER);
  // Customize solver 
  TSSetType(ts,TSBEULER);
  TSGetSNES(ts,&ts_snes);

  // Set initial conditions
  FormInitialSolution(da,x);
  TSSetTime(ts,0.0);
  TSSetTimeStep(ts,0.0001);
  TSSetSolution(ts,x);

  // Set runtime options
  TSSetFromOptions(ts);

  // Solve the system
  TSSolve(ts,x);
  TSGetSolveTime(ts,&ftime);
  TSGetStepNumber(ts,&steps);

  // Free space 
  VecDestroy(&x);
  VecDestroy(&r);
  TSDestroy(&ts);
  DMDestroy(&da);

  // Close PETSC (and MPI)
  PetscFinalize();
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "printInOrder"
void printInOrder(PetscMPIInt rank, PetscMPIInt size, const char* fmt, ...)
{
  PetscMPIInt message;
  va_list args;
  va_start(args,fmt); 
  if (rank == 0) 
  {
    MPI_Send(&message, 1, MPI_INT, 1, 0, PETSC_COMM_WORLD);
    vprintf(fmt,args);
  } 
  else 
  {
    int buffer;
    MPI_Status status;
    MPI_Probe(MPI_ANY_SOURCE, 0, PETSC_COMM_WORLD, &status);
    MPI_Get_count(&status, MPI_INT, &buffer);
    if (buffer == 1) 
    {
      vprintf(fmt,args);
      MPI_Recv(&message, buffer, MPI_INT, MPI_ANY_SOURCE, 0, PETSC_COMM_WORLD, &status);
      if (rank + 1 != size) 
      {
        MPI_Send(&message, 1, MPI_INT, ++rank, 0, PETSC_COMM_WORLD);
      }
    }
  }
  va_end(args);
}

/* ------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "FormFunction"
/*
   FormFunction - Evaluates nonlinear function, F(x).

   Input Parameters:
.  ts - the TS context
.  X - input vector
.  ptr - optional user-defined context, as set by SNESSetFunction()

   Output Parameter:
.  F - function vector
 */
PetscErrorCode FormFunction(TS ts,PetscReal ftime,Vec X,Vec F,void *ptr)
{
  DM             da = (DM)ptr;
  PetscErrorCode ierr;
  PetscInt       i,j,Mx,My,xs,ys,xm,ym;
  PetscReal      hx,hy,/*hxdhy,hydhx,*/ sx,sy;
  PetscScalar    u,uxx,uyy,v,***x,***f;
  Vec            localX;

  PetscFunctionBeginUser;
  ierr = DMGetLocalVector(da,&localX);CHKERRQ(ierr);
  ierr = DMDAGetInfo(da,PETSC_IGNORE,&Mx,&My,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,
                     PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE);

  hx = 1.0/(PetscReal)(Mx-1); sx = 1.0/(hx*hx);
  hy = 1.0/(PetscReal)(My-1); sy = 1.0/(hy*hy);
  /*hxdhy  = hx/hy;*/
  /*hydhx  = hy/hx;*/

  /*
     Scatter ghost points to local vector,using the 2-step process
        DMGlobalToLocalBegin(),DMGlobalToLocalEnd().
     By placing code between these two statements, computations can be
     done while messages are in transition.
  */
  ierr = DMGlobalToLocalBegin(da,X,INSERT_VALUES,localX);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(da,X,INSERT_VALUES,localX);CHKERRQ(ierr);

  /*
     Get pointers to vector data
  */
  ierr = DMDAVecGetArrayDOF(da,localX,&x);CHKERRQ(ierr);
  ierr = DMDAVecGetArrayDOF(da,F,&f);CHKERRQ(ierr);

  /*
     Get local grid boundaries
  */
  ierr = DMDAGetCorners(da,&xs,&ys,NULL,&xm,&ym,NULL);CHKERRQ(ierr);

  /*
     Compute function over the locally owned part of the grid
  */
  for (j=ys; j<ys+ym; j++) {
    for (i=xs; i<xs+xm; i++) {
      if (i == 0 || j == 0 || i == Mx-1 || j == My-1) {
        f[j][i][0] = x[j][i][0];
        f[j][i][1] = x[j][i][1];
        continue;
      }
      u          = x[j][i][0];
      v          = x[j][i][1];
      uxx        = (-2.0*u + x[j][i-1][0] + x[j][i+1][0])*sx;
      uyy        = (-2.0*u + x[j-1][i][0] + x[j+1][i][0])*sy;
      f[j][i][0] = v;
      f[j][i][1] = uxx + uyy;
    }
  }

  /*
     Restore vectors
  */
  ierr = DMDAVecRestoreArrayDOF(da,localX,&x);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArrayDOF(da,F,&f);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(da,&localX);CHKERRQ(ierr);
  ierr = PetscLogFlops(11.0*ym*xm);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/* ------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "FormInitialSolution"
PetscErrorCode FormInitialSolution(DM da,Vec U)
{
  PetscErrorCode ierr;
  PetscInt       i,j,xs,ys,xm,ym,Mx,My;
  PetscScalar    ***u;
  PetscReal      hx,hy,x,y,r;

  PetscFunctionBeginUser;
  ierr = DMDAGetInfo(da,PETSC_IGNORE,&Mx,&My,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,
                     PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE);

  hx = 1.0/(PetscReal)(Mx-1);
  hy = 1.0/(PetscReal)(My-1);

  /*
     Get pointers to vector data
  */
  ierr = DMDAVecGetArrayDOF(da,U,&u);CHKERRQ(ierr);

  /*
     Get local grid boundaries
  */
  ierr = DMDAGetCorners(da,&xs,&ys,NULL,&xm,&ym,NULL);CHKERRQ(ierr);

  /*
     Compute function over the locally owned part of the grid
  */
  for (j=ys; j<ys+ym; j++) {
    y = j*hy;
    for (i=xs; i<xs+xm; i++) {
      x = i*hx;
      r = PetscSqrtScalar((x-.5)*(x-.5) + (y-.5)*(y-.5));
      if (r < .125) {
        u[j][i][0] = PetscExpScalar(-30.0*r*r*r);
        u[j][i][1] = 0.0;
      } else {
        u[j][i][0] = 0.0;
        u[j][i][1] = 0.0;
      }
    }
  }

  /*
     Restore vectors
  */
  ierr = DMDAVecRestoreArrayDOF(da,U,&u);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MyTSMonitor"
PetscErrorCode MyTSMonitor(TS ts,PetscInt step,PetscReal ptime,Vec v,void *ctx)
{
  PetscErrorCode ierr;
  PetscReal      norm;
  MPI_Comm       comm;

  PetscFunctionBeginUser;
  ierr = VecNorm(v,NORM_2,&norm);CHKERRQ(ierr);
  ierr = PetscObjectGetComm((PetscObject)ts,&comm);CHKERRQ(ierr);
  ierr = PetscPrintf(comm,"timestep %D time %g norm %g\n",step,ptime,norm);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

