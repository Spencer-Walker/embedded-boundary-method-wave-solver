// Some extra information
static char help[] = "Solves the wave equation using some embedded boundary method\n\n";
/*
  u_tt - \Delta u_old = 0

  which we split into two first-order systems

  u_t -     u    = 0 
  v_t - \Delta u_old = 0 
*/

// Include some petsc stuff
#include <petscdmda.h>
#include <petscts.h>
#include <iostream>
#include <iomanip>
// Function that makes each processor print in order 
extern void printInOrder(PetscMPIInt rank, PetscMPIInt size, const char* fmt, ...);
extern PetscErrorCode TimeStep(DM da, PetscReal Lx, PetscReal Ly, Vec X, PetscReal dt2);
extern PetscErrorCode FormInitialSolution(DM,Vec, PetscReal, PetscReal, PetscReal);
extern PetscErrorCode MyMonitor(TS,PetscInt,PetscReal,Vec,void*);

// Main
int main(int argc,char **argv)
{
  //---------------------------------------------------------------------------
  // Number of global gridpoints in the x and y direction
  PetscInt       Mx = 50, My = 50;
  // MPI rank and size
  PetscMPIInt    rank,size;
  // PETSC error code
  PetscErrorCode ierr;
  // Discritization manager/distributed array
  DM             da;
  // Solution and residual vectors 
  Vec            x;                      
  // PETSC viewer for exporting/printing complicated data structures
  PetscViewer    viewer;
  // 
  PetscReal      dt = 0.005, Lx = 1.0, Ly = 1.0;
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
  
  // Extra command line arguments for the global gridpoints 
  ierr = PetscOptionsGetInt(NULL,NULL,"-Mx",&Mx,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(NULL,NULL,"-My",&My,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-dt",&dt,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-Lx",&Lx,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-Ly",&Ly,NULL);CHKERRQ(ierr);
  
  // Only rank 0 prints things (just gridpoints for now)
  if (rank == 0) 
     PetscPrintf(PETSC_COMM_SELF,"[%d/%d] Mx = %D, My = %D \n",rank,size,Mx,My); 

  // Create the 2D grid that we will use for this problem
  DMDACreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DMDA_STENCIL_STAR, \
    Mx, My, PETSC_DECIDE, PETSC_DECIDE, 2, 1, PETSC_NULL, PETSC_NULL, &da);
  DMSetFromOptions(da);
  DMSetUp(da);
  DMDASetFieldName(da,0,"u_old");
  DMDASetFieldName(da,1,"u");

  // Create solution and residual vectors
  DMCreateGlobalVector(da,&x);  

  // Set initial condition 
  FormInitialSolution(da,x,Lx,Ly,dt);
  
  // Norm
  PetscReal norm;
  ierr = VecNorm(x,NORM_2,&norm);CHKERRQ(ierr);
  if (rank == 0)
    std::cout << "Initial Norm = " << norm << std::endl;

  PetscReal dt2 = dt*dt;

  for(int i = 0; i< 500; i++)
  {
    // Timestep 
    TimeStep(da,Lx,Ly,x,dt2);

    ierr = VecNorm(x,NORM_2,&norm);CHKERRQ(ierr);
    if (rank == 0)
      std::cout << "Norm = " << std::setprecision (16) << norm << std::endl;
    std::string s = "wave"+std::to_string(i)+".out";
    PetscViewerASCIIOpen(PETSC_COMM_WORLD,s.c_str(),&viewer);
    VecView(x,viewer);

  }


  // Free space 
  VecDestroy(&x);
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
#define __FUNCT__ "TimeStep"
/*
   TimeStep - Evaluates nonlinear function, F(x).

   Input Parameters:
.  ts - the TS context
.  X - input vector
.  ptr - optional user-defined context, as set by SNESSetFunction()

   Output Parameter:
.  F - function vector
 */
PetscErrorCode TimeStep(DM da, PetscReal Lx, PetscReal Ly, Vec X, PetscReal dt2)
{
  PetscErrorCode ierr;
  PetscInt       i,j,Mx,My,xs,ys,xm,ym;
  PetscReal      hx,hy,/*hxdhy,hydhx,*/ sx,sy;
  PetscScalar    u_old,uxx,uyy,u,***x, tmp_x;
  Vec            localX;

  PetscFunctionBeginUser;
  ierr = DMGetLocalVector(da,&localX);CHKERRQ(ierr);
  ierr = DMDAGetInfo(da,PETSC_IGNORE,&Mx,&My,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,
                     PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE);

  hx = Lx/(PetscReal)(Mx-1); sx = 1.0/(hx*hx);
  hy = Ly/(PetscReal)(My-1); sy = 1.0/(hy*hy);
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
        continue;
      }
      u_old      = x[j][i][0];
      u          = x[j][i][1];
      uxx        = (-2.0*u + x[j][i-1][1] + x[j][i+1][1])*sx;
      uyy        = (-2.0*u + x[j-1][i][1] + x[j+1][i][1])*sy;
      x[j][i][0] = dt2*(uxx + uyy)+2.0*u-u_old;
    }
  }

  for (j=ys; j<ys+ym; j++) {
    for (i=xs; i<xs+xm; i++) {
      tmp_x = x[j][i][0];
      x[j][i][0] = x[j][i][1];
      x[j][i][1] = tmp_x;
    }
  }


  /*
     Restore vectors
  */
  ierr = DMDAVecRestoreArrayDOF(da,localX,&x);CHKERRQ(ierr);
  DMLocalToGlobalBegin(da,localX,INSERT_VALUES,X);
  DMLocalToGlobalEnd(da,localX,INSERT_VALUES,X);
  ierr = DMRestoreLocalVector(da,&localX);CHKERRQ(ierr);



  PetscFunctionReturn(0);
}


/* ------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "FormInitialSolution"
PetscErrorCode FormInitialSolution(DM da,Vec U, PetscReal Lx, PetscReal Ly, PetscReal dt)
{
  PetscErrorCode ierr;
  PetscInt       i,j,xs,ys,xm,ym,Mx,My;
  PetscScalar    ***u;
  PetscReal      hx,hy,x,y,r;

  PetscFunctionBeginUser;
  ierr = DMDAGetInfo(da,PETSC_IGNORE,&Mx,&My,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,
                     PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE);

  hx = Lx/(PetscReal)(Mx-1);
  hy = Lx/(PetscReal)(My-1);

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
      r = PetscSqrtScalar((x-.5*Lx)*(x-.5*Lx) + (y-.5*Ly)*(y-.5*Ly));
      if (r < Lx || r < Ly) {
        u[j][i][0] = PetscExpScalar(-30.0*r*r);
        u[j][i][1] = PetscExpScalar(-30.0*r*r)+ 0.0*dt;
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

