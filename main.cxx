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
externPetscErrorCode TimeStep(DM da_u_old,DM da_u,DM da_phi1_old,DM da_phi1,DM da_phi2_old,DM da_phi2, \
  PetscReal Lx, PetscReal Ly, Vec vec_u_old, Vec vec_u, Vec vec_phi1_old, Vec vec_phi1, \
  Vec vec_phi2_old, Vec vec_phi2, PetscReal dt, PetscReal pml_width, PetscReal pml_strength, \
  PetscReal t, PetscReal kx, PetscReal ky);
PetscErrorCode FormInitialSolution( DM da_u_old, DM da_u, DM da_phi1_old, DM da_phi1, DM da_phi2_old, DM da_phi2\
  , Vec vec_u_old, Vec vec_u, Vec vec_phi1_old, Vec vec_phi1, Vec vec_phi2_old, Vec vec_phi2 , \
  PetscReal Lx, PetscReal Ly, PetscReal dt,PetscReal width);
extern PetscErrorCode MyMonitor(TS,PetscInt,PetscReal,Vec,void*);
extern PetscReal PML(PetscReal Lx, PetscReal x, PetscReal width, PetscReal amp);

// Main
int main(int argc,char **argv)
{
  //---------------------------------------------------------------------------
  // Number of global gridpoints in the x and y direction
  PetscInt       Mx = 50, My = 50, Nt = 10;
  // MPI rank and size
  PetscMPIInt    rank,size;
  // PETSC error code
  PetscErrorCode ierr;
  // Discritization manager/distributed array
  DM             da_u_old, da_u, da_phi1_old, da_phi1, da_phi2_old, da_phi2;
  // Solution and residual vectors 
  Vec            vec_u_old, vec_u, vec_phi1_old, vec_phi1, vec_phi2_old, vec_phi2;                      
  // PETSC viewer for exporting/printing complicated data structures
  PetscViewer    viewer;
  // 
  PetscReal      dt = 0.005, Lx = 1.0, Ly = 1.0, pml_width = 0.0, pml_strength = 0.0;
  PetscReal      kx = 100.0, ky = 0.0;
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
  ierr = PetscOptionsGetInt(NULL,NULL,"-Nt",&Nt,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-dt",&dt,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-Lx",&Lx,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-Ly",&Ly,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-pml_width",&pml_width,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-pml_strength",&pml_strength,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-kx",&kx,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-ky",&ky,NULL);CHKERRQ(ierr);
  // Only rank 0 prints things (just gridpoints for now)
  if (rank == 0) 
     PetscPrintf(PETSC_COMM_SELF,"[%d/%d] Mx = %D, My = %D \n",rank,size,Mx,My); 

  // Create the 2D grid that we will use for this problem
  DMDACreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DMDA_STENCIL_STAR, \
    Mx, My, PETSC_DECIDE, PETSC_DECIDE, 1, 1, PETSC_NULL, PETSC_NULL, &da_u_old);
  DMSetFromOptions(da_u_old);
  DMSetUp(da_u_old);
  DMDASetFieldName(da_u_old,0,"u_old");

  DMDACreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DMDA_STENCIL_STAR, \
    Mx, My, PETSC_DECIDE, PETSC_DECIDE, 1, 1, PETSC_NULL, PETSC_NULL, &da_u);
  DMSetFromOptions(da_u);
  DMSetUp(da_u);  
  DMDASetFieldName(da,0,"u");

  DMDACreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DMDA_STENCIL_STAR, \
    Mx, My, PETSC_DECIDE, PETSC_DECIDE, 1, 1, PETSC_NULL, PETSC_NULL, &da_phi1_old);
  DMSetFromOptions(da_phi1_old);
  DMSetUp(da_phi1_old);
  DMDASetFieldName(da_phi1_old,0,"phi1_old");

  DMDACreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DMDA_STENCIL_STAR, \
    Mx, My, PETSC_DECIDE, PETSC_DECIDE, 1, 1, PETSC_NULL, PETSC_NULL, &da_phi1);
  DMSetFromOptions(da_phi1);
  DMSetUp(da_phi1);
  DMDASetFieldName(da_phi1,0,"phi1");

  DMDACreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DMDA_STENCIL_STAR, \
    Mx, My, PETSC_DECIDE, PETSC_DECIDE, 1, 1, PETSC_NULL, PETSC_NULL, &da_phi2_old);
  DMSetFromOptions(da_phi2_old);
  DMSetUp(da_phi2_old);
  DMDASetFieldName(da_phi2_old,0,"phi2_old");

  DMDACreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DMDA_STENCIL_STAR, \
    Mx, My, PETSC_DECIDE, PETSC_DECIDE, 1, 1, PETSC_NULL, PETSC_NULL, &da_phi2);
  DMSetFromOptions(da_phi2);
  DMSetUp(da_phi2);
  DMDASetFieldName(da_phi2,0,"phi2");
  

  // Create solution and residual vectors
  DMCreateGlobalVector(da_u_old,&vec_u_old);
  DMCreateGlobalVector(da_u,&vec_u);
  DMCreateGlobalVector(da_phi1_old,&vec_phi1_old);
  DMCreateGlobalVector(da_phi1,&vec_phi1);
  DMCreateGlobalVector(da_phi2_old,&vec_phi2_old);
  DMCreateGlobalVector(da_phi2,&vec_phi2);  

  // Set initial condition 
  FormInitialSolution(da_u_old, da_u, da_phi1_old, da_phi1, da_phi2_old, da_phi2 \
    ,vec_u_old, vec_u, vec_phi1_old, vec_phi1, vec_phi2_old, vec_phi2, Lx, Ly, dt, pml_width);
  
  // Norm
  PetscReal norm;
  ierr = VecNorm(x,NORM_2,&norm);CHKERRQ(ierr);
  if (rank == 0)
    std::cout << "Initial Norm = " << norm << std::endl;

  for(int i = 0; i< Nt; i++)
  {
    // Timestep 
    TimeStep(da_u_old, da_u, da_phi1_old, da_phi1, da_phi2_old, da_phi2\
      ,Lx,Ly,vec_u_old, vec_u, vec_phi1_old, vec_phi1, vec_phi2_old, vec_phi2,dt,pml_width,pml_strength,dt*i,kx,ky);

    ierr = VecNorm(vec_u_old,NORM_2,&norm);CHKERRQ(ierr);
    if (rank == 0)
      std::cout << "Norm = " << std::setprecision (16) << norm << std::endl;
    std::string s = "wave"+std::to_string(i)+".out";
    PetscViewerASCIIOpen(PETSC_COMM_WORLD,s.c_str(),&viewer);
    VecView(vec_u_old,viewer);

  }


  // Free space 
  VecDestroy(&vec_u_old);
  VecDestroy(&vec_u);
  VecDestroy(&vec_phi1_old);
  VecDestroy(&vec_phi1);
  VecDestroy(&vec_phi2_old);
  VecDestroy(&vec_phi2);
  DMDestroy(&da_u_old);
  DMDestroy(&da_u);
  DMDestroy(&da_phi1_old);
  DMDestroy(&da_phi1);
  DMDestroy(&da_phi2_old);
  DMDestroy(&da_phi2);
  ViewerDestroy(&viewer);
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
PetscErrorCode TimeStep(DM da_u_old,DM da_u,DM da_phi1_old,DM da_phi1,DM da_phi2_old,DM da_phi2, \
  PetscReal Lx, PetscReal Ly, Vec vec_u_old, Vec vec_u, Vec vec_phi1_old, Vec vec_phi1, \
  Vec vec_phi2_old, Vec vec_phi2, PetscReal dt, PetscReal pml_width, PetscReal pml_strength, \
  PetscReal t, PetscReal kx, PetscReal ky)
{
  PetscErrorCode ierr;
  PetscInt       i,j,Mx,My,xs,ys,xm,ym;
  PetscReal      hx,hy,/*hxdhy,hydhx,*/ sx,sy;
  PetscScalar    u_old,uxx,uyy,u,***x, tmp_x;
  Vec            local_u_old, local_u, local_phi1_old, local_phi1, local_phi2_old, local_phi2;

  PetscFunctionBeginUser;
  ierr = DMGetLocalVector(da,&localX);CHKERRQ(ierr);
  ierr = DMDAGetInfo(da,PETSC_IGNORE,&Mx,&My,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,
                     PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE);

  PetscReal dt2 = dt*dt;
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
  PetscReal PMLX;
  PetscReal PMLY;
  for (j=ys; j<ys+ym; j++) {
    PMLY = PML(Ly,j*hy,pml_width,pml_strength );
    for (i=xs; i<xs+xm; i++) {
      if (i == 0 || j == 0 || i == Mx-1 || j == My-1) {
        continue;
      }

      PMLX = PML(Lx,i*hx,pml_width,pml_strength );
      u_old      = x[j][i][0];
      u          = x[j][i][1];
      uxx        = (-2.0*u + x[j][i-1][1] + x[j][i+1][1])*sx;
      uyy        = (-2.0*u + x[j-1][i][1] + x[j+1][i][1])*sy;
      if (i*hx<pml_width||Lx-i*hx<pml_width||j*hy<pml_width||Ly-j*hy<pml_width)
      {
        x[j][i][0] = (dt2*(uxx + uyy) +(2.0-dt2*PMLX*PMLY)*u+(dt*(PMLX+PMLY)-1.0)*u_old+ dt2*(x[j][i+1][2]-x[j][i+1][2])/(2*hx) + dt2*(x[j+1][i][3]-x[j-1][i][3])/(2*hy)   )/(1.0+dt*(PMLX+PMLY));
        x[j][i][2] = dt*(-PMLX*x[j][i][2] +(PMLY-PMLX)*(x[j][i+1][1] - x[j][i-1][1])/(2*hx));
        x[j][i][3] = dt*(-PMLY*x[j][i][3] +(PMLX-PMLY)*(x[j+1][i][1] - x[j-1][i][1])/(2*hy));
      }
      else
      {
        x[j][i][0] = dt2*(uxx + uyy) + 2.0*u -1.0*u_old;
      }

      // Set all values inside of sphere = 0 (TEMP UNTIL EMBEDDED BOUNDARY SET)
      if ( (i*hx-Lx/2.0)*(i*hx-Lx/2.0) + (j*hy-Ly/2.0)*(j*hy-Ly/2.0) < 0.1*0.1)
      {
        x[j][i][0] = sin(kx*i*hx+ky*hy*j-sqrt(kx*kx+ky*ky)*t);
      }
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
PetscErrorCode FormInitialSolution( DM da_u_old, DM da_u, DM da_phi1_old, DM da_phi1, DM da_phi2_old, DM da_phi2\
  , Vec vec_u_old, Vec vec_u, Vec vec_phi1_old, Vec vec_phi1, Vec vec_phi2_old, Vec vec_phi2 , \
    PetscReal Lx, PetscReal Ly, PetscReal dt,PetscReal width)
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
  for (j=ys; j<ys+ym; j++) 
  {
    y = j*hy;
    for (i=xs; i<xs+xm; i++) 
    {
      x = i*hx;
      u[j][i][0] = 0.0;
      u[j][i][1] = 0.0;
      u[j][i][2] = 0.0;
      u[j][i][3] = 0.0;
    }
    
  }

  /*
     Restore vectors
  */
  ierr = DMDAVecRestoreArrayDOF(da,U,&u);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PML"
PetscReal PML(PetscReal Lx, PetscReal x, PetscReal width, PetscReal amp)
{
  PetscReal pml = 0.0;
  if(x<width)
  {
    pml += amp*( (width-x)/width - sin(2.0*M_PI*(width-x)/width)/(2.0*M_PI));
  }
  if(Lx-x < width)
  {
    pml += amp*( (x-Lx+width)/width - sin(2.0*M_PI*(x-Lx+width)/width)/(2.0*M_PI));
  }
  return pml;
}
