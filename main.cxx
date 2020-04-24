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
#include <fstream>
#include <stdlib.h>
// Function that makes each processor print in order 
extern void printInOrder(PetscMPIInt rank, PetscMPIInt size, const char* fmt, ...);
extern PetscErrorCode TimeStep(DM da_u_old,DM da_u,DM da_phi1_old,DM da_phi1,DM da_phi2_old,DM da_phi2, \
  PetscReal Lx, PetscReal Ly, Vec vec_u_old, Vec vec_u, Vec vec_phi1_old, Vec vec_phi1, \
  Vec vec_phi2_old, Vec vec_phi2, PetscReal dt, PetscReal pml_width, PetscReal pml_strength, \
  PetscReal t, PetscReal kx, PetscReal ky, PetscInt *m, PetscReal *f, PetscReal* d, PetscReal* alpha1\
  , PetscReal *alpha2);
extern PetscErrorCode FormInitialSolution( DM da_u_old, DM da_u, DM da_phi1_old, DM da_phi1, DM da_phi2_old, DM da_phi2\
  , Vec vec_u_old, Vec vec_u, Vec vec_phi1_old, Vec vec_phi1, Vec vec_phi2_old, Vec vec_phi2 , \
  PetscReal Lx, PetscReal Ly, PetscReal dt,PetscReal width);
extern PetscErrorCode MyMonitor(TS,PetscInt,PetscReal,Vec,void*);
extern PetscReal PML(PetscReal Lx, PetscReal x, PetscReal width, PetscReal amp);
using namespace std;
// Main
int main(int argc,char **argv)
{
  //---------------------------------------------------------------------------
  // Number of global gridpoints in the x and y direction
  PetscInt       Mx = 10, My = 10, Nt = 1;
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
  PetscInt       *m;
  PetscReal      *f, *d, *alpha1, *alpha2;
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
  ierr = PetscOptionsGetInt(NULL,NULL,"-M",&Mx,NULL);CHKERRQ(ierr);
  My = Mx;
  ierr = PetscOptionsGetInt(NULL,NULL,"-Nt",&Nt,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-dt",&dt,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-L",&Lx,NULL);CHKERRQ(ierr);
  Ly = Lx;
  ierr = PetscOptionsGetReal(NULL,NULL,"-pml_width",&pml_width,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-pml_strength",&pml_strength,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-kx",&kx,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-ky",&ky,NULL);CHKERRQ(ierr);

  m = new PetscInt[Mx*My];
  d = new PetscReal[Mx*My];
  f = new PetscReal[Mx*My];
  alpha1 = new PetscReal[Mx*My];
  alpha2 = new PetscReal[Mx*My];


  ifstream myReadFile;
  myReadFile.open("mij.out");
  PetscReal output;
  string str;
  
  for(int j = 0; j < My; j++)
  {
    for(int i = 0; i <Mx ; i++)
    {
      
      getline(myReadFile,str);
      output = atof(str.c_str());
      m[j*Mx + i] = (PetscInt) round(output);
      std::cout <<  i << "," << j <<std::endl;
      
    }
  }

  for(int j = 0; j < My; j++)
  {
    for(int i = 0; i <Mx ; i++)
    {
      std::cout <<  m[j*Mx + i]<<std::endl;
    }
  }
  
  myReadFile.close();

  myReadFile.open("alpha_1.out");
  for(int j = 0; j < My; j ++)
  {
    for(int i = 0; i <Mx ; i ++)
    {
      getline(myReadFile,str);
      output = atof(str.c_str());
      alpha1[j*Mx + i] = output*(Mx-1)/Lx;
    }
  }
  myReadFile.close();

  myReadFile.open("alpha_2.out");
  for(int j = 0; j < My; j ++)
  {
    for(int i = 0; i <Mx ; i ++)
    {
     
      getline(myReadFile,str);
      output = atof(str.c_str());
      alpha2[j*Mx + i] = output*(My-1)/Ly;
      
    }
  }
  myReadFile.close();

  for(int j =0; j<My;j++)
  {
    for(int i=0;i<Mx;i++)
    {
      d[j*Mx+i] = 0.0;
      if (not isnan((1.0-alpha1[j*Mx+i])) and not isinf((1.0-alpha1[j*Mx+i])/alpha1[j*Mx+i]))
        d[j*Mx + i] += (1.0-alpha1[j*Mx+i])/alpha1[j*Mx+i];
      if (not isnan((1.0-alpha2[j*Mx+i])/alpha2[j*Mx+i]) and not isinf((1.0-alpha2[j*Mx+i])/alpha2[j*Mx+i]))
        d[j*Mx + i] += (1.0-alpha2[j*Mx+i])/alpha2[j*Mx+i];
      
      
      if(m[j*Mx+i] == -1)
        std::cout << "one = "<< d[j*Mx + i] << (1.0-alpha1[j*Mx+i])/alpha1[j*Mx+i]<< " two = " << (1.0-alpha2[j*Mx+i])/alpha2[j*Mx+i] << std::endl;
    }
  }

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
  DMDASetFieldName(da_u,0,"u");

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
  ierr = VecNorm(vec_u_old,NORM_2,&norm);CHKERRQ(ierr);
  if (rank == 0)
    std::cout << "Initial Norm = " << norm << std::endl;

  for(int i = 0; i< Nt; i++)
  {
    // Timestep 
    TimeStep(da_u_old, da_u, da_phi1_old, da_phi1, da_phi2_old, da_phi2\
      ,Lx,Ly,vec_u_old, vec_u, vec_phi1_old, vec_phi1, vec_phi2_old, \
      vec_phi2,dt,pml_width,pml_strength,dt*i,kx,ky,m,f,d,alpha1,alpha2);
      

    ierr = VecNorm(vec_u_old,NORM_2,&norm);CHKERRQ(ierr);
    if (rank == 0)
      std::cout << "Norm = " << std::setprecision (16) << norm << std::endl;
    std::string s = "wave"+std::to_string(i)+".out";
    PetscViewerASCIIOpen(PETSC_COMM_WORLD,s.c_str(),&viewer);
    VecView(vec_u_old,viewer);

  }


  // Free space 
  delete [] m;
  delete [] d;
  delete [] f;
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
  PetscViewerDestroy(&viewer);
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
  PetscReal t, PetscReal kx, PetscReal ky, PetscInt *m, PetscReal *f, PetscReal* d, PetscReal* alpha1\
  , PetscReal *alpha2)
{
  PetscErrorCode ierr;
  PetscInt       i,j,Mx,My,xs,ys,xm,ym;
  PetscReal      hx,hy,/*hxdhy,hydhx,*/ sx,sy;
  PetscScalar    phi1_old,phi1,phi2_old,phi2,u_old,uxx,uyy,u,***array_u_old,***array_u,***array_phi1_old, ***array_phi1;
  PetscScalar    ***array_phi2_old, ***array_phi2, tmp_x;
  Vec            local_u_old, local_u, local_phi1_old, local_phi1, local_phi2_old, local_phi2;

  PetscFunctionBeginUser;
  ierr = DMGetLocalVector(da_u_old,&local_u_old);CHKERRQ(ierr);
  ierr = DMDAGetInfo(da_u_old,PETSC_IGNORE,&Mx,&My,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,
                     PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE);

  ierr = DMGetLocalVector(da_u,&local_u);CHKERRQ(ierr);

  ierr = DMGetLocalVector(da_phi1_old,&local_phi1_old);CHKERRQ(ierr);

  ierr = DMGetLocalVector(da_phi1,&local_phi1);CHKERRQ(ierr);
  
  ierr = DMGetLocalVector(da_phi2_old,&local_phi2_old);CHKERRQ(ierr);
  
  ierr = DMGetLocalVector(da_phi2,&local_phi2);CHKERRQ(ierr);
 

  PetscReal dt2 = dt*dt;
  hx = Lx/(PetscReal)(Mx-1); sx = 1.0/(hx*hx);
  hy = Ly/(PetscReal)(My-1); sy = 1.0/(hy*hy);

  ierr = DMGlobalToLocalBegin(da_u_old,vec_u_old,INSERT_VALUES,local_u_old);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(da_u_old,vec_u_old,INSERT_VALUES,local_u_old);CHKERRQ(ierr);
  ierr = DMDAVecGetArrayDOF(da_u_old,local_u_old,&array_u_old);CHKERRQ(ierr);

  ierr = DMGlobalToLocalBegin(da_u,vec_u,INSERT_VALUES,local_u);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(da_u,vec_u,INSERT_VALUES,local_u);CHKERRQ(ierr);
  ierr = DMDAVecGetArrayDOF(da_u,local_u,&array_u);CHKERRQ(ierr);

  ierr = DMGlobalToLocalBegin(da_phi1_old,vec_phi1_old,INSERT_VALUES,local_phi1_old);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(da_phi1_old,vec_phi1_old,INSERT_VALUES,local_phi1_old);CHKERRQ(ierr);
  ierr = DMDAVecGetArrayDOF(da_phi1_old,local_phi1_old,&array_phi1_old);CHKERRQ(ierr);

  ierr = DMGlobalToLocalBegin(da_phi1,vec_phi1,INSERT_VALUES,local_phi1);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(da_phi1,vec_phi1,INSERT_VALUES,local_phi1);CHKERRQ(ierr);
  ierr = DMDAVecGetArrayDOF(da_phi1,local_phi1,&array_phi1);CHKERRQ(ierr);

  ierr = DMGlobalToLocalBegin(da_phi2_old,vec_phi2_old,INSERT_VALUES,local_phi2_old);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(da_phi2_old,vec_phi2_old,INSERT_VALUES,local_phi2_old);CHKERRQ(ierr);
  ierr = DMDAVecGetArrayDOF(da_phi2_old,local_phi2_old,&array_phi2_old);CHKERRQ(ierr);

  ierr = DMGlobalToLocalBegin(da_phi2,vec_phi2,INSERT_VALUES,local_phi2);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(da_phi2,vec_phi2,INSERT_VALUES,local_phi2);CHKERRQ(ierr);
  ierr = DMDAVecGetArrayDOF(da_phi2,local_phi2,&array_phi2);CHKERRQ(ierr);

  ierr = DMDAGetCorners(da_u,&xs,&ys,NULL,&xm,&ym,NULL);CHKERRQ(ierr);

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
      u_old      = array_u_old[j][i][0];
      u          = array_u[j][i][0];
      
      if (m[j*Mx+i]==0)
      {
        PMLX = PML(Lx,i*hx,pml_width,pml_strength );
        uxx        = (-2.0*u + array_u[j][i-1][0] + array_u[j][i+1][0])*sx;
        uyy        = (-2.0*u + array_u[j-1][i][0] + array_u[j+1][i][0])*sy;
        if (i*hx<pml_width||Lx-i*hx<pml_width||j*hy<pml_width||Ly-j*hy<pml_width)
        {
          phi1_old = array_phi1_old[j][i][0];
          phi1 = array_phi1[j][i][0];
          phi2_old = array_phi2_old[j][i][0];
          phi2 = array_phi2[j][i][0];

          array_u_old[j][i][0] = (dt2*(uxx + uyy) +(2.0-dt2*PMLX*PMLY)*u+(dt*(PMLX+PMLY)-1.0)*u_old+ dt2*(array_phi1[j][i+1][0]\
            -array_phi1[j][i+1][0])/(2*hx) + dt2*(array_phi2[j+1][i][0]-array_phi2[j-1][i][0])/(2*hy))/(1.0+dt*(PMLX+PMLY));
          array_phi1_old[j][i][0] = 2.0*dt*(-PMLX*phi1 +(PMLY-PMLX)*(array_u[j][i+1][0] - array_u[j][i-1][0])/(2*hx)) + phi1_old;
          array_phi2_old[j][i][0] = 2.0*dt*(-PMLY*phi2 +(PMLX-PMLY)*(array_u[j+1][i][0] - array_u[j-1][i][0])/(2*hy)) + phi2_old;
        }
        else
        {
          array_u_old[j][i][0] = dt2*(uxx + uyy) + 2.0*u -1.0*u_old;
        }
      }

      // Set all values inside of sphere = 0 (TEMP UNTIL EMBEDDED BOUNDARY SET)
      if(m[j*Mx+i] == -1)
      {
        f[j*Mx+i] = 0.0;
        if (not isnan(((1-alpha1[j*Mx+i])*array_u[j][i][0]+alpha1[j*Mx+i]*array_u[j][i-1][0])/alpha1[j*Mx+i]) and not isinf(((1-alpha1[j*Mx+i])*array_u[j][i][0]+alpha1[j*Mx+i]*array_u[j][i-1][0])/alpha1[j*Mx+i]))
        {
          f[j*Mx + i] += ((1-alpha1[j*Mx+i])*array_u[j][i][0]+alpha1[j*Mx+i]*array_u[j][i-1][0])/alpha1[j*Mx+i];
        }
        else
        {
          f[j*Mx + i] = array_u[j][i];
        }
        if (not isnan(((1-alpha2[j*Mx+i])*array_u[j][i][0]+alpha2[j*Mx+i]*array_u[j+1][i][0])/alpha2[j*Mx+i]) and not isinf(((1-alpha2[j*Mx+i])*array_u[j][i][0]+alpha2[j*Mx+i]*array_u[j+1][i][0])/alpha2[j*Mx+i]))
          f[j*Mx + i] += ((1-alpha2[j*Mx+i])*array_u[j][i][0]+alpha2[j*Mx+i]*array_u[j+1][i][0])/alpha2[j*Mx+i];
  

        array_u_old[j][i][0] = (dt2*sy*(-4.0*u+array_u[j][i+1][0]+array_u[j-1][i][0])\
          - dt2*sx*0.5*u_old*d[j*Mx+i] + dt2*sx*f[j*Mx+i])/(1.0+dt2*sx*0.5*d[j*Mx+i]);
      }

      if (m[j*Mx+i] == 1 )
      {
        array_u_old[j][i][0] = -sin(kx*(i*hx)+ky*(hy*j)-sqrt(kx*kx+ky*ky)*t);
      }

      /*if ( (i*hx-Lx/2.0)*(i*hx-Lx/2.0) < 0.1*0.1 and (j*hy-Ly/2.0)*(j*hy-Ly/2.0) < 0.1*0.1)
      {
        array_u_old[j][i][0] = -sin(kx*i*hx+ky*hy*j-sqrt(kx*kx+ky*ky)*t);
      }*/
    }
  }

  for (j=ys; j<ys+ym; j++) {
    for (i=xs; i<xs+xm; i++) {
      tmp_x = array_u_old[j][i][0];
      array_u_old[j][i][0] = array_u[j][i][0];
      array_u[j][i][0] = tmp_x;
      if (i*hx<pml_width||Lx-i*hx<pml_width||j*hy<pml_width||Ly-j*hy<pml_width)
      {
        tmp_x = array_phi1_old[j][i][0];
        array_phi1_old[j][i][0] = array_phi1[j][i][0];
        array_phi1[j][i][0] = tmp_x;
        tmp_x = array_phi2_old[j][i][0];
        array_phi2_old[j][i][0] = array_phi2[j][i][0];
        array_phi2[j][i][0] = tmp_x;
      }
    }
  }


  /*
     Restore vectors
  */
  ierr = DMDAVecRestoreArrayDOF(da_u_old,local_u_old,&array_u_old);CHKERRQ(ierr);
  DMLocalToGlobalBegin(da_u_old,local_u_old,INSERT_VALUES,vec_u_old);
  DMLocalToGlobalEnd(da_u_old,local_u_old,INSERT_VALUES,vec_u_old);
  ierr = DMRestoreLocalVector(da_u_old,&local_u_old);CHKERRQ(ierr);

  ierr = DMDAVecRestoreArrayDOF(da_u,local_u,&array_u);CHKERRQ(ierr);
  DMLocalToGlobalBegin(da_u,local_u,INSERT_VALUES,vec_u);
  DMLocalToGlobalEnd(da_u,local_u,INSERT_VALUES,vec_u);
  ierr = DMRestoreLocalVector(da_u,&local_u);CHKERRQ(ierr);

  ierr = DMDAVecRestoreArrayDOF(da_phi1_old,local_phi1_old,&array_phi1_old);CHKERRQ(ierr);
  DMLocalToGlobalBegin(da_phi1_old,local_phi1_old,INSERT_VALUES,vec_phi1_old);
  DMLocalToGlobalEnd(da_phi1_old,local_phi1_old,INSERT_VALUES,vec_phi1_old);
  ierr = DMRestoreLocalVector(da_phi1_old,&local_phi1_old);CHKERRQ(ierr);

  ierr = DMDAVecRestoreArrayDOF(da_phi1,local_phi1,&array_phi1);CHKERRQ(ierr);
  DMLocalToGlobalBegin(da_phi1,local_phi1,INSERT_VALUES,vec_phi1);
  DMLocalToGlobalEnd(da_phi1,local_phi1,INSERT_VALUES,vec_phi1);
  ierr = DMRestoreLocalVector(da_phi1,&local_phi1);CHKERRQ(ierr);

  ierr = DMDAVecRestoreArrayDOF(da_phi2_old,local_phi2_old,&array_phi2_old);CHKERRQ(ierr);
  DMLocalToGlobalBegin(da_phi2_old,local_phi2_old,INSERT_VALUES,vec_phi2_old);
  DMLocalToGlobalEnd(da_phi2_old,local_phi2_old,INSERT_VALUES,vec_phi2_old);
  ierr = DMRestoreLocalVector(da_phi2_old,&local_phi2_old);CHKERRQ(ierr);

  ierr = DMDAVecRestoreArrayDOF(da_phi2,local_phi2,&array_phi2);CHKERRQ(ierr);
  DMLocalToGlobalBegin(da_phi2,local_phi2,INSERT_VALUES,vec_phi2);
  DMLocalToGlobalEnd(da_phi2,local_phi2,INSERT_VALUES,vec_phi2);
  ierr = DMRestoreLocalVector(da_phi2,&local_phi2);CHKERRQ(ierr);
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
  PetscScalar    **array_u_old,**array_u,**array_phi1_old, **array_phi1;
  PetscScalar    **array_phi2_old, **array_phi2;
  
  PetscFunctionBeginUser;
  ierr = DMDAGetInfo(da_u_old,PETSC_IGNORE,&Mx,&My,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,
                     PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE);

  /*
     Get pointers to vector data
  */
  ierr = DMDAVecGetArrayDOF(da_u_old,vec_u_old,&array_u_old);CHKERRQ(ierr);
  ierr = DMDAVecGetArrayDOF(da_u,vec_u,&array_u);CHKERRQ(ierr);
  ierr = DMDAVecGetArrayDOF(da_phi1_old,vec_phi1_old,&array_phi1_old);CHKERRQ(ierr);
  ierr = DMDAVecGetArrayDOF(da_phi1,vec_phi1,&array_phi1);CHKERRQ(ierr);
  ierr = DMDAVecGetArrayDOF(da_phi2_old,vec_phi2_old,&array_phi2_old);CHKERRQ(ierr);
  ierr = DMDAVecGetArrayDOF(da_phi2,vec_phi2,&array_phi2);CHKERRQ(ierr);

  /*
     Get local grid boundaries
  */
  ierr = DMDAGetCorners(da_u_old,&xs,&ys,NULL,&xm,&ym,NULL);CHKERRQ(ierr);

  /*
     Compute function over the locally owned part of the grid
  */
  for (j=ys; j<ys+ym; j++) 
  {
    for (i=xs; i<xs+xm; i++) 
    {
      array_u_old[j][i] = 0.0;
      array_u[j][i] = 0.0;
      array_phi1_old[j][i] = 0.0;
      array_phi1[j][i] = 0.0;
      array_phi2_old[j][i] = 0.0;
      array_phi2[j][i] = 0.0;
    }
    
  }

  /*
     Restore vectors
  */
  ierr = DMDAVecRestoreArrayDOF(da_u_old,vec_u_old,&array_u_old);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArrayDOF(da_u,vec_u,&array_u);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArrayDOF(da_phi1_old,vec_phi1_old,&array_phi1_old);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArrayDOF(da_phi1,vec_phi1,&array_phi1);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArrayDOF(da_phi2_old,vec_phi2_old,&array_phi2_old);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArrayDOF(da_phi2,vec_phi2,&array_phi2);CHKERRQ(ierr);
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
