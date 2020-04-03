// Solving 1D heat equation
// dphi/dt = k(d^2 phi/ dx^2)
// k is a coefficient of diffusion (a real number)
// dt = time step; dx = space between 2 nodes              
//!------------------------------------------!
//  bc[0] 1 2  ... bc[11]
//        1 2 ... 10
// boundary conditions phi[0] = 0; phi[nx-1] = 0
//                     phi_new[0] = 0; phi_new[nx-1] = 0
//!------------------------------------------!

//               kdt/(dx*dx) <= a dimensionless number dictating whether our sim will blow 
//                              By experiment, kdt/(dx*dx) <~ 0.5 for numerical stability
// 
//  Physically, we're seeking for numerical solution of diffusing ink (cross-sectionally uniform)
//  through a 1D rod of size L = dx*(nx-1) = 11 meters with coefficient k = 1 m^2/s

#include <stdlib.h> // malloc
#include <stdio.h> // printf
#include <iostream> // cout
#include <cmath> // math lib
#include <fstream> // file I/O
using namespace std;

void initialize(double **var, double **var_new, int nx, int ny){
  // initialize phi
  for(int i = 0; i <= nx-1; i++){
    for(int j = 0; j <= ny-1; j++){
      var[i][j] = 0.; var_new[i][j] = 0.;
    }
  }
  var[nx/2][ny/2] = 5000.;
}

void simulation(double **var, double **var_new, int nx, int ny,  double c, double d, double dt){
  double RHS;
  // simulation
  for(int i = 1; i <= nx-2; i++){
    for(int j = 1; j <= ny-2; j++){
      RHS = (c/(d*d))*(var[i][j-1] + var[i-1][j] - 4.*var[i][j] + var[i+1][j]+var[i][j+1]);
      var_new[i][j] = var[i][j] + dt*RHS;
    }
  }
}

void visualization(double **var, int nx, int ny){
  // visualization
  for(int i = 0; i <= nx-1; i++){
    for(int j = 0; j <= ny-1; j++){
      cout << var[i][j] << " ";
    }
  cout << "\n";
  }
  cout << "\n";
}

void set_equal(double **var, double **var_new, int nx, int ny){
  for(int i = 0; i <= nx-1; i++){
    for(int j = 0; j <= ny-1; j++){
      var[i][j] = var_new[i][j];
    }
  }
}

void write_result(double **var, int nx, int ny, int ts){
  ofstream myfile;
  string file_name = "test2_s" + to_string(ts) + ".vtk";
  myfile.open(file_name);
  /* .......................................................... */
  // Paraview header                                                                                                                                                                     
  myfile << "# vtk DataFile Version 2.0\n";
  myfile << "FlowField\n";
  myfile << "ASCII\n";

  // Grid                                                                                                                                                                                
  myfile << "DATASET STRUCTURED_GRID\n";
  myfile << "DIMENSIONS " << nx << " " << 1 << " " << ny << "\n";
  myfile << "POINTS " << nx*1*ny << " float\n";
  for(int i = 0; i <= nx-1; i++){
    for(int j = 0; j <= ny-1; j++){
      myfile << i << " " << j << " 0\n";
    }
  }

  // Data header                                                                                                                                                                         
  myfile << "\n";
  myfile << "POINT_DATA";
  myfile << " " << nx*ny << "\n";

  // Data                                                                                                                                                                                
  myfile << "\n";
  myfile << "SCALARS U float 1\n";
  myfile << "LOOKUP_TABLE default\n";
  for(int i = 0; i <= nx-1; i++){
    for(int j = 0; j <= ny-1; j++){
      myfile << var[i][j] << "\n";
    }
  }
  /* .......................................................... */
  myfile.close();


}

int main(){

  int nx = 21; //int = integer 0 1 2 n
  int ny = 21;
  double k = 1.; // double := real number
  double dx = 0.05;
  double dt = 0.0005;

  //......................................
  double** phi ;   
  phi = (double**)  malloc (nx*sizeof(double*));
  //phi =  new int* [nx];
  for(int row = 0;row < nx;row++){
    phi[row] = (double*) malloc(ny*sizeof(double));
  }
  double** phi_new ;   
  phi_new = (double**)  malloc (nx*sizeof(double*));
  //phi =  new int* [nx];
  for(int row = 0;row < nx;row++){
    phi_new[row] = (double*) malloc(ny*sizeof(double));
  }
  // ,.....................x..........
  //double *ptr[ny];*/
  //double (*phi_new)[ny];
  /*phi = new int* [nx];
  for(int row = 0;row < nx;row++){
    phi[row] = new int[ny];
    }*/
  //phi_new = (double *)  malloc (ny*nx*sizeof(double));
  
  initialize(phi, phi_new,nx,ny);
  visualization(phi, nx,ny);
  
  // ptr = &phi[0];

  // Begin time loop
  for(int n = 0; n <= 999; n++){
    write_result(phi,nx,ny,n);
    simulation(phi, phi_new, nx,ny, k, dx, dt);
    visualization(phi_new, nx,ny);
    set_equal(phi, phi_new, nx,ny);
    }
  
  

}
