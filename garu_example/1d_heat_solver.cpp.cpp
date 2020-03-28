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

void initialize(double *var, double *var_new, int nx){
  // initialize phi
  for(int i = 0; i <= nx-1; i++){
    var[i] = 0.; var_new[i] = 0.;
  }
  var[nx/2] = 1.;
}

void simulation(double *var, double *var_new, int nx, double c, double d, double dt){
  double RHS;
  // simulation
  for(int i = 1; i <= nx-2; i++){
    RHS = (c/(d*d))*(var[i+1] - 2.*var[i] + var[i-1]);
    var_new[i] = var[i] + dt*RHS;
  }
}

void visualization(double *var, int nx){
  // visualization
  for(int i = 0; i <= nx-1; i++){
    cout << var[i] << " ";
  }
  cout << "\n";
}

void set_equal(double *var, double *var_new, int nx){
  for(int i = 0; i <= nx-1; i++){
    var[i] = var_new[i];
  }
}

int main(){

  int nx = 11; //int = integer 0 1 2 n
  double k = 1.; // double := real number
  double dx = 1.;
  double dt = 0.01;

  //......................................
  double *phi;   
  phi = (double *) malloc (nx*sizeof(double));
  // ,...............................
  double *ptr;
  double *phi_new;

  phi_new = (double *) malloc (nx*sizeof(double));
  
  initialize(phi, phi_new,nx);
  visualization(phi, nx);
  
  ptr = &phi[0];

  // Begin time loop
  for(int n = 1; n <= 10000; n++){
    simulation(phi, phi_new, nx, k, dx, dt);
    visualization(phi_new, nx);
    set_equal(phi, phi_new, nx);
  }
  

}
