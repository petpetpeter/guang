#include <stdlib.h> // malloc
#include <stdio.h> // printf
#include <iostream> // cout
#include <cmath> // math lib
#include <fstream> // file I/O
using namespace std;


int find_max(double mat[51][51],int nx,int ny) 
{ // Initializing max element as INT_MIN 
    double maxElement = 0.;
    for (int i = 0; i < nx; i++) { 
        for (int j = 0; j < ny; j++) { 
            if (abs(mat[i][j]) > maxElement) { 
                maxElement = abs(mat[i][j]); 
            } 
        } 
    } 
  
    // finally return maxElement 
    return maxElement; 
} 

void initialize(double **var, double **var_new,double **pressure,double **pressure_new, int nx, int ny){
  // initialize phi
  for(int i = 0; i <= nx-1; i++){
    for(int j = 0; j <= ny-1; j++){
      var[i][j] = 0.; var_new[i][j] = 0.;
    }
  }
  var[nx/2][ny/2] = 10.;
  // initialize pressure
}

int poisson(double **pressure,double **pressure_new,double **RHS,double **residule,int nx,int ny,double delx,double dely,double eps,int itermax,double omg){
  double res[nx][ny];
  double r_max;
  for(int i = 0; i <= nx-1; i++){
    for(int j = 0; j <= ny-1; j++){
      res[i][j] = 0.;
    }
  }
  for(int it = 1; it <= itermax;it++){
    for(int i = 1; i <= nx-2; i++){
      for(int j = 1; j <= ny-2; j++){
	      pressure_new[i][j] = (1-omg)*pressure[i][j] + (omg/((2/(delx*delx))+(2/(dely*dely)))) *
                           ( ((pressure[i+1][j]+pressure_new[i-1][j])/(delx*delx)) 
                           + ((pressure[i+1][j+1]+pressure_new[i][j-1])/(dely*dely))
                           - RHS[i][j] );
        }
    }
  set_equal(pressure,pressure_new, nx,ny);
  for(int i = 1; i <= nx-2; i++){
    for(int j = 1; j <= ny-2; j++){
      res[i][j] = ((pressure[i+1][j]-pressure[i][j]) - (pressure[i][j]-pressure[i-1][j]))/(delx*delx)
                  +  ((pressure[i][j+1]-pressure[i][j]) - (pressure[i][j]-pressure[i][j-1]))/(dely*dely)
                  - RHS[i][j];
    }
  }
  r_max = find_max(res,nx,ny);
  if (r_max > eps){
      cout << "converged iteration";
      break;
  }
  }
  cout << "full iteration";
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
  myfile << " " << nx*1 << "\n";

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

  int nx = 51; //int = integer 0 1 2 n
  int ny = 51;
  double k = 1.; // double := real number
  double dx = 1.;
  double dt = 0.1;

  //......................................
  double** phi ;   
  phi = (double**)  malloc (nx*sizeof(double*));
  for(int row = 0;row < nx;row++){
    phi[row] = (double*) malloc(ny*sizeof(double));
  }
  double** phi_new ;   
  phi_new = (double**)  malloc (nx*sizeof(double*));
  for(int row = 0;row < nx;row++){
    phi_new[row] = (double*) malloc(ny*sizeof(double));
  }
  double** pressure ;   
  phi = (double**)  malloc (nx*sizeof(double*));
  for(int row = 0;row < nx;row++){
    phi[row] = (double*) malloc(ny*sizeof(double));
  }
  double** pressure_new ;   
  phi = (double**)  malloc (nx*sizeof(double*));
  for(int row = 0;row < nx;row++){
    phi[row] = (double*) malloc(ny*sizeof(double));
  }
  double** residule ;   
  phi = (double**)  malloc (nx*sizeof(double*));
  for(int row = 0;row < nx;row++){
    phi[row] = (double*) malloc(ny*sizeof(double));
  }
  
  initialize(phi, phi_new,pressure,pressure_new,nx,ny);
  visualization(phi, nx,ny);
  
  // ptr = &phi[0];

  // Begin time loop
  for(int n = 1; n <= 10; n++){
    simulation(phi, phi_new, nx,ny, k, dx, dt);
    visualization(phi_new, nx,ny);
    write_result(phi_new,nx,ny,n);
    set_equal(phi, phi_new, nx,ny);
    }
  
  

}
