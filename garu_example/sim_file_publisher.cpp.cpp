#include <stdlib.h> // malloc
#include <stdio.h> //printf
#include <iostream> //cout
#include <cmath> //use pow & sqrt
#include <fstream> // file I/O
using namespace std;

int i,j;
const int nx = 51;
const int ny = 51;
int i_c, j_c; 
double tmp;
int phi[nx][ny];

void intialize(){
  // Always Initialize your variables
  for(i = 0; i <= nx-1; i++){
    for(j = 0; j <= ny-1; j++){
      phi[i][j] = 0; 
   }
  }

}
 
void set_phi(){
  // Assign phi = 1 inside a circle around i_c, j_c of "radius" = 15.0
  // letting dx=dy=1 as, e.g., x_spacing is dx*diff(index) 
  //                                     not diff(index) itself
  i_c = (nx-1)/2;
  j_c = (ny-1)/2;
  for(i = 0; i <= nx-1; i++){
    for(j = 0; j <= ny-1; j++){
      if( sqrt( pow(i-i_c,2)  + pow(j-j_c,2) ) < 15.0 ){
	phi[i][j] = 1; } else { phi[i][j] = 0; }
    }
  }


}

void visualize(){
  
  for(j = 0; j <= ny-1; j++){
    for(i = 0; i <= nx-1; i++){
      cout << phi[i][j] << " ";
   }
    cout << "\n";
  }
}

int main(){
  intialize();
  set_phi();
  visualize();

  ofstream myfile;
  myfile.open("test2.vtk");
  /* .......................................................... */
  // Paraview header
  myfile << "# vtk DataFile Version 2.0\n";
  myfile << "FlowField\n";
  myfile << "ASCII\n";

  // Grid
  myfile << "DATASET STRUCTURED_GRID\n";
  myfile << "DIMENSIONS " << nx << " " << 1 << " " << ny << "\n";
  myfile << "POINTS " << nx*1*ny << " float\n";
  for(j = 0; j <= ny-1; j++){
    for(i = 0; i <= nx-1; i++){
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
  for(j = 0; j <= ny-1; j++){
    for(i = 0; i <= nx-1; i++){
      myfile << phi[i][j] << "\n";
    }
  }
  /* .......................................................... */
  myfile.close();

}  
