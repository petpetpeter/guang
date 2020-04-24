#include <stdlib.h> // malloc
#include <stdio.h> // printf
#include <iostream> // cout
#include <cmath> // math lib
#include <fstream> // file I/O
#include <iomanip> //forsomething
using namespace std;

void visualization_gen(double** gen, int nx, int ny){
  cout << setprecision(4) << fixed;
  for ( int i = 0; i<=ny-1; i++){
    for ( int j = 0; j<=nx-1; j++){
      cout << gen[i][j] << " ";
    }
    cout << "\n";
  }
  cout << "\n\n";
}


void initualization(double**u, double**u_new, double**v, double**v_new, double**p, double**p_new, double**F, double**G, double**p_rhs, int nx, int ny){
  for (int i = 0; i <= ny-1; i++){
    for (int j = 0; j <= nx-1; j++){
      u[i][j] = 0;
      u_new[i][j] = 0;
      v[i][j] = 0;
      v_new[i][j] = 0;
      p[i][j] = 0;
      p_new[i][j] = 0;
      F[i][j] = 0;
      G[i][j] =0;
      p_rhs[i][j] =0;
    }
  }
}


void setbound(double**u, double**u_new, double**v, double**v_new, double** p,double** p_new, int nx, int ny, double u_in){
  for (int i = 1; i <= ny-2; i++){
    u_new[i][0] = u_in;
    v_new[i][0] = -v_new[i][1];// inflow for west bound
    u_new[i][nx-1] = u_new[i][nx-2];
    v_new[i][nx-1] = v_new[i][nx-2];// outflow for east bound
  }
  
  for (int j = 0; j <= nx-1;j++){
    u_new[0][j] = -u_new[1][j]; // no slip north
    v_new[0][j] =  0;// no slip condition for north bound
    u_new[ny-1][j]= -u_new[ny-2][j];// no slup south u_new[ny-2][j];
    v_new[ny-1][j]=  0;// free slip conition for south bound for symmetry
  }

}


void alittlestepforhumankind ( double**u, double**u_new, double**v, double**v_new, int nx, int ny){
  for (int i = 0; i<=ny-1; i++){
    for (int j = 0; j<=nx-1;j++){
      u[i][j] = u_new[i][j];
      v[i][j] = v_new[i][j];
    }
  }
}


void showinipls ( double**u, double**v, double**p, int nx, int ny){
  cout << "t=0 \n";
  cout << "u" << "\n";
  for ( int i = 0; i<=ny-1; i++){
    for ( int j = 0; j<=nx-1; j++){
      cout << u[i][j] << " ";
    }
    cout << "\n";
  }
  cout << "\n\n";

  cout << "v" << "\n";
  for ( int i = 0; i<=ny-1; i++){
    for ( int j = 0; j<=nx-1; j++){
      cout << v[i][j] << " ";
    }
    cout << "\n";
  }
  cout << "\n\n";

  cout << "p" << "\n";
  for ( int i = 0; i<=ny-1; i++){
    for ( int j = 0; j<=nx-1; j++){
      cout << p[i][j] << " ";
    }
    cout << "\n";
  }
  cout << "\n\n\n";
}


void comp_FG(double**u, double**v, double**F, double**G, double dx, double dy, double donor, double Re, double gx, double gy, double ts, int nx, int ny){
  //upkob FG
  for (int i = 0; i <= ny-1;i++){
    F[i][0] = u[i][0]; //Fkobeast
    F[i][nx-1] = u[i][nx-1]; // Fkobwest
  }
  for(int j =0;j<=nx-1;j++){
    G[0][j] = v[0][j];
    G[nx-1][j] = v[nx-1][j];
  }
  for ( int i = ny-2; i >= 1; i--){
    for ( int j = 1; j <= nx-2; j++){
      
      double duu_dx = ( ( pow(((u[i][j]+u[i][j+1])/2),2) - pow(((u[i][j-1]+u[i][j])/2),2) ) / dx )  + (donor/dx)*( abs(u[i][j]+u[i][j+1])*(u[i][j]-u[i][j+1])/4 - abs(u[i][j-1]+u[i][j])*(u[i][j-1]-u[i][j])/4 ) ;

      double duv_dy = ( ( (v[i][j]+v[i][j+1])*(u[i][j]+u[i-1][j])/4 - (v[i+1][j]+v[i+1][j+1])*(u[i+1][j]+u[i][j])/4 ) / dy )  +  (donor/dy)*( abs(v[i][j]+v[i][j+1])*(u[i][j]-u[i-1][j])/4 - abs(v[i+1][j] + v[i+1][j+1])*(u[i+1][j]-u[i][j])/4 ) ;

      double ddu_dxdx = (u[i][j+1]-2*u[i][j]+u[i][j-1])/ pow(dx,2);

      double ddu_dydy = (u[i-1][j]-2*u[i][j]+u[i+1][j])/pow( dy,2);

      F[i][j] = u[i][j] + ts*((ddu_dxdx+ddu_dydy)/Re - duu_dx - duv_dy + gx); //Fchecked


      double dvv_dy = ( (pow( (v[i][j]+v[i-1][j])/2,2) - pow((v[i+1][j]+v[i][j])/2,2) ) / dy)   +  (donor/dy)*( abs(v[i][j]+v[i-1][j])*(v[i][j]-v[i-1][j])/4 - abs(v[i+1][j]+v[i][j])*(v[i+1][j]-v[i][j])/4 ) ;

      double duv_dx = ( ( (u[i][j]+u[i-1][j])*(v[i][j]+v[i][j+1])/4 - (u[i][j-1]+u[i-1][j-1])*(v[i][j-1]+v[i][j])/4 ) / dx )  +  (donor/dx)*( abs(u[i][j]+u[i-1][j])*(v[i][j]-v[i][j+1])/4 - abs(u[i][j-1] + u[i-1][j-1])*(v[i][j-1]-v[i][j])/4 ) ;

      double ddv_dxdx = (v[i][j+1]-2*v[i][j]+v[i][j-1])/ pow(dx,2);

      double ddv_dydy = (v[i-1][j]-2*v[i][j]+v[i+1][j])/ pow(dy,2);

      G[i][j] = v[i][j] + ts*((ddv_dxdx+ddv_dydy)/Re - duv_dx - dvv_dy + gy);
    }
  }
}


void comp_p_rhs( double**F, double**G, double**p_rhs,  double ts, double dx, double dy, int nx, int ny){
  //cout << "rhs" << "\n";
  for ( int i = ny-2; i>=1; i--){
    for ( int j = 1; j<=nx-2; j++){
      p_rhs[i][j] = ( (F[i][j]-F[i][j-1])/dx + (G[i][j]-G[i+1][j])/dy )/ts;
    }
  }

  //visualization_gen(p_rhs,nx,ny);
  
}

double find_max(double** mat,int nx,int ny)
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


void comp_p( double**p, double**p_new, double**p_rhs, double dx, double dy,  int nx, int ny, int it_max, double eps){
  int it = 1;
  double r = 9999;
  double** res ;
  res = (double**)  malloc (nx*sizeof(double*));
  for(int row = 0;row < nx;row++){
    res[row] = (double*) malloc(ny*sizeof(double));
  }

  while ( it <= it_max && abs(r) >= eps){
    //get p_new boundary
    for ( int i=ny-1; i>=0; i--){
      p_new[i][0] = p[i][1];
      p_new[i][nx-1] = p[i][nx-2]; // presure strip for west and east bound
    }
    for ( int j=0; j<=nx-1; j++){
      p_new[0][j] = p[1][j];
      p_new[ny-1][j] = p[ny-2][j]; // pressure strip for north and south bound 
    }
    
    double sum_rsqr = 0;
    //cout << "pressure_lang_up_kob" << "\n";
    //visualization_gen(p_new,nx,ny);
    for ( int i = ny-2; i>=1; i--){
      for ( int j = 1; j<=nx-2; j++){
	p_new[i][j] = (-0.7)*p[i][j] + (1.7/( 2/pow(dx,2) + 2/pow(dy,2) ))*( (p[i][j+1]+p_new[i][j-1])/pow(dx,2) + (p[i-1][j]+p_new[i+1][j])/pow(dy,2) - p_rhs[i][j]);
      }
    }
    //cout << "pressure_lang_up_value" << "\n";
    //visualization_gen(p_new,nx,ny);
    
    /*
    for ( int i = 1; i<=ny-2; i++){
      for ( int j = 1; j<=nx-2; j++){
	double rsqr = pow(((p_new[i][j+1] - 2*p_new[i][j] + p_new[i][j-1])/pow(dx,2) + (p_new[i-1][j] - 2*p_new[i][j] + p_new[i+1][j])/pow(dy,2) - p_rhs[i][j]),2);
	sum_rsqr = sum_rsqr + rsqr;
	cout << rsqr << " "; 
      }
      cout << "\n";
    }

    r = pow((sum_rsqr/ ((nx-2)*(ny-2))),0.5);
    */ 
    //cout << "res" << "\n";
    for(int i = 1; i <= ny-2; i++){
      for(int j = 1; j <= nx-2; j++){
	res[i][j] = ((p_new[i][j+1]-p_new[i][j]) - (p_new[i][j]-p_new[i][j-1]))/(dx*dx)
                  +  ((p_new[i+1][j]-p_new[i][j]) - (p_new[i][j]-p_new[i-1][j]))/(dy*dy)
                  - p_rhs[i][j];
	 }
    }

    //visualization_gen(res,nx,ny);
    
    //cout << "\n";
    r = find_max(res,nx,ny);
    
    it++;
    //cout << "residual = " << r << "\n" << "interation =" << it << "\n";
    //update p
    for ( int i=0; i<=ny-1; i++){
      for ( int j=0; j<=nx-1; j++){
	p[i][j] = p_new[i][j];
      }
    }
  }
  cout << "residual = " << abs(r) << "\n" << "interation =" << it << "\n";
  
}
  

void finalcomp( double**u_new, double**v_new, double**F, double**G, double**p_new, double ts, double dx, double dy, int nx, int ny){
  for (int i=ny-2; i>=0; i--){
    for ( int j =1; j<=nx-2; j++){
      u_new[i][j] = F[i][j] - (ts/dx)*(p_new[i][j+1]-p_new[i][j]);
      v_new[i][j] = G[i][j] - (ts/dy)*(p_new[i-1][j]-p_new[i][j]);
    }
  }
}



void visualization_u(double** u_new, int nx, int ny, int t){
  cout << setprecision(4) << fixed;
  cout << "t= " << t<< "\n";
  cout << "u_new" << "\n";
  for ( int i = 0; i<=ny-1; i++){
    for ( int j = 0; j<=nx-1; j++){
      cout << u_new[i][j] << " ";
    }
    cout << "\n";
  }
  cout << "\n\n";
}

void visualization_v(double **v_new,  int nx, int ny, int t){
  cout << setprecision(4) << fixed;
  cout << "t= " << t<< "\n";
  cout << "v_new" << "\n";
  for ( int i = 0; i<=ny-1; i++){
    for ( int j = 0; j<=nx-1; j++){
      cout << v_new[i][j] << " ";
    }
    cout << "\n";
  }
  cout << "\n\n";

}

void visualization_p(double **p_new, int nx, int ny, int t){
  cout << setprecision(5) << fixed;
  cout << "t= " << t<< "\n";
  cout << "p_new" << "\n";
  for ( int i = 0; i<=ny-1; i++){
    for ( int j = 0; j<=nx-1; j++){
      cout << p_new[i][j] << " ";
    }
    cout << "\n";
  }
  cout << "\n\n";
}
 
       
int main(){

  int nx = 11;
  int ny = 7;
  double dx = 0.5;
  double dy = 0.1;
  double donor = 0;
  double gx = 0;
  double gy = 0;
  int t_max = 1000;
  int it_max = 1000; //for poisson iteration
  double eps = 0.00001; //for poisson iteration


  
  double u_in = 0.1;
  double Re = (1*u_in*ny*dy)/0.1;//1500;
  
  //double ts = Re/((1/pow(dx,2))+(1/pow(dy,2)))/2*1.1; // choose ts based on 3.5 with a sf of 1.1
  double ts = 0.0001;
  ////////////////////////////////////////////////////CREATE ARRAYS///////////////////////////////////////////
  double*u[ny];
  for ( int i = 0; i <= ny-1; i++){
    u[i] = (double*)malloc(nx*sizeof(double));
  }
  double*ptr;
  double*u_new[ny];
  for ( int i = 0; i<= ny-1; i++){
    u_new[i] = (double*)malloc(nx*sizeof(double));
  }
  double*v[ny];
  for ( int i = 0; i <= ny-1; i++){
    v[i] = (double*)malloc(nx*sizeof(double));
  }
  double*v_new[ny];
  for ( int i = 0; i<= ny-1; i++){
    v_new[i] = (double*)malloc(nx*sizeof(double));
  }
  double*p[ny];
  for ( int i = 0; i <= ny-1; i++){
    p[i] = (double*)malloc(nx*sizeof(double));
  }
  double*p_new[ny];
  for ( int i = 0; i<= ny-1; i++){
    p_new[i] = (double*)malloc(nx*sizeof(double));
  }
  double*F[ny];
  for ( int i = 0; i <= ny-1; i++){
    F[i] = (double*)malloc(nx*sizeof(double));
  }
  double*G[ny];
  for ( int i = 0; i<= ny-1; i++){
    G[i] = (double*)malloc(nx*sizeof(double));
  }
  double*p_rhs[ny];
  for ( int i = 0; i <= ny-1; i++){
    p_rhs[i] = (double*)malloc(nx*sizeof(double));
  }
  ////////////////////////////////////////////////////CREATE ARRAYS///////////////////////////////////////////////////////

  initualization(u, u_new, v, v_new, p, p_new, F, G, p_rhs, nx, ny);

  setbound(u, u_new, v, v_new,p, p_new, nx, ny, u_in);

  alittlestepforhumankind (u, u_new, v, v_new, nx, ny);

  showinipls ( u, v, p, nx, ny);

  ptr=&u[0][0];



  for ( int t = 0; t <= t_max-1; t++){
    
    comp_FG(u, v, F, G, dx, dy, donor, Re, gx, gy,ts, nx, ny);

    comp_p_rhs( F, G, p_rhs, ts, dx, dy, nx, ny);

    comp_p( p, p_new, p_rhs, dx, dy, nx, ny, it_max, eps);

    finalcomp( u_new, v_new, F, G, p_new, ts, dx, dy, nx, ny);

    setbound( u, u_new, v, v_new,p,p_new, nx, ny, u_in);

    alittlestepforhumankind ( u, u_new, v, v_new, nx, ny);

    visualization_u(u_new, nx, ny,t);
    //visualization_v(v_new,nx,ny,t);
    //visualization_p(p_new,nx,ny,t);
  }
}
