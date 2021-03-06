#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <time.h>
using namespace std;

inline double square(double x){
  return(x*x);
}

int main(int argc, char *argv[]){

//start clock

  double elapsed_t = clock()/double(CLOCKS_PER_SEC);

  if (argc != 2) {
    printf("USAGE: %s <nx>\n", argv[0]);
    exit(1);
  }

  const int nx = atoi(argv[1]);
  const double kappa = 1;  //not specified by assignment
  const double dx = M_PI/nx;
  const double dt = square(dx)/(5*kappa); //to satisfy numerical stability 
  const double num_steps = 0.5*square(M_PI)/(kappa*dt);

//initialize T^n
  double** T_n = new double*[nx];
  for (int i=0; i < nx; i++){
    T_n[i] = new double[nx];
    for(int j=0; j < nx; j++){ 
      T_n[i][j] = 0;
    }
  } 

//boundary conditions 
  for (int i=0; i < nx; i++){
    T_n[i][0] = square(cos(i * M_PI/double(nx)));  //T(x,0) = cos^2(x)
    T_n[i][nx-1] = square(sin(i * M_PI/double(nx))); //T(x,\pi) = sin^2(x)
  }

//initialize T^n+1
  double** T_n1 = new double*[nx];
  for (int i=0; i < nx; i++){
    T_n1[i] = new double[nx];
    for(int j=0; j < nx; j++){ 
      T_n1[i][j] = T_n[i][j];
    }
  }

//finite difference simulation
  int t = 0;
  while(t < num_steps){
    for(int i=1; i < nx-1; i++){
      for(int j=1; j < nx-1; j++){ 
        T_n1[i][j] = T_n[i][j] + kappa * dt/square(dx) * (T_n[i-1][j]+T_n[i+1][j]+T_n[i][j-1]+T_n[i][j+1]-4*T_n[i][j]);
      }
    }

  //periodicity
    for(int j=1; j < nx-1; j++){
      T_n1[nx-1][j] = T_n[nx-1][j] + kappa * dt/square(dx) * (T_n[nx-2][j]+T_n[0][j]+T_n[nx-1][j-1]+T_n[nx-1][j+1]-4*T_n[nx-1][j]);
      T_n1[0][j] = T_n1[nx-1][j];
    }

    for(int i=0; i < nx; i++){
      for(int j=0; j < nx; j++){ 
        T_n[i][j] = T_n1[i][j];
      }
    }
    t += 1;
  }

//stop time
  elapsed_t = clock()/double(CLOCKS_PER_SEC) - elapsed_t;


//print output
  double tot = 0.0;
  char filename[50];
  sprintf(filename,"heatmap_serial_%d.txt",nx);
  ofstream fout(filename);

  for(int i = 0; i < nx; i++){
    for(int j = 0; j < nx; j++){ 
      fout<< i<<" "<<j<<" "<< T_n[i][j]<<endl;
      tot += T_n[i][j]; //average
    }
    fout<<endl;
  }
  fout.close();

  cout<<"Total time elapsed "<< elapsed_t <<endl;
  cout<<"Volume average temperature is "<<tot/(square(nx))<<endl;
  
  //free memory
  for(int i = 0; i < nx; i++){
    delete[] T_n[i];
    delete[] T_n1[i];
  }
  delete[] T_n;
  delete[] T_n1;
}
