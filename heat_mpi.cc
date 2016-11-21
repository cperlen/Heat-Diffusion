#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <mpi.h>

using namespace std;

inline double square(double a){
  return(a*a);
}


int main(int argc, char *argv[]){
  const int nx = atoi(argv[1]);
  const double kappa = 1;  //not specified by assignment
  const double dx = M_PI/nx;
  const double dt = square(dx)/(5*kappa); //to satisfy numerical stability 
  const double num_steps = 0.5*square(M_PI)/(kappa*dt);

  int num_tasks;
  int rank;
  struct timeval t_start,t_end;
  gettimeofday(&t_start, NULL);
  
  MPI_Status stat[4];
  MPI_Request reqs[4];

  int rc = MPI_Init(&argc,&argv);
  if (rc != MPI_SUCCESS) {
    printf("Error starting MPI program. Terminating.\n");
    MPI_Abort(MPI_COMM_WORLD, rc);
  }

  MPI_Comm_size(MPI_COMM_WORLD, &num_tasks);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  //initialize T^n  
  int nrows = nx/num_tasks + 2; //2 ghost rows

  double T_n[nrows][nx]; //transpose representation aligns with memory structure, static allocation to avoid segfaulting
  for (int i=0;i < nrows; i++){
    for(int j=0; j < nx; j++){  
      T_n[i][j] = 0;
    }
  }

  //boundary conditions
  if(rank == 0){
    for(int i=0; i < nx; i++){
      T_n[1][i] = square(cos(i*M_PI/double(nx)));
    }
  }

  if(rank == (num_tasks-1)){
    for(int i=0; i < nx; i++){
      T_n[nx/num_tasks][i] = square(sin(i*M_PI/double(nx)));
    }
  }

  //initialize T^n+1
  double T_n1[nrows][nx]; 
  for (int i=0; i < nrows; i++){
    for(int j=0; j < nx; j++){  
      T_n1[i][j] = T_n[i][j];
    }
  }

  int t = 0;
  while(t < num_steps){
    int rowstart = 1;
    if(rank == 0){
      rowstart = 2; //avoid updating boundary conditions
    }
    int rowend = nrows - 1;
    if(rank == num_tasks -1){
      rowend = nrows - 2;
    }

    for(int i = rowstart; i < rowend; i++){
      for(int j=1; j < nx-1; j++){  
        T_n1[i][j] = T_n[i][j]+ kappa * dt/square(dx) * (T_n[i-1][j] + T_n[i+1][j] + T_n[i][j-1] + T_n[i][j+1] - 4 * T_n[i][j]);
      }
    }
      //periodicity
    for(int j = rowstart; j < rowend; j++){
      T_n1[j][nx-1] = T_n[j][nx-1] + kappa * dt/square(dx) * (T_n[j][nx-2]+T_n[j][0] + T_n[j-1][nx-1] + T_n[j+1][nx-1]-4*T_n[j][nx-1]);
      T_n1[j][0] = T_n1[j][nx-1];
    }

    for(int i=0; i < nrows; i++){
      for(int j=0; j < nx; j++){  
        T_n[i][j] = T_n1[i][j];
      }
    }

    int fwd =(rank+1) % num_tasks;
    int back =(rank-1) % num_tasks;
    while(back < 0){
      back = back + num_tasks;
    }

    rc = MPI_Isend(&T_n[nx/num_tasks][0], nx, MPI_DOUBLE, fwd , 1, MPI_COMM_WORLD, &reqs[0]);
    if (rc != MPI_SUCCESS) {
      printf("Error sending from %d to %d \n", rank, fwd);
      MPI_Abort(MPI_COMM_WORLD, rc);
    }

    rc = MPI_Isend(&T_n[1][0], nx, MPI_DOUBLE, back, 2, MPI_COMM_WORLD, &reqs[1]);
    if (rc != MPI_SUCCESS) {
      printf("Error sending from %d to %d \n", rank, back);
      MPI_Abort(MPI_COMM_WORLD, rc);
    }

    rc = MPI_Irecv(&T_n[nx/num_tasks+1][0], nx, MPI_DOUBLE, fwd, 2, MPI_COMM_WORLD, &reqs[2]);
    if (rc != MPI_SUCCESS) {
      printf("Error receiving from %d at %d \n", rank, fwd);
      MPI_Abort(MPI_COMM_WORLD, rc);
    }

    rc = MPI_Irecv(&T_n[0][0], nx, MPI_DOUBLE, back, 1, MPI_COMM_WORLD, &reqs[3]);
    if (rc != MPI_SUCCESS) {
      printf("Error receiving from %d at %d \n", rank, back);
      MPI_Abort(MPI_COMM_WORLD, rc);
    }
    MPI_Waitall(4, reqs, stat);
    t += 1;
  }
  double tot = 0.0;
  for(int i = 1; i < nrows - 1; i++){ //not including ghost columns
    for(int j =0; j < nx; j++){
      tot += T_n[i][j];
    }
  }
  MPI_Request reqs2[2*(num_tasks-1)];
  MPI_Status stat2[num_tasks];
  if(rank != 0){
    rc = MPI_Isend(&tot, 1, MPI_DOUBLE, 0, rank, MPI_COMM_WORLD, &reqs2[rank-1]);
    if (rc != MPI_SUCCESS) {
        printf("Error receiving total value from %d", rank);
        MPI_Abort(MPI_COMM_WORLD, rc);
      }
  }
  
  if(rank == 0){
    double totals[num_tasks];
    for(int i = 1; i < num_tasks; i++){
      rc = MPI_Irecv(&totals[i], 1, MPI_DOUBLE, i, i, MPI_COMM_WORLD, &reqs2[(num_tasks-1)+i]);
      if (rc != MPI_SUCCESS) {
        printf("Error receiving total value from %d", i);
        MPI_Abort(MPI_COMM_WORLD, rc);
      }
    }
    MPI_Waitall(num_tasks-1, &reqs2[num_tasks], stat2);
    for(int i=1; i < num_tasks; i++){
      tot += totals[i];
    }
  
    gettimeofday(&t_end,NULL);
    float runtime = ((t_end.tv_sec -t_start.tv_sec) * 1000000u + t_end.tv_usec - t_start.tv_usec) / 1.e6;

    cout <<  "Total time elapsed "   << runtime << endl;
    cout << "The average temperature is " << tot/square(nx) << endl;
  }

  char filename[50];
  sprintf(filename,"map_mpi%d.txt",rank);
  ofstream fout(filename);
  for(int i=0;i<nrows;i++){
   for(int j=0;j<nx;j++){ 
     fout<< i<<" "<<j<<" "<< T_n[i][j]<<endl;
   }fout<<endl;
 }

 fout.close();
 MPI_Finalize();
}