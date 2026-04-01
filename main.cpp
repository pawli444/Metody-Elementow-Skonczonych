int send[4] = {0,1,2,3};
int recv[4];
MPI_Alltoall(send, 1, MPI_INT,recv, 1, MPI_INT,MPI_COMM_WORLD);
