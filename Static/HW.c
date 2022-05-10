#define FILE_NAME "points.txt"
#include<mpi.h>
#include<stdio.h>
#include <stdlib.h>
#include <math.h>
#define TAG 0

double heavy(double x, double y, int param) {
	double center[2] = { 0.4, 0.2 };
	int i, loop, size = 1, coeff = 10000;
	double sum = 0, dx, dy, radius = 0.2*size;
	int longLoop = 1000, shortLoop = 1;
	double pi = 3.14;
	dx = (x - center[0]) * size;
	dy = (y - center[1]) * size;
	loop = (sqrt(dx * dx + dy * dy) < radius) ? longLoop : shortLoop;

	for (i = 1; i < loop * coeff; i++)
		sum += cos(2*pi * dy * dx + 0.1) * sin(exp(10*cos(pi* dx))) / i;

	return  sum;
}

double* readFromFile(const char* fileName, int* numberOfPoints, int *param) {
	FILE* fp;
	double* points;

	// Open file for reading points
	if ((fp = fopen(fileName, "r")) == 0) {
		printf("cannot open file %s for reading\n", fileName);
		exit(0);
	}

	// Param
	fscanf(fp, "%d", param);

	// Number of points
	fscanf(fp, "%d", numberOfPoints);

	// Allocate array of points end Read data from the file
	points = (double*)malloc(2 * *numberOfPoints * sizeof(double));
	if (points == NULL) {
		printf("Problem to allocate memory\n");
		exit(0);
	}
	for (int i = 0; i < *numberOfPoints; i++) {
		fscanf(fp, "%le %le", &points[2 * i], &points[2 * i + 1]);
	}

	fclose(fp);

	return points;
}

int main(int argc, char* argv[]) {
	double *points, x, y;
	double t1, t2; 
	double answer = 0;
	double xyArray[2];
	int param;
	int rank, size;
	int numberOfPoints = 10;
	int master = 0;
	int pi = 0;
	
	MPI_Status status;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	
	// Read points from the file
	points = readFromFile(FILE_NAME, &numberOfPoints, &param);		
	x = points[0];
	y = points[1];
	answer = heavy(x, y, param);
		
	// Process 0 is Master
	if(rank == master) {
		double tmpAnswer = 0;
				
		// Perform heavy sequential computation
		for (int i = 1; i < numberOfPoints; i++) {
			// sent to process pi = ( pi + 1 ) % size
			// if pi == master => pi++
			pi = ( pi + 1 ) % size;
			if(pi == master)
				pi = ( pi + 1 ) % size;
			
			xyArray[0] = points[2 * i];
			xyArray[1] = points[2 * i + 1];
			MPI_Send(xyArray, 2, MPI_DOUBLE, pi, TAG, MPI_COMM_WORLD);
		}
		
		t1 = MPI_Wtime();
		for(int i = 1; i < size; i++)
		{
			MPI_Recv(&tmpAnswer, 1, MPI_DOUBLE, i, TAG, MPI_COMM_WORLD, &status);		
			answer = fmax(answer, tmpAnswer);
		}
	}
	
	if(rank != master) {
		for (int i = 1; i < numberOfPoints; i++) {
			pi = ( pi + 1 ) % size;
			if(pi == master)
				pi = ( pi + 1 ) % size;
				
			if(rank == pi)
			{
				MPI_Recv(xyArray, 2, MPI_DOUBLE, master, TAG, MPI_COMM_WORLD, &status);
				// try xyArray[0] or *xyArray[0]
				answer = fmax(answer, heavy(xyArray[0], xyArray[1], param));
			}
		}
		MPI_Send(&answer, 1, MPI_DOUBLE, master, TAG, MPI_COMM_WORLD);
	}
	if(rank == master)
	{
		t2 = MPI_Wtime();
		printf("MPI_Wtime measured execution time: %1.2f\n",  t2-t1); fflush(stdout);
		printf("answer = %e\n", answer);
	}	
	MPI_Finalize( );
	return 0;
}
