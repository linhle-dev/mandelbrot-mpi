//////////////////////////////////////////////////////////////////////////////////////
// mandelbrot.c program: Mandelbort Set Fractal (Color Serial Code Implementation).
// --------------------------------
//  1. Draws Mandelbrot set for Fc(z) = z*z +c
//  using Mandelbrot algorithm ( boolean escape time )
//	This code is modified from the original version as available at:
//	http://rosettacode.org/wiki/Mandelbrot_set#PPM_non_interactive
// -------------------------------         
// 2. Technique of creating ppm file is  based on the code of Claudio Rocchini
// http://en.wikipedia.org/wiki/Image:Color_complex_plot.jpg
// create 24 bit color graphic file ,  portable pixmap file = PPM 
// see http://en.wikipedia.org/wiki/Portable_pixmap
// to see the file use external application ( graphic viewer)
//////////////////////////////////////////////////////////////////////////////////////

/*
COMMENT HEADER
Name: Linh Hoang LE
ID: 28127668
Email: llee0010@student.monash.edu
*/
#include <stdio.h> 
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <mpi.h>

int main (int argc, char** argv) {
    int rank, size;
    MPI_Status status;
	double startTime, endTime, time;

    /* screen ( integer) coordinate */
	int iX,iY;
	const int iXmax = 8000; // default
	const int iYmax = 8000; // default

	/* world ( double) coordinate = parameter plane*/
	double Cx, Cy;
	const double CxMin = -2.5;
	const double CxMax = 1.5;
	const double CyMin = -2.0;
	const double CyMax = 2.0;

	/* */
	double PixelWidth = (CxMax - CxMin)/iXmax;
	double PixelHeight = (CyMax - CyMin)/iYmax;

	/* color component ( R or G or B) is coded from 0 to 255 */
	/* it is 24 bit color RGB file */
	const int MaxColorComponentValue = 255; 
	FILE * fp;
	char *filename = "Mandelbrotmpi_2.ppm";
	char *comment = "# ";	/* comment should start with # */

	/* Z = Zx + Zy*i;	Z0 = 0 */
	double Zx, Zy;
	double Zx2, Zy2; /* Zx2 = Zx*Zx;  Zy2 = Zy*Zy  */
	/*  */
	int Iteration;
	const int IterationMax = 1000; // default

	/* bail-out value , radius of circle ;  */
	const double EscapeRadius = 400;
	double ER2 = EscapeRadius * EscapeRadius;

    // Get current clock time.
	startTime = MPI_Wtime();

    // Initialize MPI
    MPI_Init(&argc , & argv);
    MPI_Comm_size ( MPI_COMM_WORLD , & size );
    MPI_Comm_rank ( MPI_COMM_WORLD , & rank );

    //local buffer
	char *localBuffer = malloc(iXmax * 3);

    // global buffer
	char *globalBuffer = NULL;
	if (rank == 0){
		globalBuffer = malloc(iXmax * 3 * size);

		/*create new file,give it a name and open it in binary mode  */
		fp = fopen(filename, "wb"); /* b -  binary mode */

		/*write ASCII header to the file (PPM file format)*/
		fprintf(fp,"P6\n %s\n %d\n %d\n %d\n", comment, iXmax, iYmax, MaxColorComponentValue);

		printf("File: %s successfully opened for writing.\n", filename);
		printf("Computing Mandelbrot Set. Please wait...\n");

	}

    
    /* compute and write image data bytes to the file */
	for(iY = rank; iY < iYmax; iY+=size)    // 
	{
		//printf("%d\n", iY);
		Cy = CyMin + (iY * PixelHeight);
		if (fabs(Cy) < (PixelHeight / 2))
		{
			Cy = 0.0; /* Main antenna */
		}

		for(iX = 0; iX < iXmax; iX++)
		{
			unsigned char* color = localBuffer + (iX * 3); 

			Cx = CxMin + (iX * PixelWidth);
			/* initial value of orbit = critical point Z= 0 */
			Zx = 0.0;
			Zy = 0.0;
			Zx2 = Zx * Zx;
			Zy2 = Zy * Zy;

			for(Iteration = 0; Iteration < IterationMax && ((Zx2 + Zy2) < ER2); Iteration++)
			{
				Zy = (2 * Zx * Zy) + Cy;
				Zx = Zx2 - Zy2 + Cx;
				Zx2 = Zx * Zx;
				Zy2 = Zy * Zy;
			};

			/* compute  pixel color (24 bit = 3 bytes) */
			if (Iteration == IterationMax)
			{
				// Point within the set. Mark it as black
				color[0] = 0;
				color[1] = 0;
				color[2] = 0;
			}
			else 
			{
				// Point outside the set. Mark it as white
				double c = 3*log((double)Iteration)/log((double)(IterationMax) - 1.0);
				if (c < 1)
				{
					color[0] = 0;
					color[1] = 0;
					color[2] = 255*c;
				}
				else if (c < 2)
				{
					color[0] = 0;
					color[1] = 255*(c-1);
					color[2] = 255;
				}
				else
				{
					color[0] = 255*(c-2);
					color[1] = 255;
					color[2] = 255;
				}
			}

		}

		MPI_Gather(localBuffer, (iXmax*3), MPI_BYTE, globalBuffer, (iXmax*3), MPI_BYTE, 0, MPI_COMM_WORLD);

		if(rank == 0){
			fwrite(globalBuffer, iXmax*3*size, 1, fp);
		}
	}


	endTime = MPI_Wtime();
    time = endTime - startTime;

	// fclose(fp);

    if (rank == 0){	

		printf("Completed Computing Mandelbrot Set.\n");
		printf("File: %s successfully closed.\n", filename);
		printf("Mandelbrot parallel computational process time: %lf\n", time);

    }

    // Finalize the MPI
    MPI_Finalize ();


    return (0);
}