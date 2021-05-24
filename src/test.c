#include <stdio.h>
#include <sys/time.h>
#include <time.h>
#include <stdlib.h>
#include <ctype.h>
#include <getopt.h>
#include <string.h>

#include <poissonsolver.h>
#include <example.h>

int main(int argc, char **argv){
   
	unsigned int use_multigrid = 0;	// Whether to use multigrids or not
	double reltol = 1e-5;			// Relative tolerance
    unsigned int from = 100;
    unsigned int to = 1000;
    unsigned int step = 100;

	int nm_flags = NM_NONE;
	unsigned int p = 0;				// Problem to solve: 0=Mixed, 1=Dirichlet1, 2=Dirichlet2, 3=Neumann1, 4=Neumann2
	double (*dirichlet)(double x, double y);
	double (*g)(double x, double y);

	char filename[50] = "";
  	int c;
	// Process the options
	while ((c = getopt(argc, argv, "f:t:s:mr:n:p:")) != -1)
		switch (c)
		{
        case 'f':
			from = atoi(optarg);
			break;
        case 't':
			to = atoi(optarg);
			break;
        case 's':
			step = atoi(optarg);
			break;
		case 'p':
			p = atoi(optarg);
			break;
		case 'm':
			use_multigrid = 1;
			break;
		case 'n':
			strcpy(filename,optarg);
			break;
		case 'r':
			reltol = atof(optarg);
			if(reltol <= 0.0){
				fprintf(stderr,"Invalid argument %s for reltol. Needs to be a double > 0\n", optarg);
				return 1;
			}
			break;
		
		case '?':
			if (optopt == 'r')
				fprintf (stderr, "Option -%c requires an argument.\n", optopt);
			else if (isprint (optopt))
				fprintf (stderr, "Unknown option `-%c'.\n", optopt);
			else
				fprintf (stderr,"Unknown option character `\\x%x'.\n", optopt);
			return 1;
		
		}
	
	switch (p)
	{
	case 0:
		dirichlet = &mix;
		g = NULL;
		nm_flags = NM_X0 | NM_Y0;
		break;
	case 1:
		dirichlet = &phi;
		g = &g1;
		break;
	case 2:
		dirichlet = &phi;
		g = &g2;
		break;
	case 3:
		dirichlet = NULL;
		g = &g1;
		nm_flags = NM_ONLY;
		break;
	case 4:
		dirichlet = NULL;
		g = &g2;
		nm_flags = NM_ONLY;
		break;
	
	default:
	return 1;
		break;
	}

	FILE *f;
	c = strcmp(filename, "");
	if(c != 0){
			f = fopen(filename, "w");
	}
   
    for (int n = from; n <= to; n+=step)
    {
        /* Create a poisson boundary value problem */
        bvp_t bvp = bvp_create(n, dirichlet, g, nm_flags);
        
        /* Start timer */
        struct timespec start, end;
        clock_gettime(CLOCK_MONOTONIC, &start);
        
        /* Solve the BVP */
        int i = solve_poisson_bvp(bvp, use_multigrid, reltol);
        
        /* End timer */
        clock_gettime(CLOCK_MONOTONIC, &end);
        double secs = (double) (end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec) / 1000000000;
		printf("%d %f %d\n",n,secs,i);
		fflush(stdout);
        
		if(c != 0){
			fprintf(f,"%d %f %d\n",n,secs,i);
			fflush(f);
		}

		bvp_destroy(bvp);
    }
	if(c != 0){
		 fclose(f);
	}

	return 0;
}