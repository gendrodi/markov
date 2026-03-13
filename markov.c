/* Stochastic integration in two dimensions	*/
/* written by Gergely Endrodi			*/

#define CONTROL

#include "includes.h"

double density_function ( double x, double y )
{
	switch ( distr )
	{
		case NORMAL:
		case NARROW: 
			return ( exp( - pow((x-meanx)/sigx,2.0)/2. ) * 
				 exp( - pow((y-meany)/sigy,2.0)/2. ) );
			break;
		case MULTIMODAL:
			return ( 0.5 * exp( - pow((x-meanx)/sigx,2.0)/2. ) * 
				       exp( - pow((y-meany)/sigy,2.0)/2. ) +
				 0.5 * exp( - pow((x-meanx2)/sigx2,2.0)/2. ) * 
				       exp( - pow((y-meany2)/sigy2,2.0)/2. ) );
			break;
		case WALL:
			if ( (x>0.47)&&(x<0.57) ) return ( 0. );
			return ( 0.5 * exp( - pow((x-meanx)/sigx,2.0)/2. ) * 
				       exp( - pow((y-meany)/sigy,2.0)/2. ) +
				 0.5 * exp( - pow((x-meanx2)/sigx2,2.0)/2. ) * 
				       exp( - pow((y-meany2)/sigy2,2.0)/2. ) );
			break;
	}
}

double fill_density_functions ( )
{
	double res, resZ, x, dx, xmin, xmax;
	double dy,y,ymin,ymax;
	ymin = xmin = 0.;
	ymax = xmax = 1.;
	dx = dy = 1.e-3;
	int i;

	dfx = (double*)malloc( NDENSE * sizeof(double) );
	dfy = (double*)malloc( NDENSE * sizeof(double) );
	double m1=-1, m2=-1;

	fprintf(stderr," Calculating 1d density functions ...");

	x = xmin;
	for (i=0; i<NDENSE; i++)
	{
		res = 0.0;
		for (y=ymin; y<=ymax; y+=dy) res += dy * density_function( x , y );
		dfx[i] = res;
		x += (xmax-xmin)/(double)NDENSE;
		if (m1<dfx[i]) m1 = dfx[i];
	}

	y = ymin;
	for (i=0; i<NDENSE; i++)
	{
		res = 0.0;
		for (x=xmin; x<=xmax; x+=dx) res += dx * density_function( x , y );
		dfy[i] = res;
		y += (ymax-ymin)/(double)NDENSE;
		if (m2<dfy[i]) m2 = dfy[i];
	}

	for (i=0; i<NDENSE; i++) { dfx[i] *= (1./m1); dfy[i] *= (1./m2); }

	/* average values */
	res = resZ = 0.;
	for (i=0; i<NDENSE; i++) { res += dfx[i] * (xmin+(xmax-xmin)*(double)i/(double)NDENSE); resZ += dfx[i]; }
	Xavg_analytic = res / resZ;

	res = resZ = 0.;
	for (i=0; i<NDENSE; i++) { res += dfy[i] * (ymin+(ymax-ymin)*(double)i/(double)NDENSE); resZ += dfy[i]; }
	Yavg_analytic = res / resZ;

	fprintf(stderr," done\n");
}

double density_function_x ( double x )
{
	return ( dfx [ (int)floor(x * (double)NDENSE + 0.5) ] );

	/* the below are incorrect for MULTIMODAL	*/
	/* which does not factorize!			*/
	switch ( distr )
	{
		case NORMAL: 
		case NARROW: 
			return ( exp( - pow((x-meanx)/sigx,2.0)/2. ) );
			break;
		case MULTIMODAL:
			return ( 0.85 * exp( - pow((x-meanx)/sigx,2.0)/2. ) + 
				 0.85 * exp( - pow((x-meanx2)/sigx2,2.0)/2. ) );
			break;
	}
}

double density_function_y ( double y )
{
	return ( dfy [ (int)floor(y * (double)NDENSE + 0.5) ] );

	/* the below are incorrect for MULTIMODAL	*/
	/* which does not factorize!			*/
	switch ( distr )
	{
		case NORMAL:
		case NARROW: 
			return ( exp( - pow((y-meany)/sigy,2.0)/2. ) );
			break;
		case MULTIMODAL:
			return ( 0.85 * exp( - pow((y-meany)/sigy,2.0)/2. ) + 
				 0.85 * exp( - pow((y-meany2)/sigy2,2.0)/2. ) );
			break;
	}
}

double get_max_distr ()
{
	switch ( distr )
	{
		case NORMAL: 
		case NARROW: 
			return ( density_function( meanx, meany) );
			break;
		case MULTIMODAL:
		case WALL: 
			return ( density_function( meanx, meany) );
			break;
	}
}

void update_histogram ( double x, double *hist, int hist_n, double delta )
{
	int n = (int) (x * hist_n);
	if (n>=hist_n) n = hist_n-1;
	if (n<0) n = 0;
	hist[ n ] += delta;
}

// type = 0 : uniform
// type = 1 : gaussian
// type = 2 : corner
void pick_random_point ( int type )
{
	if ( Npoints < VISIBLE_POINTS ) points = (struct point*) realloc( points, (Npoints+1) * sizeof(struct point));

	int np = Npoints % VISIBLE_POINTS;
	Npoints++;

	points[np].x = (type==0) ? ran3() : ( ( type==1) ? (meanx + rng_gauss_single( sigx ) ) : 0.05 );
	points[np].y = (type==0) ? ran3() : ( ( type==1) ? (meany + rng_gauss_single2( sigy ) ) : 0.05 );

	update_histogram ( points[np].x, hist_x, hist_nx, (type==0) ? density_function ( points[np].x, points[np].y ) : 1.0 );
	update_histogram ( points[np].y, hist_y, hist_ny, (type==0) ? density_function ( points[np].x, points[np].y ) : 1.0 );

	Xtot += points[np].x * ( (type==0) ? density_function ( points[np].x, points[np].y ) : 1.0 );
	Ytot += points[np].y * ( (type==0) ? density_function ( points[np].x, points[np].y ) : 1.0 );
	Norm += (type==0) ? density_function ( points[np].x, points[np].y ) : 1.0;
}

void fill_history ()
{
	int i;
	struct History event;
	event.Xavg = Xavg;
	event.Yavg = Yavg;
	event.acceptance = (double)N_acc/(double)(N_acc+N_rej);
	
	/* fill graph for plotting */
	if ( Npoints < HISTORY_LENGTH ) 
	{
		history[Npoints] = event;
	}
	else
	{
		if (!paused) for (i=0; i<HISTORY_LENGTH-1; i++) 
		{
			history[i] = history[i+1];
		}
		history[HISTORY_LENGTH-1] = event;
	}
}

void metropolis_step ( )
{
	if ( Npoints < VISIBLE_POINTS ) points = (struct point*) realloc( points, (Npoints+1) * sizeof(struct point));

	int np = Npoints % VISIBLE_POINTS;
	Npoints++;

	//double x = ran3();
	//double y = ran3();

	double x = (current_point.x + (ran3()-.5) * step_size );
	double y = (current_point.y + (ran3()-.5) * step_size );

	// reflection
	if (x<0) x=-x; if (x>MAX_X) x=2*MAX_X-x;
	if (y<0) y=-y; if (y>MAX_Y) y=2*MAX_Y-y;

	double expmindeltaH = density_function ( x, y ) / density_function ( current_point.x, current_point.y );

	/* jump with probability P */
	if ( ran3() < expmindeltaH ) { current_point.x = x; current_point.y = y; N_acc++; } else N_rej++;
	proposed_point.x = x; proposed_point.y = y;

	points[np] = current_point;
	update_histogram ( points[np].x, hist_x, hist_nx, 1.0 );
	update_histogram ( points[np].y, hist_y, hist_ny, 1.0 );

	Xtot += points[np].x;
	Ytot += points[np].y;
	Norm++;
}

void free_arrays ()
{
	free ( points );
	free ( hist_x );
	free ( hist_y );
	free ( dfx );
	free ( dfy );
}

int init_parameters (int argc, char *argv[] )
{
	int usage()
	{
		fprintf(stderr,"\nVisualization of sampling methods\n");
		fprintf(stderr,"  USAGE: %s mode [distribution]\n", argv[0]);
		fprintf(stderr,"  mode = 0: uniformly distributed random points\n");
		fprintf(stderr,"       = 1: Gaussian distributed random points\n");
		fprintf(stderr,"       = 2: Markov chain (Metropolis algorithm)\n");
		fprintf(stderr,"  distribution = 0: normal (default)\n");
		fprintf(stderr,"               = 1: narrow\n");
		fprintf(stderr,"               = 2: multi-modal\n");
		fprintf(stderr,"               = 3: multi-modal with wall\n");
		fprintf(stderr,"  [mode=1 only works with normal/narrow distribution]\n\n");
                fprintf(stderr,"  press 'h' for help\n\n");
		fflush(stderr);
		return(-1);
	}
	
	if ( (argc<2) || ( (argc==2) && (strcmp(argv[1],"--help")==0) ) ) return( usage() );
	mode = atoi ( argv[1] );
	distr = NORMAL;
	if (argc>2) distr = atoi( argv[2] );

	switch ( mode )
	{
		case NAIVE_UNIFORM: 
			sprintf(title,"Sampling with uniformly distributed random points" ); 
			break;
		case NAIVE_DISTR:
			sprintf(title,"Sampling with Gaussian distributed random points" ); 
			if ( ( distr != NORMAL ) && ( distr != NARROW ) ) {
				usage();
				exit(0);
			}
			break;
		case MARKOV:
			sprintf(title,"Sampling with a Markov chain (Metropolis algorithm)" ); 
			break;
	}

	switch ( distr )
	{
		case NORMAL: 
			meanx = 0.5;  meany = 0.6;  sigx = 0.075; sigy = 0.15;
			Xavg_analytic = meanx;
			Yavg_analytic = meany;
			break;
		case NARROW: 
			meanx = 0.5;  meany = 0.6;  sigx = 0.01; sigy = 0.04;
			Xavg_analytic = meanx;
			Yavg_analytic = meany;
			break;
		case MULTIMODAL:
		case WALL:
			//meanx = 0.25;  meany = 0.2;  sigx = 0.1;   sigy = 0.05;
			//meanx2 = 0.7; meany2 = 0.5; sigx2 = 0.05; sigy2 = 0.15;
			meanx = 0.25;  meany = 0.6;  sigx = 0.04;   sigy = 0.1;
			meanx2 = 0.75; meany2 = 0.6; sigx2 = 0.04; sigy2 = 0.1;
			Xavg_analytic = 0.;
			Yavg_analytic = 0.;
			break;
	}

	max_distr = get_max_distr();

	hist_nx = NHIST;
	hist_x = (double*) malloc(hist_nx * sizeof(double) );
	hist_ny = NHIST;
	hist_y = (double*) malloc(hist_ny * sizeof(double) );

	Npoints = 0;
	points = (struct point*) malloc( (Npoints+1) * sizeof(struct point));

	Norm = Xtot = Ytot = 0.0;
	N_acc = N_rej = 0;
	Xavg = Yavg = 0.0;

	paused = 1;
	refresh_period = 512;
		
	step_size = 0.2;
	
	invert_color = 0;
	
	return(0);
}

void help()
{
	fprintf(stderr," h         : output help\n");
	fprintf(stderr," a         : increase stepsize by 0.05\n");
	fprintf(stderr," z         : decrease stepsize by 0.05\n");
	fprintf(stderr," s         : faster");
	fprintf(stderr," x         : slower");
	fprintf(stderr," n         : reset measurements\n");
	fprintf(stderr," CTRL+n    : reset measurements and histograms\n");
	fprintf(stderr," CTRL+m    : reset measurements and histograms, start from corner\n");
	fprintf(stderr," i         : invert colors\n");
	fprintf(stderr," Space     : pause\n");
	fprintf(stderr," Esc/Enter : quit\n");
	fflush(stderr);
}


int main(int argc, char *argv[])
{
	/* set parameters */
	if ( init_parameters (argc, argv ) ) exit(-1);

	/* init random number generator */
	rlxd_init(1,time(NULL));

	/* calculate expected probability densities */
	fill_density_functions ( );
	
	/* pick starting point randomly */
	pick_random_point ( 0 );
	current_point.x = points[0].x; 
	current_point.y = points[0].y;

	/* init openGL stuff */
	init_openGL(argc, argv );

	/* register cleanup function at exit */
	atexit( free_arrays );

	/* start graphics */
	start_openGL();

	return(0);
}

