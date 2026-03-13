#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "ranlxd.h"

/* Definition of globals */
#ifdef CONTROL
#define EXTERN
#else
#define EXTERN extern
#endif

/* graphics details */
#define WINDOW1_SIZEX 1420			/* window sizes in pixels */
#define WINDOW1_SIZEY 1024
#define WINDOW2_SIZEX 500
#define WINDOW2_SIZEY 900
#define PIXEL_SIZE 0.002
#define MAX_X 1.0				/* window coordinate */
#define MAX_Y 1.0
#define BAR_X 0.2				/* histogram bar size */
#define BAR_Y 0.2
#define VISIBLE_POINTS 100

EXTERN int win1,win2;					/* window identifier */
EXTERN clock_t press_time;				/* mouse press time */
EXTERN int paused;
EXTERN int refresh_period;				/* in microseconds */
EXTERN int invert_color;

#define NORMAL 0
#define NARROW 1
#define MULTIMODAL 2
#define WALL 3
EXTERN int distr;
EXTERN double max_distr;
EXTERN double meanx, meany, sigx, sigy;		/* distributions */
EXTERN double meanx2, meany2, sigx2, sigy2;

#define NAIVE_UNIFORM 0
#define NAIVE_DISTR 1
#define MARKOV 2
EXTERN int mode;
EXTERN char title [256];

/* points */
struct point {
	double x,y;
};
EXTERN struct point *points;
EXTERN int Npoints;

EXTERN struct point current_point, proposed_point;

EXTERN double Norm;
EXTERN double Xtot, Ytot, Xavg, Yavg;
EXTERN int N_acc, N_rej;

EXTERN double step_size;				/* Metropolis step size */

/* histograms */
EXTERN int hist_nx, hist_ny;
EXTERN double *hist_x, *hist_y;
#define NHIST 150

/* history */
#define HISTORY_LENGTH 200			/* plotting history */
struct History {
	double Xavg, Yavg, acceptance;
};
EXTERN struct History history[HISTORY_LENGTH];

/* 1d density functions and average values */
#define NDENSE 500
EXTERN double *dfx, *dfy;
EXTERN double Xavg_analytic, Yavg_analytic;

/* graphics_utils.c */
void init_openGL(int argc, char *argv[] );
void start_openGL();
void keyboard(unsigned char key, int x, int y);
void mouse(int button, int state, int x, int y);
void display_graph ();
void display_distro ();
void update_displays ();

/* markov.c */
void pick_random_point ( int type );
double density_function ( double x, double y );
double density_function_x ( double x );
double density_function_y ( double y );
void metropolis_step ( );
void fill_history ();
void help();

/* rng_gauss.c */
double rng_gauss_single(double s);
double rng_gauss_single2(double sigma);
double rng_gauss_array(double g[], int n, double s);
double ran3();

