#include "includes.h"
#include <GL/gl.h>
#include <GL/glut.h>
#include <time.h>

void init_openGL(int argc, char *argv[]) 
{
	/* openGL initialization */
	glutInit(&argc, argv); 
	glutInitDisplayMode (GLUT_SINGLE | GLUT_RGB);
	
	/* first window, displaying the distributions */
	glutInitWindowSize ( WINDOW1_SIZEX, WINDOW1_SIZEY); 
	glutInitWindowPosition (00, 00);
	glutCreateWindow(""); 
	win1 = glutGetWindow();

	/* initialize viewing values  */
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(-BAR_X, MAX_X, -BAR_Y, MAX_Y, -1.0, 1.0);

	/* register keyboard and mouse utilities */
	glutKeyboardFunc (keyboard);

	/* assign display function */
	glutDisplayFunc(display_distro); 

	/* needed for transparent colors */
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	/* nice dots */
	glEnable(GL_POINT_SMOOTH);
	glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);

	/* second window, displaying the graph */
	glutInitWindowSize ( WINDOW2_SIZEX, WINDOW2_SIZEY); 
	glutInitWindowPosition (00+WINDOW1_SIZEX, 00);
	glutCreateWindow("");
	win2 = glutGetWindow();
	
	/* initialize viewing values  */
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(-0.11, 1.0, 0.0, 3.01, -1.0, 1.0);
	
	/* register keyboard and mouse utilities */
	glutKeyboardFunc (keyboard);

	/* assign display function */
	glutDisplayFunc(display_graph); 

	/* updating routine for both windows */
	glutTimerFunc( refresh_period, update_displays, 0 );
}

void start_openGL()
{
	glutMainLoop();
}

void keyboard(unsigned char key, int x, int y)
{
	int mod = glutGetModifiers();
	if (mod==2) key += 96;
        else if (mod!=0) return;

	switch (key)
	{
		case 27: case 13:		// ESCAPE or ENTER key
			fprintf(stderr,"\n");
			exit (0);
			break;
		case 32:			// SPACE
			paused = !paused;
			break;
		case 'a':
			if (step_size<0.95) step_size += 0.05;
			break;
		case 'z':
			if (step_size>0.051) step_size -= 0.05;
			break;
		case 'n':
			Npoints = 0;
			free(points);
			points = (struct point*) malloc( (Npoints+1) * sizeof(struct point));
			Norm = Xtot = Ytot = 0.0;
			N_acc = N_rej = 0;
			if (mod==2) {
				for (int i=0; i<hist_nx; i++) hist_x[i] = 0.0;
				for (int i=0; i<hist_ny; i++) hist_y[i] = 0.0;
				pick_random_point ( 0 );
				current_point.x = points[0].x; 
				current_point.y = points[0].y;
			}
			break;
		case 'm':
			Npoints = 0;
			free(points);
			points = (struct point*) malloc( (Npoints+1) * sizeof(struct point));
			Norm = Xtot = Ytot = 0.0;
			N_acc = N_rej = 0;
			if (mod==2) {
				for (int i=0; i<hist_nx; i++) hist_x[i] = 0.0;
				for (int i=0; i<hist_ny; i++) hist_y[i] = 0.0;
				pick_random_point ( 2 );
				current_point.x = points[0].x; 
				current_point.y = points[0].y;
			}
			break;
		case 'h':
			help();
			break;
		case 's':
			refresh_period /= 2;
			break;
		case 'x':
			refresh_period *= 2; if (refresh_period==0) refresh_period=1;
			break;
		case 'i':
			invert_color = ( invert_color + 1 ) % 2;
			break;
	}
}

/* font size */
#define GL_FONT GLUT_BITMAP_9_BY_15
//#define GL_FONT GLUT_BITMAP_HELVETICA_18
#define GL_FONT_WIDTH 9
#define GL_FONT_HEIGHT 15

void drawBitmapText(char *string, double x, double y, double z )
{  
	char *c;
	glRasterPos3f( x, y, z );
	for (c=string; *c != 0; c++) glutBitmapCharacter(GL_FONT, *c);
}

void draw_background()
{
	/* clear background */
	glClearColor ((invert_color)?1.0:0.0, (invert_color)?1.0:0.0, (invert_color)?1.0:0.0, 1.0);
	glClear (GL_COLOR_BUFFER_BIT);
	
	/* axes */
	glColor3f((invert_color)?1-0.7:0.7,(invert_color)?1-0.7:0.7,(invert_color)?1-0.7:0.7);
	glLineWidth(2.0); 
	
	glBegin(GL_LINES);
		glVertex2f (-0.04*BAR_X, MAX_Y);
		glVertex2f (-0.04*BAR_X, -0.04*BAR_Y);
	glEnd();
	glBegin(GL_LINES);
		glVertex2f (-0.04*BAR_X, -0.04*BAR_Y);
		glVertex2f (MAX_X, -0.04*BAR_Y);
	glEnd();

	/* expected distribution functions */
	glColor3f ((invert_color)?1-.5:.5, (invert_color)?1-0.5:0.5, (invert_color)?1-0.6:0.6);
	GLfloat x, y, xp, yp, dx, dy;
	dx = dy = PIXEL_SIZE;

	x = 0.;
	y = -0.05*BAR_X -0.95*BAR_X * density_function_x ( x );
	for (x=0; x<1; x+=dx)
	{
		yp = y;
		y = -0.05*BAR_X -0.95*BAR_X * density_function_x ( x );
		
		glBegin(GL_LINES);
			glVertex3f (x-dx, yp, 0.0);
			glVertex3f (x, y, 0.0);
		glEnd();
	}

	y = 0.;
	x = -0.05*BAR_Y -0.95*BAR_Y * density_function_y ( y );
	for (y=0; y<1; y+=dy)
	{
		xp = x;
		x = -0.05*BAR_Y -0.95*BAR_Y * density_function_y ( y );
		
		glBegin(GL_LINES);
			glVertex3f (xp, y-dy, 0.0);
			glVertex3f (x, y, 0.0);
		glEnd();
	}
}

double max_entry ( int hist_n, double *hist )
{
	int i;
	double k=1.0;
	for (i=0; i<hist_n; i++) if (hist[i]>k) k=hist[i];
	return(k);
}

void draw_histograms ()
{
	int i;
	double max_x, max_y;
	GLfloat x,y,dx,dy;
	dx = 1.0 / (double)hist_nx;
	dy = 1.0 / (double)hist_ny;
	max_x = max_entry(hist_nx, hist_x);
	max_y = max_entry(hist_ny, hist_y);

	glColor4f((invert_color)?1-0.7:0.7,(invert_color)?1-0.4:0.4,(invert_color)?1-0.1:0.1,.8);
	for (i=0; i<hist_nx; i++)
	{
		x = i * dx;
		y = -0.05*BAR_X -0.95*BAR_X * hist_x[i] / max_x;
		glRectf(x, -.01, x+dx, y);
	}

	for (i=0; i<hist_ny; i++)
	{
		y = i * dy;
		x = -0.05*BAR_Y -0.95*BAR_Y * hist_y[i] / max_y;
		glRectf(-.01, y, x, y+dy );
	}
}

void draw_points ()
{
	int p, pp;

	glColor3f((invert_color)?1-0.7:0.7,(invert_color)?1-0.4:0.4,(invert_color)?1-0.1:0.1);
	glPointSize(5.);

	int Np = (Npoints<VISIBLE_POINTS) ? Npoints : VISIBLE_POINTS;
	for (p=0; p<Np; p++)
	{
		glBegin(GL_POINTS);
			glVertex2f(points[p].x, points[p].y );
		glEnd();
	}

	/* path up to current point */
	if ( ( mode == MARKOV ) && ( Npoints>2 ) )
	{
		glColor3f((invert_color)?1:0.0,(invert_color)?1-0.8:0.8,(invert_color)?1:0.);
		for (p=Npoints-2; (p>Npoints-VISIBLE_POINTS)&&(p>0); p--)
		{
			pp = p-1;
			glBegin(GL_LINES);
				glVertex2f (points[pp % VISIBLE_POINTS].x, points[pp % VISIBLE_POINTS].y );
				glVertex2f (points[p  % VISIBLE_POINTS].x, points[p  % VISIBLE_POINTS].y);
			glEnd();
		}
	}

	/* arrow to last point */
	if ( ( mode == MARKOV ) && ( Npoints>1 ) )
	{
		p = (Npoints-1) % VISIBLE_POINTS;
		pp = (p==0) ? (VISIBLE_POINTS-1) : (p-1);

		if ( pow( proposed_point.x - current_point.x , 2.0 ) + pow( proposed_point.y - current_point.y , 2.0 ) > 1.e-4 )
			glColor3f((invert_color)?1-0.8:0.8,(invert_color)?1:0.,(invert_color)?1:0.);		//rejected
		else
			glColor3f((invert_color)?1:0.,(invert_color)?1-0.8:0.8,(invert_color)?1:0.);		//accepted

		glBegin(GL_LINES);
			glVertex2f (points[pp].x, points[pp].y );
			glVertex2f (proposed_point.x, proposed_point.y);
		glEnd();

		double ax = points[pp].x + (proposed_point.x-points[pp].x)*.9;
		double ay = points[pp].y + (proposed_point.y-points[pp].y)*.9;
		double dx = (proposed_point.x-ax);
		double dy = (proposed_point.y-ay);
		double ux = ax - dy;
		double uy = ay + dx;
		glBegin(GL_LINES);
			glVertex2f (ux, uy );
			glVertex2f (proposed_point.x, proposed_point.y);
		glEnd();
		ux = ax + dy;
		uy = ay - dx;
		glBegin(GL_LINES);
			glVertex2f (ux, uy );
			glVertex2f (proposed_point.x, proposed_point.y);
		glEnd();
	}
}

void display_distro()
{
	/* set window title */
	char s[128]={}, s2[1024];
	if ( mode == MARKOV ) sprintf(s,"Step size %.2f",step_size);
	sprintf(s2,"%s %s Refresh period: %d microsec.", title, s, refresh_period );
	glutSetWindowTitle(s2);

	draw_background();

	GLfloat x,y,xp,yp, dx, dy;
	dx = PIXEL_SIZE; dy = PIXEL_SIZE;

	GLfloat r,g,b;
	r = 0; g = 0;
	double b1, b2;

	for (x=0.; x<1.; x+=dx) for (y=0.; y<1.; y+=dy)
	{
		b = density_function ( x, y ) * 0.9 / max_distr;
		//glColor3f ((invert_color)?1-r:r, (invert_color)?1-(b*0.8):(b*.8), (invert_color)?1-b:b);
		//glColor3f ((invert_color)?1-(b):(b), (invert_color)?1-r:r, (invert_color)?1-b:b);
		glColor3f ((invert_color)?1-(b):(b), (invert_color)?1-b:b, (invert_color)?1-b:b);
		glRectf(x, y, x+dx, y+dy);
	}

	draw_points();

	draw_histograms();

	glFlush();
}

void draw_grid ( double offset )
{
	GLfloat x;
	glClearColor ((invert_color)?1:0.0, (invert_color)?1:0.0, (invert_color)?1:0.0, 1.0);
	glClear (GL_COLOR_BUFFER_BIT);
	
	glColor3f((invert_color)?1-0.7:0.7,(invert_color)?1-0.7:0.7,(invert_color)?1-0.7:0.7);
	glLineWidth(8.0); 
	
	for (x=0.+(mode!=MARKOV); x<=3.0; x+=1.0) 
	{
		glBegin(GL_LINES);
			glVertex3f (0., x, 0.0);
			glVertex3f (1., x, 0.0);
		glEnd();
	}

	glLineWidth(1.0); 

	for (x=0.+(mode!=MARKOV); x<=3.0; x+=0.25) 
	{
		glBegin(GL_LINES);
			glVertex3f (0., x, 0.0);
			glVertex3f (1., x, 0.0);
		glEnd();
	}

	for (x=0.-offset; x<=1.0; x+=0.2) 
	{
		if (x<0.0) continue;
		glBegin(GL_LINES);
			glVertex3f (x, 0.+(mode!=MARKOV), 0.0);
			glVertex3f (x, 3., 0.0);
		glEnd();
	}

	/* averages */
	glLineWidth(1.0); 
	glColor3f((invert_color)?1-0.4:0.4,(invert_color)?1-0.4:0.4,(invert_color)?1-0.4:0.4);
	glBegin(GL_LINES);
		glVertex3f (0., Xavg_analytic+2., 0.0);
		glVertex3f (1., Xavg_analytic+2., 0.0);
	glEnd();
	glBegin(GL_LINES);
		glVertex3f (0., Yavg_analytic+1., 0.0);
		glVertex3f (1., Yavg_analytic+1., 0.0);
	glEnd();
	

	/* labels */
	glColor3f((invert_color)?1-0.7:0.7,(invert_color)?1-0.7:0.7,(invert_color)?1-0.7:0.7);
	drawBitmapText("<x>", -.1, 2.75 , 0. );
	drawBitmapText("0", -.05, 2.02 , 0. );
	drawBitmapText("0.5", -.09, 2.5 , 0. );
	drawBitmapText("1", -.05, 2.98 , 0. );
	drawBitmapText("<y>", -.1, 1.75 , 0. );
	drawBitmapText("0", -.05, 1.02 , 0. );
	drawBitmapText("0.5", -.09, 1.5 , 0. );
	drawBitmapText("1", -.05, 1.98 , 0. );
	if (mode==MARKOV) {
	drawBitmapText("acc", -.1, 0.75 , 0. );
	drawBitmapText("0", -.05, 0.02 , 0. );
	drawBitmapText("0.5", -.09, 0.5 , 0. );
	drawBitmapText("1", -.05, 0.98 , 0. );
	}
}

void display_graph ()
{
	/* set window title */
	char title[256];
	sprintf(title,"History of averages, #points = %05d",Npoints);
	glutSetWindowTitle(title);

	int i;
	GLfloat x, y, xp,yp;

	double offset = (Npoints>HISTORY_LENGTH) ? 
		((double)((Npoints-HISTORY_LENGTH)%HISTORY_LENGTH) / (double)HISTORY_LENGTH) : 0.;
	draw_grid ( offset );
	
	glLineWidth(2.5); 

	for (i=0; ( (i<HISTORY_LENGTH-1) && (i<Npoints) ); i++)
	{
		x  = (GLfloat)((float)i / (float)HISTORY_LENGTH);
		xp = (GLfloat)((float)(i+1) / (float)HISTORY_LENGTH);
		y  = (GLfloat)(history[i].Xavg);
		yp = (GLfloat)(history[i+1].Xavg);
		
		glColor3f((invert_color)?1:0,(invert_color)?1-0.8:0.8,(invert_color)?1:0.0);
		glBegin(GL_LINES);
			glVertex3f (x, y+2., 0.0);
			glVertex3f (xp, yp+2., 0.0);
		glEnd();

		y  = (GLfloat)(history[i].Yavg);
		yp = (GLfloat)(history[i+1].Yavg);
		
		glColor3f((invert_color)?1-0.8:0.8,(invert_color)?1-0.1:0.1,(invert_color)?1:0.0);
		glBegin(GL_LINES);
			glVertex3f (x, y+1., 0.0);
			glVertex3f (xp, yp+1., 0.0);
		glEnd();

		if (mode==MARKOV) {
		y  = (GLfloat)(history[i].acceptance);
		yp = (GLfloat)(history[i+1].acceptance);
		
		glColor3f((invert_color)?1:0.0,(invert_color)?1-0.1:0.1,(invert_color)?1-0.8:0.8);
		glBegin(GL_LINES);
			glVertex3f (x, y+0., 0.0);
			glVertex3f (xp, yp+0., 0.0);
		glEnd();
		}
	}

	glFlush();
}

void update_displays ()
{
	/* call this update every 'refresh_period' microsecs */
	glutTimerFunc( refresh_period, update_displays, 0 );

	if (!paused)
	{
		switch ( mode )
		{
			case NAIVE_UNIFORM: 
				pick_random_point( 0 );
				break;
			case NAIVE_DISTR: 
				pick_random_point( 1 );
				break;
			case MARKOV:
				metropolis_step ( );
				break;
		}
		//paused=1;

		Xavg = Xtot / Norm; 
		Yavg = Ytot / Norm;
		//fprintf(stderr,"Npts = %d\tXavg = %lg\tYavg = %lg\tXtot = %lg\tYtot=%lg\tNorm = %lg\n", Npoints, Xavg, Yavg,Xtot,Ytot,Norm );

		fill_history();

		/* refresh window */
		glutSetWindow (win1);
		glutPostRedisplay ();
		glutSetWindow( win2);
		glutPostRedisplay ();
	}
}

