
/*******************************************************************************
*
* File rng_gauss.c
*
* Copyright (C) 2016 Sebastian Schmalzbauer
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Random number generator "rng_gauss"
*
* Use the routines from "ranlxd" of Martin Luescher to initialize, set the seed, etc.
*
* The externally accessible functions are
*
*   double rng_gauss_single(double sigma)
*     Returns a single Gaussian random number with standard deviation sigma
*
*   void rng_gauss_array(double g[], int n, double sigma)
*     Fills the first n elements of the array g with Gaussian random numbers with standard deviation sigma
*
*******************************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include "ranlxd.h"
#include <math.h>

double uniform_lookup[2];
int gaussian_state = 0;
double gaussian_lookup[2];
double r, phi;

double rng_gauss_single(double sigma)
{
   if(gaussian_state == 0)
   {
      ranlxd(uniform_lookup, 2);
      r = sqrt(-2.0*pow(sigma,2)*log(uniform_lookup[0]));
      phi = 2.0*M_PI*uniform_lookup[1];
      gaussian_lookup[0] = r*cos(phi);
      gaussian_lookup[1] = r*sin(phi);
      gaussian_state = 1;
      return gaussian_lookup[0];
   }
   else
   {
      gaussian_state = 0;
      return gaussian_lookup[1];
   }
}

int gaussian_state2 = 0;
double gaussian_lookup2[2];

double rng_gauss_single2(double sigma)
{
   if(gaussian_state2 == 0)
   {
      ranlxd(uniform_lookup, 2);
      r = sqrt(-2.0*pow(sigma,2)*log(uniform_lookup[0]));
      phi = 2.0*M_PI*uniform_lookup[1];
      gaussian_lookup2[0] = r*cos(phi);
      gaussian_lookup2[1] = r*sin(phi);
      gaussian_state2 = 1;
      return gaussian_lookup2[0];
   }
   else
   {
      gaussian_state2 = 0;
      return gaussian_lookup2[1];
   }
}


void rng_gauss_array(double g[], int n, double sigma)
{
   int i;   
   for(i = 0; i < n;)
   {
      ranlxd(uniform_lookup, 2);
      r = sqrt(-2.0*pow(sigma,2)*log(uniform_lookup[0]));
      phi = 2.0*M_PI*uniform_lookup[1];
      g[i++] = r*cos(phi);
      if (i < n)
         g[i++] = r*sin(phi);
   }
}

double ran3()
{
	double r;
	ranlxd( &r, 1);
	return (r);
}

