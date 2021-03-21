/* $Id: paltex.c,v 1.6 2000/10/05 07:17:43 joukj Exp $ */

/*
 * Paletted texture demo.  Written by Brian Paul.
 * This program is in the public domain.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <float.h>
#define GL_GLEXT_LEGACY
#include <GL/glut.h>


static float Rot = 0;
static GLboolean Anim = 1;

static float *T;
static int n;
static float *x;
static float *y;
static float Tmax = FLT_MIN;
static float Tmin = FLT_MAX;

static float z = -1;
static float fx = 0;
static float fy = 0;

static void Idle( void )
{
   float t = glutGet(GLUT_ELAPSED_TIME) * 0.001;  /* in seconds */
   //Rot = t * 360 / 4;  /* 1 rotation per 4 seconds */
   glutPostRedisplay();
}


static void Initialize(char *filename){
  FILE *f;
  f = fopen(filename, "r");

  if (f == NULL)
  {
      printf("Error opening file!\n");
      exit(1);
  }

  int nx, ny;
  fscanf(f,"%d %d", &nx, &ny);
  n = nx;
  int nT = nx*ny;
  
  T = malloc(nT*sizeof(float));

  for(int i = 0; i < nT; i++){
    float t;
    fscanf(f,"%f",&t);
    T[i] = t;

    if (t > Tmax) Tmax = t;
    else if (t < Tmin) Tmin = t;
   
   }

   x = malloc((n)*sizeof(float));
   y = malloc((n)*sizeof(float));
   x[0] = 0;
   y[0] = 0;
   x[n-1] = 1;
   y[n-1] = 1;

   double h = 1.0/((double) n - 1);
   for(int i = 1; i < (n-1); i++){
    x[i] = x[i-1] + h;
    y[i] = y[i-1] + h;
   }

  fclose(f);

  glClearColor (0.5,0,0.5,0);
}

static void setColour(float t, float Tmax, float Tmin, float *min_colour, float *max_colour){
   float dt = Tmax - Tmin;
   float tn = (t-Tmin)/dt;
   float colours[3];
   
   for(int i = 0; i < 3; i++){
      float c = (1-tn)*min_colour[i] + tn*max_colour[i];
      colours[i] = c;
   }
   glColor3f(colours[0],colours[1],colours[2]);
}

static void Display( void )
{

   glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
   /* draw background gradient */
   float min_color[3] = {0,0,0};
   float max_color[3] = {1,1,1};

   
  
   glPushMatrix();
   
   glMatrixMode( GL_PROJECTION );
   glLoadIdentity();
   
   glMatrixMode( GL_MODELVIEW );
   glLoadIdentity();
   glTranslatef(-0.5,-0.5,0);
  
   glBegin(GL_QUADS);
   
    for(int j = 0; j < n-1; j++){
      for(int i = 0; i < n-1; i++){
        float t1 = 2*T[i*n+j];
        float t2 = 2*T[i*n+j+1];
        float t3 = 2*T[(i+1)*n+j+1];
        float t4 = 2*T[(i+1)*n+j];

        
        setColour(t1, Tmax, Tmin, min_color, max_color); glVertex2f(x[i], y[j]);
        
        setColour(t2, Tmax, Tmin, min_color, max_color); glVertex2f(x[i], y[j+1]);
        
        setColour(t3, Tmax, Tmin, min_color, max_color); glVertex2f(x[i+1], y[j+1]);
        
        setColour(t4, Tmax, Tmin, min_color, max_color); glVertex2f(x[i+1], y[j]);
    }
    }

   glEnd();
   

   glPopMatrix();

   glutSwapBuffers();
}


static void Reshape( int width, int height )
{
   glViewport( 0, 0, width, height );
   glMatrixMode( GL_PROJECTION );
   glLoadIdentity();
   glOrtho( -1.5, 1.5, -1.0, 1.0, -1.0, 1.0 );
   glMatrixMode( GL_MODELVIEW );
   glLoadIdentity();
}


static void Key( unsigned char key, int x, int y )
{
   (void) x;
   (void) y;
   switch (key) {
      case 27:
         exit(0);
         break;
      case 'w':
         fy += 1;
         break;
      case 's':
         fy -= 1;
         break;
      case 'a':
         fx += 1;
         break;
      case 'd':
         fx -= 1;
         break;
      case ' ':
         Anim = !Anim;
         if (Anim)
            glutIdleFunc( Idle );
         else
            glutIdleFunc( NULL );
         break;
   }
   glutPostRedisplay();
}







int main( int argc, char *argv[] )
{
   glutInit( &argc, argv );
   glutInitWindowPosition( 0, 0 );
   glutInitWindowSize( 1000,  1000);

   glutInitDisplayMode( GLUT_RGB | GLUT_DOUBLE );

   glutCreateWindow(argv[0]);


   Initialize(argv[1]);
   glEnable(GL_DEPTH_TEST);
  
   glutKeyboardFunc( Key );
   glutDisplayFunc( Display );
   if (Anim)
      glutIdleFunc( Idle );

   glutMainLoop();

   free(T);
   free(x);
   free(y);

   return 0;
}
