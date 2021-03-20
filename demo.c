/* $Id: paltex.c,v 1.6 2000/10/05 07:17:43 joukj Exp $ */

/*
 * Paletted texture demo.  Written by Brian Paul.
 * This program is in the public domain.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#define GL_GLEXT_LEGACY
#include <GL/glut.h>


static float Rot = 0;
static GLboolean Anim = 1;

static float *T;
static int n;
static float *x;
static float *y;
static float Tmax;
static float Tmin;

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
    fscanf(f,"%f",T+i);
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

  
}


static void Display( void )
{
   /* draw background gradient */


   
   glBegin(GL_POLYGON);
   glColor3f(1.0, 1.0, 1.2); glVertex2f(-1.5, -1.0);
   glColor3f(1.0, 1.0, 1.2); glVertex2f( 1.5, -1.0);
   glColor3f(1.0, 1.0, 1.0); glVertex2f( 1.5,  1.0);
   glColor3f(1.0, 1.0, 1.0); glVertex2f(-1.5,  1.0);
   glEnd();

   glPushMatrix();
   glRotatef(Rot, 0, 0, 1);

   
   glBegin(GL_QUADS);
   
    for(int j = 0; j < n-1; j++){
      for(int i = 0; i < n-1; i++){
        float t1 = 2*T[i*n+j];
        float t2 = 2*T[i*n+j+1];
        float t3 = 2*T[(i+1)*n+j+1];
        float t4 = 2*T[(i+1)*n+j];
        glColor3f(t1, t1, t1); glVertex2f(0.7*x[i], 0.7*y[j]);
        glColor3f(t2, t2, t2); glVertex2f(0.7*x[i], 0.7*y[j+1]);
        glColor3f(t3, t3, t3); glVertex2f(0.7*x[i+1], 0.7*y[j+1]);
        glColor3f(t4, t4, t4); glVertex2f(0.7*x[i+1], 0.7*y[j]);
    }
    }


   /*
   glTexCoord2f(0, 1);  glVertex2f(-1, -0.5);
   glTexCoord2f(1, 1);  glVertex2f( 1, -0.5);
   glTexCoord2f(1, 0);  glVertex2f( 1,  0.5);
   glTexCoord2f(0, 0);  glVertex2f(-1,  0.5);
   */
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
      case 's':
         Rot += 0.5;
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
