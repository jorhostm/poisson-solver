/**
 * @file demo.c
 * @author Jørgen Høstmark (jorgenhostm@gmail.com)
 * @brief This file is used to plot the solutions from the poisson-solver API, using a colour plot.
 * @version 1.0
 * @date 2021-05-03
 * 
 * 
 */

#include <stdio.h>
#include <stdlib.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <string.h>
#include <float.h>
#define GL_GLEXT_LEGACY
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/freeglut.h>

static float *T;
static int n;
static float *x;
static float *y;
static float Tmax = FLT_MIN;
static float Tmin = FLT_MAX;

static float angleX = 0.0;
static float angleY = 0.0;
static float angleZ = 0.0;

static float dist = 2;

static unsigned int use_2d = 0;

/**
 * @brief Load the solution
 * 
 * @param filename The file to import from
 */
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
    
  // Scan the n^2 values and store them into T.  
  for(int i = 0; i < nT; i++){
    float t;
    fscanf(f,"%f",&t);
    T[i] = t;

   // Update the maximum and minimum values
    if (t > Tmax) Tmax = t;
    else if (t < Tmin) Tmin = t;
   
   }

   fclose(f);

   // Initialise the x- and y-values.
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

  glClearColor (1,1,1,1); // Set the background to white
}

/**
 * @brief Compute and set the colour of a point
 * 
 * @param t The value to compute the colour for. Tmin <= t <= Tmax
 */
static void setColour(float t){
   
   // Normalise the value t
   float dt = Tmax - Tmin;
   float tn = (t-Tmin)/dt;
   
   float r,g,b;
   
   // Set red value
   r = sqrt(tn);

   // Set green value
   g = powf(tn,3);
   
   // Set blue value
   b = sinf(2*M_PI*tn);
   
   
   glColor3f(r, g, b);
}

// Draw the solution in 2D
static void Display2D(){

   glTranslatef(-0.5,-0.5,0);
   glBegin(GL_QUADS);
   
    for(int j = 0; j < n-1; j++){
      for(int i = 0; i < n-1; i++){
        float t1 = T[i*n+j];
        float t2 = T[i*n+j+1];
        float t3 = T[(i+1)*n+j+1];
        float t4 = T[(i+1)*n+j];

        setColour(t4); glVertex2f(x[i+1], y[j]);
        setColour(t3); glVertex2f(x[i+1], y[j+1]);
        setColour(t2); glVertex2f(x[i], y[j+1]);
        setColour(t1); glVertex2f(x[i], y[j]);
    }
    }
    glEnd();

}

// Draw the solution in 3D
static void Display3D(){
   
   glRotatef(angleX,1,0,0);
   glRotatef(angleY,0,1,0);
   glRotatef(angleZ,0,0,1);
   glTranslatef(-0.5,-0.5,-Tmax/2);
   glBegin(GL_QUADS);
   
    for(int j = 0; j < n-1; j++){
      for(int i = 0; i < n-1; i++){
        float t1 = T[i*n+j];
        float t2 = T[i*n+j+1];
        float t3 = T[(i+1)*n+j+1];
        float t4 = T[(i+1)*n+j];

        setColour(t4); glVertex3f(x[i+1], y[j], t4);
        setColour(t3); glVertex3f(x[i+1], y[j+1], t3);
        setColour(t2); glVertex3f(x[i], y[j+1], t2);
        setColour(t1); glVertex3f(x[i], y[j], t1);
    }
    }

   glEnd();
   
   // Draw box around tha solution
   glBegin(GL_LINE_STRIP);
   glColor3f(0,0,0);
   glVertex3f(0,0,0);
   glVertex3f(1,0,0);
   glVertex3f(1,0,Tmax);
   glVertex3f(1,1,Tmax);
   glVertex3f(1,1,0);
   glVertex3f(0,1,0);
   glVertex3f(0,1,Tmax);
   glVertex3f(0,0,Tmax);
   glVertex3f(0,0,0);
   glVertex3f(0,1,0);
   glVertex3f(0,1,Tmax);
   glVertex3f(1,1,Tmax);
   glVertex3f(1,1,0);
   glVertex3f(1,0,0);
   glVertex3f(1,0,Tmax);
   glVertex3f(0,0,Tmax);
   glEnd();
}

/**
 * @brief The display function
 * 
 */
static void Display( void )
{

   glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
   
   glPushMatrix();

   if (use_2d)
   {
      Display2D();
   }
   else
   {
      Display3D();
   }

   glPopMatrix();
   glutSwapBuffers();
}


static void Reshape( int w, int h )
{
   glViewport (0, 0, (GLsizei) w, (GLsizei) h); 
   glMatrixMode (GL_PROJECTION);
   glLoadIdentity ();
   gluPerspective(60.0, (GLfloat) w/(GLfloat) h, 0.1, 20.0);
   glMatrixMode(GL_MODELVIEW);
   glLoadIdentity();
   gluLookAt (0.0, 0.0, dist, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0);
}


static void Key( unsigned char key, int x, int y )
{
   
   switch (key) {
      case 27:
         exit(0);
         break; 
      
      case 'w':
         angleX += 10;
         break; 
      case 'W':
         dist -= 0.2;
         Reshape(glutGet(GLUT_WINDOW_WIDTH), glutGet(GLUT_WINDOW_HEIGHT));
         break; 
      case 's':
         angleX -= 10;
         break;
      case 'S':
         dist += 0.2;
         Reshape(glutGet(GLUT_WINDOW_WIDTH), glutGet(GLUT_WINDOW_HEIGHT));
         break;
      case 'a':
         angleZ += 10;
         break; 
      case 'd':
         angleZ -= 10;
         break;
      case ' ':
         use_2d = !use_2d;
         angleX = 0;
         angleY = 0;
         angleZ = 0;
         break;
   }
   glutPostRedisplay();
}

int main( int argc, char *argv[] )
{
   glutInit( &argc, argv );
   glutInitWindowPosition( 0, 0 );
   glutInitWindowSize( 800,  800);

   glutInitDisplayMode( GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
   

   glutCreateWindow(argv[0]);
   glEnable(GL_DEPTH_TEST);

   glutReshapeFunc(Reshape);

   Initialize(argv[1]);
  
   glutKeyboardFunc( Key );
   glutDisplayFunc( Display );

   glutMainLoop();

   free(T);
   free(x);
   free(y);

   return 0;
}
