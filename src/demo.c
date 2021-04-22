
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <float.h>
#define GL_GLEXT_LEGACY
#include <GL/gl.h>
#include <GL/freeglut.h>

static float *T;
static int n;
static float *x;
static float *y;
static float Tmax = FLT_MIN;
static float Tmin = FLT_MAX;


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

   fclose(f);

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


  glClearColor (0.9,0.9,0.9, 1);
}

static void setColour(float t){
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

static void Display( void )
{

   glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
   glPushMatrix();
   
   glMatrixMode( GL_PROJECTION );
   glLoadIdentity();
   
   glMatrixMode( GL_MODELVIEW );
   glLoadIdentity();
   glTranslatef(-0.5,-0.5,0);
  
   glBegin(GL_QUADS);
   
    for(int j = 0; j < n-1; j++){
      for(int i = 0; i < n-1; i++){
        float t1 = T[i*n+j];
        float t2 = T[i*n+j+1];
        float t3 = T[(i+1)*n+j+1];
        float t4 = T[(i+1)*n+j];

        
        setColour(t1); glVertex2f(x[i], y[j]);
        
        setColour(t2); glVertex2f(x[i], y[j+1]);
        
        setColour(t3); glVertex2f(x[i+1], y[j+1]);
        
        setColour(t4); glVertex2f(x[i+1], y[j]);
    }
    }

   glEnd();
   

   glPopMatrix();

   glutSwapBuffers();
}


static void Reshape( int w, int h )
{
   glViewport(0, 0, (GLsizei) w, (GLsizei) h);
   glMatrixMode(GL_PROJECTION);
   glLoadIdentity();
   
   glMatrixMode(GL_MODELVIEW);
   glLoadIdentity();
   glTranslatef (-0.5f, -0.5f, -5.0f);
}


static void Key( unsigned char key, int x, int y )
{
   switch (key) {
      case 27:
         exit(0);
         break; 
   }
   
}







int main( int argc, char *argv[] )
{
   glutInit( &argc, argv );
   glutInitWindowPosition( 0, 0 );
   glutInitWindowSize( 1000,  1000);

   glutInitDisplayMode( GLUT_RGB | GLUT_DOUBLE );

   glutCreateWindow(argv[0]);
   glutReshapeFunc(Reshape);

   Initialize(argv[1]);
   glEnable(GL_DEPTH_TEST);
  
   glutKeyboardFunc( Key );
   glutDisplayFunc( Display );

   
   glutMainLoop();

   free(T);
   free(x);
   free(y);

   return 0;
}
