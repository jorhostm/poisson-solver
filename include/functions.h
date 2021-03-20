#ifndef FUNCTIONS_H
#define FUNCTIONS_H

static double zero(double x, double y){

	return 0;
	
}

static double mix(double x, double y){
    if(y == 1.0){
        return 1.0;
    }

    return 0;
}

static double phi(double x, double y){

	return 0.25*(x*x+y*y);
	
}

static double g1(double x, double y){

	return 12-12*x-12*y;

}
static double g2(double x, double y){

	return (6-12*x)*(3*y*y-2*y*y*y) + (3*x*x-2*x*x*x)*(6-12*y);

}
#endif