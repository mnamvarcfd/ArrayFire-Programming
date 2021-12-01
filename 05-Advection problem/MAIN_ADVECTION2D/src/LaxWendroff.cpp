#include "LaxWendroff.h"

LaxWendroff::LaxWendroff()
{
}

LaxWendroff::~LaxWendroff()
{
}

LaxWendroff::LaxWendroff(Field &field_, double a_, double b_):Solver<D2Q9Stencil>(field_)
{

	stencil = D2Q9Stencil(field_);
	
	a = a_;
	b = b_;
}


void LaxWendroff::setCoeff()
{
	double a2 = a * a;
	double b2 = b * b;
	double dx2 = dx * dx;
	double dy2 = dy * dy;
	double dt2 = dt * dt;

	coeff2D[1][1] = 1.0 - dt2 * (a2/dx2 + b2/dy2);
	coeff2D[1][2] = 0.5 * (b2 * dt2/dy2 - b * dt/dy);
	coeff2D[2][2] = 0.25 * a * b * dt2 / (dx*dy);
	coeff2D[2][1] = 0.5 * (a2 * dt2/dx2 - a * dt / dx);
	coeff2D[2][0] = -0.25 * a * b * dt2 / (dx * dy);
	coeff2D[1][0] = 0.5 * (b2 * dt2 / dy2 + b * dt / dy);
	coeff2D[0][0] = 0.25 * a * b * dt2 / (dx * dy);
	coeff2D[0][1] = 0.5 * (a2 * dt2 / dx2 + a * dt / dx);
	coeff2D[0][2] = -0.25 * a * b * dt2 / (dx * dy);


	stencil.set_coeff2D(coeff2D);

}


void LaxWendroff::timeStep() {

	double dx = field.get_dx();
	double dy = field.get_dy();

	dt = (CFL * min(dx , dy)) / sqrt(2 * (a * a + b * b));

}


