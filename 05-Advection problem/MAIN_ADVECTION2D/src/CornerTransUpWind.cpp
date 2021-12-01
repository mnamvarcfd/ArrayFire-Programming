#include "CornerTransUpWind.h"

CornerTransUpWind::CornerTransUpWind()
{
}

CornerTransUpWind::~CornerTransUpWind()
{
}

CornerTransUpWind::CornerTransUpWind(Field &field_, double a_, double b_):Solver<D2Q9Stencil>(field_)
{

	stencil = D2Q9Stencil(field_);
	
	a = a_;
	b = b_;
}


void CornerTransUpWind::setCoeff()
{

	double dxdy = dx*dy;
	double dt2dxdy = dt*dt / (dx * dy);
	double ap = max(a , 0.0);
	double am = min(a , 0.0);
	double bp = max(b , 0.0);
	double bm = min(b , 0.0);

	coeff2D[1][1] = 1.0 - dt * (ap - am) / dx - dt * (bp - bm) / dy + dt2dxdy*(ap*bp-ap*bm-am*bp+am*bm);
	coeff2D[1][2] = -dt * bm / dy + dt2dxdy * (ap * bm - am * bm);
	coeff2D[2][2] = dt2dxdy * am * bm;
	coeff2D[2][1] = -dt * am / dx + dt2dxdy * (am * bp - am * bm);
	coeff2D[2][0] = -dt2dxdy * am * bp;
	coeff2D[1][0] = dt * bp / dy + dt2dxdy * (am * bp - ap * bp);
	coeff2D[0][0] = dt2dxdy * ap * bp;
	coeff2D[0][1] = dt * ap / dy + dt2dxdy * (ap * bm - ap * bp);
	coeff2D[0][2] = -dt2dxdy * ap * bm;

	stencil.set_coeff2D(coeff2D);

}


void CornerTransUpWind::timeStep() {

	dt = CFL / max(abs(a/dx) , abs(b/dy));

}