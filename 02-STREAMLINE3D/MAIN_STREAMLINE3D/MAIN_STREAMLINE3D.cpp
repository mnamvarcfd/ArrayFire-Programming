//////////////////////////////////////////////////////////////
// DATE: 2020 - 07 - 31
// Code written by Sébastien Leclaire(sebastien.leclaire@polymtl.ca)
// This code solves the 3D advection of massless particles along 
// the streamlines of a chosen analytical velocity field.
//////////////////////////////////////////////////////////////
//https://atozmath.com/example/CONM/RungeKutta.aspx?he=e&q=rk4
//https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods
//https://en.wikipedia.org/wiki/Linear_multistep_method#Adams%E2%80%93Bashforth_methods


#include "arrayfire.h"
#include "IO/IO.h"

#undef min
#undef max

typedef double T;

#include <iostream>
#include <iomanip> // std::setprecision
# define M_PI 3.14159265358979323846  /* pi */

const T EPS = std::numeric_limits<T>::epsilon();


af::array pos_forwardEuler(const af::array& nPos, T dt, const T tPos, const T omega)
{
	return nPos + dt * nPos * sin(omega * tPos) * exp(-tPos);
}

af::array pos_SecondOrderRungeKutta(const af::array& nPos, T dt, const T tPos, const T omega)
{
	////T t = tPos;
	////af::array k1 = dt * nPos * sin(omega * t) * exp(-t);

	////t = tPos + dt * 0.5;
	////af::array k2 = dt * (nPos + 0.5 * k1) * sin(omega * t) * exp(-t);

	////return nPos + k2;


	T t = tPos;
	af::array k1 = nPos * sin(omega * t) * exp(-t);

	t = tPos + dt;
	af::array k2 = (nPos + dt * k1) * sin(omega * t) * exp(-t);

	return nPos + dt * 0.5 * (k1 + k2);
}

af::array pos_FourthOrderRungeKutta(const af::array& nPos, T dt, const T tPos, const T omega)
{
	T t = tPos;
	af::array k1 = dt * nPos * sin(omega * t) * exp(-t);

	t = tPos + dt * 0.5;
	af::array k2 = dt * (nPos + 0.5 * k1) * sin(omega * t) * exp(-t);

	t = tPos + dt * 0.5;
	af::array k3 = dt * (nPos + 0.5 * k2) * sin(omega * t) * exp(-t);

	t = tPos + dt;
	af::array k4 = dt * (nPos + k3) * sin(omega * t) * exp(-t);

	return nPos + k1 / 6.0 + k2 / 3.0 + k3 / 3.0 + k4 / 6.0;
}

af::array pos_AdamsBashforth(const af::array& nm1Pos, const af::array& nPos, T dt, const T tPos, const T omega)
{
	T t = tPos - dt;
	af::array Xnp1 = nPos
		+ 1.5 * dt * (nPos * sin(omega * tPos) * exp(-tPos))
		- 0.5 * dt * (nm1Pos * sin(omega * t) * exp(-t));

	return Xnp1;
}

af::array init_AdamsBashforth(const af::array& nPos, T dt, const T tPos, const T omega)
{
	T delT = 0.5 * dt;
	af::array xPredict = pos_forwardEuler(nPos, delT, tPos, omega);

	T t = tPos + 0.5 * dt;
	af::array Xnp1 = xPredict
		+ 1.5 * delT * (xPredict * sin(omega * t) * exp(-t))
		- 0.5 * delT * (nPos * sin(omega * tPos) * exp(-tPos));

	return Xnp1;
}

af::array pos_analytical(const af::array& posInit,  const T tPos, const T omega)
{
	return posInit * std::exp(-(std::exp(-tPos) * (std::sin(omega * tPos) + omega * cos(omega * tPos))) / (omega * omega + 1))*std::exp(omega / (omega * omega + (T)1));
}

af::array pos_analytical_steadyState(const af::array& posInit, const T omega)
{
	return posInit * std::exp(omega / (omega * omega + (T)1));
}

int main(int argc, char *argv[])
{
	af::setBackend(AF_BACKEND_CPU);
	af::setDevice(0);
	af::info(); std::cout << std::endl;

	int nParticle = 100000; // This line can be change only for the performance analysis, but should not be change for the accuracy analysis (streamline3D.xlsx)

	/////////////// Do not change code between thoses lines below ///////////////

	T omegaX = (T)1 * M_PI / (T)3;
	T omegaY = (T)2 * M_PI / (T)3;
	T omegaZ = (T)3 * M_PI / (T)3;
	
	// Initial position of all particles
	af::setSeed(0);
	af::array xInit = af::randu(nParticle, 1, type::TYPE_AF<T>());
	af::array yInit = af::randu(nParticle, 1, type::TYPE_AF<T>());
	af::array zInit = af::randu(nParticle, 1, type::TYPE_AF<T>());

	/////////////// Do not change code between thoses lines above ///////////////

	af::array xPosNm1;
	af::array yPosNm1;
	af::array zPosNm1;

	// Select the constant time step.
	T dt = 1./5;

	for (int times = 0; times < 5; times++) {
		dt *= 0.5;

		// Main time marching loop
		bool isSteadyStateReached = false;
		int iT = 0;
		T tPos = 0;
		af::array xPosLast = xInit;
		af::array yPosLast = yInit;
		af::array zPosLast = zInit;
		af::array xPosNew, yPosNew, zPosNew;

		while (isSteadyStateReached == false)
		{
			// Analytical solver for the position of the particles (at tPos)
			//xPosNew = pos_analytical(xInit, tPos, omegaX);
			//yPosNew = pos_analytical(yInit, tPos, omegaY);
			//zPosNew = pos_analytical(zInit, tPos, omegaZ);

			// Forward Euler scheme
			//////xPosNew = pos_forwardEuler(xPosLast, dt, tPos, omegaX);
			//////yPosNew = pos_forwardEuler(yPosLast, dt, tPos, omegaY);
			//////zPosNew = pos_forwardEuler(zPosLast, dt, tPos, omegaZ);

			// Second Order Runge-Kutta scheme (RK2 aka midpoint) 
			//////xPosNew = pos_SecondOrderRungeKutta(xPosLast, dt, tPos, omegaX);
			//////yPosNew = pos_SecondOrderRungeKutta(yPosLast, dt, tPos, omegaY);
			//////zPosNew = pos_SecondOrderRungeKutta(zPosLast, dt, tPos, omegaZ);

			// Fourth Order Runge-Kutta scheme (RK4)
			//xPosNew = pos_FourthOrderRungeKutta(xPosLast, dt, tPos, omegaX);
			//yPosNew = pos_FourthOrderRungeKutta(yPosLast, dt, tPos, omegaY);
			//zPosNew = pos_FourthOrderRungeKutta(zPosLast, dt, tPos, omegaZ);

			//Two - step Adams - Bashforth scheme with an initial half time step and a forward Euler’s predictor.
			if (iT == 0) {
				xPosNew = init_AdamsBashforth(xPosLast, dt, tPos, omegaX);
				yPosNew = init_AdamsBashforth(yPosLast, dt, tPos, omegaY);
				zPosNew = init_AdamsBashforth(zPosLast, dt, tPos, omegaZ);

				xPosNm1 = xPosNew;
				yPosNm1 = yPosNew;
				zPosNm1 = zPosNew;
			}
			else
			{
				xPosNew = pos_AdamsBashforth(xPosNm1, xPosLast, dt, tPos, omegaX);
				yPosNew = pos_AdamsBashforth(yPosNm1, yPosLast, dt, tPos, omegaY);
				zPosNew = pos_AdamsBashforth(zPosNm1, zPosLast, dt, tPos, omegaZ);

				xPosNm1 = xPosLast;
				yPosNm1 = yPosLast;
				zPosNm1 = zPosLast;
			}


			// Next time step
			tPos += dt;
			iT += 1;

			// Check if steady state is reached every 50 time steps.
			if (iT % 50 == 0 && iT > 0)
			{
				// Calculate the difference L1 between the current position and the last position of all particles.
				T L1sum = af::sum<T>(af::abs(xPosNew - xPosLast) + af::abs(yPosNew - yPosLast) + af::abs(zPosNew - zPosLast));
				//std::cout << std::setprecision(16) << "L1sum: " << L1sum << std::endl;

				// Activate stopping criterion if steady state is reached
				if (L1sum < (T)100 * EPS) isSteadyStateReached = true;
			}

			// Eval ArrayFire JIT tree for the main variables because a recursion loop may cause OOM or performance issue.
			xPosNew.eval();
			yPosNew.eval();
			zPosNew.eval();

			// Also advance in time the last position of the particles
			xPosLast = xPosNew;
			yPosLast = yPosNew;
			zPosLast = zPosNew;

		}

		// Analytical solution at steady state (i.e. at t_\inf)
		af::array xFinal = pos_analytical_steadyState(xInit, omegaX);
		af::array yFinal = pos_analytical_steadyState(yInit, omegaY);
		af::array zFinal = pos_analytical_steadyState(zInit, omegaZ);

		// Error at steady state between computed and analytical positions
		T L2error = std::sqrt(af::sum<T>((xPosLast - xFinal) * (xPosLast - xFinal) + (yPosLast - yFinal) * (yPosLast - yFinal) + (zPosLast - zFinal) * (zPosLast - zFinal)));
		std::cout << std::setprecision(16) << "L2error: " << L2error << std::endl;

	}


	return 0;
}
