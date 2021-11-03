#include "AnalyticalSolution.h"

AnalyticalSolution::AnalyticalSolution()
{
}

AnalyticalSolution::~AnalyticalSolution()
{
}

AnalyticalSolution::AnalyticalSolution(int nx, int ny, double xMin, double xMax, double yMin, double yMax) : Field(nx, ny, xMin, xMax, yMin, yMax)
{
	for (int j = 0; j < ny; j++) {
		double y = j * dy;
		for (int i = 0; i < nx; i++) {
			double x = i * dx;

			int iNode = getGlobInx(i, j);
			var[iNode] = sin(x) + cos(y);

		}
	}

}

