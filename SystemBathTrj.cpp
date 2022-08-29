// SystemBathTrj.cpp : Este archivo contiene la función "main". La ejecución del programa comienza y termina ahí.
//

#include <iostream>			// std::cout, std::fixed
#include <sstream>			// std::stringstream
#include <iomanip>			// std::setprecision
#include <fstream>			// std::ofstream

#include "utilities.h"
#include "Oscillator.h"
#include "Dynamics.h"

utilities ut;

void programUsage()
{
	std::cout << "Program Usage" << std::endl;
	std::cout << ".exe [bathMass] [TrjTime] [Trj Number] [Direction {-1,1}] [Initial Point 1 (0 for saddle point)] [Initial Point 2 (0 for saddle point)]" << std::endl;
}

/* Potential and Gradient of the system */
double DOS_V(std::vector<double> coeff, std::vector<double> r)
{
	double ep = 0.;
	for (int i = 0; i < 5; ++i)
	{
		ep += coeff[i] * ut.powerd(r[0], i);
	}
	ep += coeff[5] * ut.powerd(r[1] - coeff[6], 2);
	ep += coeff[7] / ut.powerd(r[1] - r[0], 12);

	return ep;
}
std::vector<double> DOS_G(std::vector<double> coeff, std::vector<double> r)
{
	std::vector<double> grad(2, 0.);
	grad[0] = 0;
	for (int j = 1; j < 5; ++j)
	{
		grad[0] += double(j) * coeff[j] * ut.powerd(r[0], j - 1);
	}
	grad[0] += 12. * coeff[7] / ut.powerd(r[1] - r[0], 13);
	grad[1] = 0;
	grad[1] = 2. * coeff[5] * (r[1] - coeff[6]) - 12. * coeff[7] / ut.powerd(r[1] - r[0], 13);
	return grad;
}

int main(int argc, char* argv[])
{
		/* Check the correct Parameters */
#pragma region Program Usage
	if ((argc != 7))
	{
		programUsage();
		return 1;
	}
	
#pragma endregion


	/* Variables introduced by the User */
#pragma region Input variables
	double bathMass = strtof(argv[1], NULL);					// Mass of the Bath
	double TrjTime = strtof(argv[2], NULL);						// Time of the Trj
	int I = std::stoi(argv[3], NULL);							// Number of the trajectory
	int direction = std::stoi(argv[4], NULL);					// Direction of the time Step
	
	std::vector<double> R0;										// Initial Position
	if (std::stoi(argv[5], NULL) == 0)
	{
		//std::cout << "Saddle point" << std::endl;
		R0 = { 1.38795 ,2.26838 };
	}
	else
	{
		R0 = { strtof(argv[5], NULL) ,strtof(argv[6], NULL) };
	}

#pragma endregion

/* Set the Initial Conditions */
#pragma region System Initial Conditions

	/* Initiate the Oscillator */
	Oscillator Oscf({ 321.904484,-995.713452,1118.689573,-537.856726,92.976121,1.0,1.0,0.01 },	/* Coefficients from the oscillator */\
	{ 1., bathMass },																/* Mass Values */\
		DOS_V, DOS_G);																/* Potential and Gradient Functions */
	Oscillator Oscb({ 321.904484,-995.713452,1118.689573,-537.856726,92.976121,1.0,1.0,0.01 },	/* Coefficients from the oscillator */\
	{ 1., bathMass },																/* Mass Values */\
		DOS_V, DOS_G);																/* Potential and Gradient Functions */

	srand(time(NULL));												// Initialize the random seed


	double Energy = 3.691966889;									// Energy of the system
	double timeStep = 1.e-3;										// Set the Time Step/Precision of the Dynamic multiply by the direction

#pragma endregion

/* Initialize the Variables */
#pragma region Process Required Variables
	/* We are going to perform two Dynamics at the same time,
	one in Foward and the other in Backward time.
	To do so, we need to initialize 2 Oscillators and 2 Dynamics*/

	/* Initialize the Oscillators */
#pragma region Oscillator Variables
	Oscf.setInitP(R0, true);										// Initiate at the shadle point
	Oscf.setInitialVel_NVE(Energy, I, 100);							// Set initial Random Velocities to keep energy
#pragma endregion

	/* Initialize the dynamics */
#pragma region Dynamics Variables
	Dynamics Dynf;													// Set the forward Dynamic							
	Dynf.setTimeStep(timeStep * direction);							// Set the Time step
	Dynf.setTime(TrjTime);											// and the total time (nsteps = totalTime/timeTtep)

	/* Variables fot the Stochastic Dynamics */
#if KEY_DETERM==0													
	double beta = 4.;												// Effective temperature
	double gamma = 4.;												// Disipation factor
	Dynf.setLangevin(Oscf._redMass.back(), beta, gamma);

	long T = time(0);														// Seed for the random 
	std::vector<double> gasdev = ut.box_muller(Dynf._numberStep, &T);		// Gas deviation
	double rtherm;
#endif

#pragma endregion

#pragma endregion


	std::ofstream outputFile;										//Saving file
	std::string filename = "trj_M" + std::to_string(int(bathMass)) + "T" + std::to_string(int(TrjTime)) + "_" + std::to_string(direction) + "_n" + std::to_string(I) + "_0-0.txt";
	outputFile.open(filename, std::ios::out | std::ios::trunc);
	/* Save the Initial properties of the trajectory */
	outputFile << Oscf._position[0] << ' ' << Oscf._position[1] << ' ' << ' ' << Oscf.calcMomenta()[0] << ' ' << Oscf.calcMomenta()[1] << std::endl;

	/* Initate the trayectory */
	for (size_t j = 0; j < Dynf._numberStep; j++)
	{
		/* Save Osc Position */
		double position0 = Oscf._position[0];
		double position1 = Oscf._position[1];
		double momenta0 = Oscf._velocity[0];
		double momenta1 = Oscf._velocity[1];

		/* Perform a Dinamic Step */
#pragma region Dynamic Step
#if KEY_DETERM==1
		Dynf.DynamicStep(Oscf);										//Dynamic forward step
#else
		rtherm = Dynf.sigma * gasdev[j];								//Bath Random
		Dynf.DynamicStep(Oscf, rtherm);								//Dynamic foward step
#endif
#pragma endregion

		/* Save the properties of the trajectory */
		outputFile << Oscf._position[0] << ' ' << Oscf._position[1] << ' ' << ' ' << Oscf.calcMomenta()[0] << ' ' << Oscf.calcMomenta()[1] << std::endl;

#if KEY_BOX==1
		if (Oscf._position[0] > 1.5 || Oscf._position[0] < 1.3)
		{
			outputFile << std::endl;
			
			Oscf._position[0] = position0;
			Oscf._position[1] = position1;
			Oscf._velocity[0] = momenta0;
			Oscf._velocity[1] = momenta1;

			Oscf._velocity[0] = Oscf._velocity[0] * -1.;

			outputFile << Oscf._position[0] << ' ' << Oscf._position[1] << ' ' << ' ' << Oscf.calcMomenta()[0] << ' ' << Oscf.calcMomenta()[1] << std::endl;
		}

#endif // KEY_BOX==1

	}
}
