/* =====================================================================================
New Sundials version
===================================================================================== */
// #include "TChem_CommandLineParser.hpp"
#include <sundials/sundials_config.h>
#include <sundials/sundials_nvector.h>
#include <sundials/sundials_math.h>
#include <nvector/nvector_serial.h>

//CVODE includes
#include <cvode/cvode.h>               /* prototypes for CVODE fcts., consts.  */
#include <nvector/nvector_serial.h>    /* access to serial N_Vector            */
#include <sunmatrix/sunmatrix_dense.h> /* access to dense SUNMatrix            */
#include <sunlinsol/sunlinsol_dense.h> /* access to dense SUNLinearSolver      */
#include <sundials/sundials_types.h>   /* defs. of realtype, sunindextype      */
#include <sunlinsol/sunlinsol_spgmr.h>  /* SPGMR SUNLinearSolver                */
#include <sunlinsol/sunlinsol_spbcgs.h>       /* access to SPBCGS SUNLinearSolver            */
#include <sunnonlinsol/sunnonlinsol_newton.h> /* access to Newton SUNNonlinearSolver         */

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <chrono>
#include "Epi3V.h"
#include "ProblemClass.h"

using namespace std;
#define NEQ 3
#define EPS 1e-2

N_Vector InitialConditions();
int RHS(realtype, N_Vector, N_Vector, void *);
int Jtv(N_Vector, N_Vector, realtype, N_Vector, N_Vector , void *, N_Vector);
int ComputeJac(N_Vector, void*);
int CheckStep(realtype, realtype);
int CVodeComputeJacWrapper(realtype t, N_Vector u, N_Vector fy, SUNMatrix Jac,
               void * pb, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);


//using value_type = realtype;

#define BAR "===================="

using namespace std;


int main(int argc, char* argv[])
{
  sundials::Context sunctx;   // SUNDIALS context

  static const realtype FinalTime = 1e1;
	static const int NumBands = 3;
	string MyFile = "OutPutFile.txt";
	realtype StepSize=1e-2;
	int ProgressDots=0; //From 0 to 100; only care about percentage
	long int StepCount = 0;
	realtype PercentDone=0;
	void *userData=nullptr;
	realtype TNow=0;
	realtype TNext=0;
	N_Vector y = N_VNew_Serial(NEQ, sunctx); //The data y
	N_VScale(1.0, InitialConditions(), y);
	realtype *data = NV_DATA_S(y);
	string opt = "EPI3V" ;
	bool LastStepAdjusted = 0;
	realtype KTime = 0;
	realtype PrintStepSize = 0;
	realtype KrylovTol = 1e-14;
  	realtype absTol    = 1e-10;
  	realtype relTol    = 1e-6;
	//Input Parsing
	if(argc == 1)
	{
		//CheckStep(FinalTime, StepSize);
		StepSize=stof(argv[1]);
	}
	else if(argc== 2)
	{
		StepSize=stof(argv[1]);
		std :: cout << "Step-size: " << StepSize << endl;
		//CheckStep(FinalTime, StepSize);
		opt = "EPI2";
	}
	else if(argc == 3)
	{
		StepSize=stof(argv[1]);
		std :: cout << "Step-size: " << StepSize << endl;
		opt = argv[2];
		cout << "Method: " << opt << endl;
	}
	else if(argc == 4 )
	{
		StepSize = stof(argv[1]);
		std :: cout << "Step-size: " << StepSize << endl;
		opt = argv[2];
		cout << "Method: " << opt << endl;
		MyFile = argv[3];
		cout << "Output file: " << MyFile << endl;
	}
	else if(argc == 5)
	{
		StepSize = stof(argv[1]);
		opt = argv[2];
		MyFile = argv[3];
		absTol = stof(argv[4]);
	}
	else if(argc == 6)
	{
		StepSize = stof(argv[1]);
		opt = argv[2];
		MyFile = argv[3];
		absTol = stof(argv[4]);
		relTol = stof(argv[5]);
	}
	else
	{
		cout << "!!!Failure:  Incorrect number of args!!!" << endl;
		exit(EXIT_FAILURE);
	}
	PrintStepSize=StepSize;
	//===================================================
	// Create the cvode integrator
  	//===================================================
	void * cvode_mem;
  	void * UserData;
  	SUNMatrix A						    = SUNDenseMatrix(NEQ, NEQ, sunctx);
  	SUNLinearSolver LS 					= SUNLinSol_SPGMR(y, PREC_NONE, 20, sunctx);
  	SUNNonlinearSolver NLS 				= SUNNonlinSol_Newton(y,sunctx);
	int retVal = 0;
	realtype tret=0;
	N_Vector AbsTol= N_VNew_Serial(NEQ, sunctx);
	for ( int i = 0 ; i < NEQ ; i++)
    NV_Ith_S(AbsTol,i)=absTol;

	cvode_mem = CVodeCreate (CV_BDF, sunctx);
	//Give CVODE the user problem
	retVal = CVodeSetUserData(cvode_mem, UserData);
	//Set intial state
	retVal = CVodeInit(cvode_mem, RHS, 0, y);
	//Set intial stepsize guess
	retVal = CVodeSetInitStep(cvode_mem, StepSize);
	//Set max stepsize
  	retVal = CVodeSetMaxStep(cvode_mem, 1);
  	//retVal = CVodeSetMaxStep(cvode_mem, StepSize);
	//retVal = CVodeSetMinStep(cvode_mem, StepSize);
	//Set tolerances
	retVal = CVodeSVtolerances(cvode_mem, relTol, AbsTol);
	//Set linear solver
	retVal = CVodeSetLinearSolver(cvode_mem, LS, A);
	retVal = CVodeSetNonlinearSolver(cvode_mem, NLS);
	retVal = CVodeSetJacTimes(cvode_mem, NULL, Jtv);
  	retVal = CVodeSetMaxNumSteps(cvode_mem, max(10000,NEQ*NEQ*NEQ));
	N_VDestroy_Serial(AbsTol);
  
  	const int MaxKrylovIters = 10;

	myPb problem;
	void *pbptr= &problem;


	IntegratorStats *Stats = NULL;
	Epi3VChem_KIOPS * EPI3V = new Epi3VChem_KIOPS(RHS, Jtv, 
		CVodeComputeJacWrapper, 
		pbptr, MaxKrylovIters, y, NEQ);

	//========================
	// Set integrator parameters
	//========================
	int startingBasis[] = {1, 3};
	//cout<<"=======================Initial Data======================\n";
	//Run the cvode integrator
	auto Start=std::chrono::high_resolution_clock::now();
	if(opt == "CVODE")
	{
		CVode(cvode_mem, FinalTime, y, &TNext, CV_NORMAL);
	}
	else if(opt == "EPI3V")
	{
		Stats = EPI3V->NewIntegrateNoTChem(StepSize, 1, absTol, relTol, 
			0.0, FinalTime, startingBasis, y);
	}
	auto Stop=std::chrono::high_resolution_clock::now();
	auto Pass = std::chrono::duration_cast<std::chrono::nanoseconds>(Stop-Start);
	KTime+=Pass.count()/1e9;
	realtype error_inf = 0;
 	realtype error_tmp = 0;
	
	// if(EPS == 1e-1)
	// {
	// 	error_inf = max( 
	// 		data[0] - 0.000000002988028,
	// 		max(
	// 			data[1] - 0.032702315935000,
	// 			data[2] - 1.967297681076978
	// 		));
    //   error_tmp = abs(data[2] - 1.967297681076978);
	// }
	// else if(EPS == 1e-2)
	// {
	// 	error_inf = max( 
	// 		abs(data[0] - 0.000000000000000),
	// 		max(
	// 			abs(data[1] - 0.005848132117459),
	// 			abs(data[2] - 1.994151867882557)
	// 		));
    // error_tmp = abs(data[2] - 1.994151867882557);
	// }

	//=======================================
	//Console Output
	//=======================================
	if(opt == "EPI3V")
	{
		Stats->PrintStats();
		StepCount 		= Stats->numTimeSteps;
	}
	else if(opt == "CVODE")
	{
		CVodeGetNumSteps(cvode_mem, &StepCount);
	}
	cout<<"===========================Data========================\n";
	cout << setprecision(17);
	cout <<"y=" << data[0] << "\t\t z=" << data[1] <<"\t\t Temp=" << data[2]<<endl;
	if (LastStepAdjusted == 1)
		cout<<"Altered final ";
	cout << "Simulation ended at time: " << setprecision(16) << TNext << endl;
	cout << "Steps taken: " << StepCount <<endl;
	cout << "Integration time: " <<KTime << endl;
	cout << "StepSize: " << StepSize << endl;
	cout <<"Error in sum=" << abs(data[0]+data[1]+data[2]-2) <<endl;
	cout <<"===================Printing data to " <<MyFile <<"==============\n";
	//=======================================
	//Data file output
	//=======================================
	ofstream myfile(MyFile, std::ios_base::app);
	myfile << setprecision(17) << fixed << data[0] << "\t\t" << data[1] << "\t\t";
	myfile << data[2] << "\t\t" << StepCount << "\t\t" << KTime;
	myfile << "\t\t" << absTol << "\t\t" << relTol << "\t\t" << abs(data[0]+data[1]+data[2]-2) <<endl;
	myfile.close();
	cout << "======================Simulation complete=============\n\n\n\n";


}



/*
//=========================================================
//  //====  ==   ==  ==     ==    _.==_    ========  _===_
//  ||      ||   ||  ||\\   ||  ./     \\     ||    //   \\
//  ||====  ||   ||  || \\  ||  ||            ||    \\___
//  ||      ||   ||  ||  \\ ||	~\            ||         \\
//  ||      \\___//  ||   \\||   ~=___-//     ||    \\___//
//=========================================================
*/

/*
 * ===========================================================================================
 * 
 * Function InitialConditions
 * 
 * Computes initial condition. 
 * 
 * Output:
 * y0   data
 * 
 * ===========================================================================================
 */

N_Vector InitialConditions()
{
  sundials::Context sunctx;   // SUNDIALS context
  N_Vector y0 = N_VNew_Serial(NEQ, sunctx);
  realtype *data = NV_DATA_S(y0);
	data[0]=1.0-.5*EPS;
	data[1]=.5*EPS;
	data[2]=1.0;
  return y0;
}


/*
 * ===========================================================================================
 * 
 * Function RHS
 * 
 * If y' = f(t, y) then RHS is function f. 
 * Function RHS is written is such way that it can handle arrays or arbitrarly length. 
 * It is important that all data is in form of an array. 
 * 
 * Inputs:
 * t          time
 * u          input vector
 * udot       result
 * userData   
 * 
 * ===========================================================================================
 */

int RHS(realtype t, N_Vector u, N_Vector udot, void *userData)
{
	realtype *y = NV_DATA_S(u);
	realtype *dy = NV_DATA_S(udot);
	double Omega= .5 * y[0] * y[1] * exp ( (y[2]-1.0 ) /  (EPS * y[2]  ) );
	dy[0] = -Omega;
	dy[1] = Omega-y[1];
	dy[2] = y[1];
	return 0;
}

/*
 * ===========================================================================================
 * 
 * Function Jtv
 * 
 * This function computes Jacobian matrix times some vector v. In Epirk Jacobian is never
 * stored and it is computed every time we need it and every time we compute it we acctualy
 * compute it we compute its product with some vector v. 
 * Function Jtv is written is such way that it can handle arrays or arbitrarly length. 
 * It is important that all data is in form of an array. 
 * It is used in JTimesv class. 
 * 
 * Inputs:
 * v          vector being multiplied with Jacobian Matrix
 * Jv         result
 * t          time
 * u          vector used ot compute Jacobian Matrix
 * fu         f(u), i.e., RHS of u
 * userData
 * tmp
 * 
 * ===========================================================================================
 */
int Jtv(N_Vector v, N_Vector Jv, realtype t, N_Vector u, N_Vector fu, void *userData, N_Vector tmp)
{
	realtype *JV = NV_DATA_S(Jv);
	realtype *Y = NV_DATA_S(u);
	realtype *V = NV_DATA_S(v);
  double EXP= exp ( (Y[2] - 1.0 ) / (EPS * Y[2]));
  double DEXP=1/(EPS * Y[2]*Y[2]);
	JV[0] = -.5*EXP* (V[0]*Y[1]+V[1]*Y[0]+V[2]*Y[0]*Y[1]*DEXP);
	JV[1] =-JV[0]-V[1];
	JV[2] = V[1];
  return 0;
}

int CVodeComputeJacWrapper(realtype t, N_Vector u, N_Vector fy, SUNMatrix Jac,
               void * pb, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
	ComputeJac(u, pb);
	return 0;
}

int ComputeJac(N_Vector y, void* UserData)
{
	return 0;
}