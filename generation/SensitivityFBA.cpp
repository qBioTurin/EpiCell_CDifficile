
#include <iostream>
#include "class.hpp"


#include <iostream>
#include "/home/SensitivityFBA.hpp"

using namespace SDE;
extern double epsilon;

 string places[]={"P0"};
 string transitions[]={"T0"};

int main(int argc, char **argv) {

 time_t time_1,time_4;
 int who = RUSAGE_SELF;
 struct rusage usage;
 int SOLVE = 7, runs=1;
 long int seed = 0;
 bool OUTPUT=false;
 std::string fbound="", finit="", fparm="";
 double hini = 1e-6, atolODE = 1e-6, rtolODE=1e-6, ftime=1.0, stime=0.0, itime=0.0, epsTAU=0.1;

 cout<<"\n\n =========================================================\n";
 cout<<"|	              ODE/SDE Solver                       |\n";
 cout<<" =========================================================\n";
 cout<<"\n If you find any bug, send an email to beccuti@di.unito.it\n";

 if (argc<2)
	{
 	 std::cerr<<"\n\nUSE:SensitivityFBA_solve <out_file> [OPTIONS]";
	std::cerr<<"\n\n\tOPTIONS\n";
	 std::cerr<<"\n\t -type <type>:\t\t ODE-E or ODE-RKF or ODE45 or LSODA or HLSODA or (H)SDE or HODE or SSA or TAUG or STEP. Default: LSODA ";
	 std::cerr<<"\n\t -hini <double>:\t Initial step size. Default: 1e-6";
	 std::cerr<<"\n\t -atol <double>:\t Absolute error tolerance. Dafault: 1e-6";
	 std::cerr<<"\n\t -rtol <double>:\t Relative error tolerance. Dafault: 1e-6";
	 std::cerr<<"\n\t -taueps <double>:\t Epsilon value for Tau-leaping algorithm. Dafault: 0.1";
	 std::cerr<<"\n\t -runs <int>:\t\t Integer number corresponding to runs (only used in SSA,TAUG, HODE,HLSODA). Default: 1";
	 std::cerr<<"\n\t -ftime <double>:\t Double number used to set the upper bound of the evolution time. Dafault: 1";
	 std::cerr<<"\n\t -stime <double>:\t Double number used to set the step in the output. Default: 0.0 (no output)";
	 std::cerr<<"\n\t -itime <double>:\t Double number used to set the initial simulation time. Default: 0.0 ";
	 std::cerr<<"\n\t -b <bound_file>:\t Soft bound are defined in the file <bound_file>";
	 std::cerr<<"\n\t -seed <double>:\t Seed of random number generator";
	 std::cerr<<"\n\t -init <init_file>:\t The file <initial_file> contains the initial marking. Default:  initial marking in the orginal net";
	 std::cerr<<"\n\t -parm <parm_file>:\t The file <parm_file> contains a set of pairs with format <transition name> <value> or <place name> <value>.\n\t\t\t\t For transition  the value is used to set a new rate value, while for place  it is used to set a new initial marking.";
	 std::cerr<<endl<<endl;
	 exit(EXIT_FAILURE);
	}

 int ii=2;
 for (; ii<argc; ii++){
	 if (strcmp("-type", argv[ii])==0){
		 if (++ii<argc){
			 if ((strcmp(argv[ii],"ODE-E")==0)||(strcmp(argv[ii],"ode-e")==0)) SOLVE = 1;
			 else if ((strcmp(argv[ii],"ODE-RKF")==0)||(strcmp(argv[ii],"ode-rkf")==0)) SOLVE = 5;
			 else if ((strcmp(argv[ii],"ODE45")==0)||(strcmp(argv[ii],"ode45")==0)) SOLVE = 6;
			 else if ((strcmp(argv[ii],"LSODA")==0)||(strcmp(argv[ii],"lsoda")==0)) SOLVE = 7;
			 else if ((strcmp(argv[ii],"STEP")==0)||(strcmp(argv[ii],"step")==0)) SOLVE = 4;
			 else if ((strcmp(argv[ii],"SSA")==0)||(strcmp(argv[ii],"ssa")==0)){

				 SOLVE=3;
			 }
			 else if ((strcmp(argv[ii],"HODE")==0)||(strcmp(argv[ii],"hode")==0)) SOLVE = 2;
			 else if ((strcmp(argv[ii],"HLSODA")==0)||(strcmp(argv[ii],"hlsoda")==0)) SOLVE = 8;
			 else if ((strcmp(argv[ii],"SDE")==0)||(strcmp(argv[ii],"sde")==0) || (strcmp(argv[ii],"HSDE")==0)||(strcmp(argv[ii],"hsde")==0) ) SOLVE = 0;
			 else if ((strcmp(argv[ii],"TAUG")==0)||(strcmp(argv[ii],"taug")==0)) SOLVE = 9;
			 else{
				 std::cerr<<"\n\tError:  -type  <value>\n\n";
			 exit(EXIT_FAILURE);
		 }
		 }
		 continue;
	 }
	 if (strcmp("-hini", argv[ii])==0){
		 if (++ii<argc){
			 hini=atof(argv[ii]);
		 }
		 else{
			 std::cerr<<"\nError:  -hini  <value>\n";
			 exit(EXIT_FAILURE);
		 }
		 continue;
	 }
	 if (strcmp("-atol", argv[ii])==0){
		 if (++ii<argc){
			 atolODE=atof(argv[ii]);
		 }
		 else{
			 std::cerr<<"\nError:  -atol  <value>\n";
			 exit(EXIT_FAILURE);
		 }
		 continue;
	 }
	 if (strcmp("-rtol", argv[ii])==0){
		 if (++ii<argc){
			 rtolODE=atof(argv[ii]);
		 }
		 else{
			 std::cerr<<"\nError:  -rtol  <value>\n";
			 exit(EXIT_FAILURE);
		 }
		 continue;
	 }
	 if (strcmp("-runs", argv[ii])==0){
		 if (++ii<argc){
			 runs=atoi(argv[ii]);
		 }
		 else{
			 std::cerr<<"\nError:  -runs  <value>\n";
			 exit(EXIT_FAILURE);
		 }
		 continue;
	 }
	 if (strcmp("-ftime", argv[ii])==0){
		 if (++ii<argc){
			 ftime=atof(argv[ii]);
		 }
		 else{
			 std::cerr<<"\nError:  -ftime  <value>\n";
			 exit(EXIT_FAILURE);
		 }
		 continue;
	 }
	 if (strcmp("-stime", argv[ii])==0){
		 if (++ii<argc){
			 OUTPUT=true;
			 stime=atof(argv[ii]);
		 }
		 else{
			 std::cerr<<"\nError:  -stime  <value>\n";
			 exit(EXIT_FAILURE);
		 }
		 continue;
	 }
	 if (strcmp("-itime", argv[ii])==0){
		 if (++ii<argc){
			 itime=atof(argv[ii]);
		 }
		 else{
			 std::cerr<<"\nError:  -itime <value>\n";
			 exit(EXIT_FAILURE);
		 }
		 continue;
	 }
	 if (strcmp("-b", argv[ii])==0){
		 if (++ii<argc){
			 fbound=string(argv[ii]);
		 }
		 else{
			 std::cerr<<"\nError:  -b  <file_name>\n";
			 exit(EXIT_FAILURE);
		 }
		 continue;
	 }
	 if (strcmp("-taueps", argv[ii])==0){
		 if (++ii<argc){
			 epsTAU=atof(argv[ii]);
		 }
		 else{
			 std::cerr<<"\nError:  -taueps  <value>\n";
			 exit(EXIT_FAILURE);
		 }
		 continue;
	 }
	 if (strcmp("-seed", argv[ii])==0){
		 if (++ii<argc){
			 seed=atol(argv[ii]);
		 }
		 else{
			 std::cerr<<"\nError:  -seed  <value>\n";
			 exit(EXIT_FAILURE);
		 }
		 continue;
	 }
	 if (strcmp("-init", argv[ii])==0){
		 if (++ii<argc){
			 finit=string(argv[ii]);
		 }
		 else{
			 std::cerr<<"\nError:  -init  <file_name>\n";
			 exit(EXIT_FAILURE);
		 }
		 continue;
	 }
	 if (strcmp("-parm", argv[ii])==0){
		 if (++ii<argc){
			 fparm=string(argv[ii]);
		 }
		 else{
			 std::cerr<<"\nError:  -parm  <file_name>\n";
			 exit(EXIT_FAILURE);
		 }
		 continue;
	 }
			 std::cerr<<"\nError:  unknown parameter "<<argv[ii]<<"\n\n";
			 exit(EXIT_FAILURE);
 }


 if (stime==0.0)  stime=ftime;

 time(&time_1);

 cout<<"\n=====================INPUT PARAMETERS======================\n";
 cout<<"\n\tCompact CPP code: OFF\n";
 cout<<"\tType solution: "<<SOLVE<<"\n";
 cout<<"\tTransition policy: Genelarized Mass Action policy\n";
 cout<<"\tSolution final time: "<<ftime<<"\n";
 cout<<"\tInitial size step: "<<hini<<"\n";
 cout<<"\tInitial  time: "<<itime<<"\n";
 cout<<"\tAbosolute tolerance: "<<atolODE<<"\n";
 cout<<"\tRelative tolerance: "<<rtolODE<<"\n";
 cout<<"\tEpsilon value for TAU-leaping: "<<epsTAU<<"\n";
 cout<<"\tSolution runs: "<<runs<<"\n";
 if (fbound!="") cout<<"\tBound file: "<<fbound<<"\n";
 if (finit!="") cout<<"\tInitial marking file: "<<finit<<"\n";
 if (fparm!="") cout<<"\tInitial parameter file: "<<fparm<<"\n";
 cout<<"\tDetailed output: "<<stime<<"\n";
 cout<<"\n===========================================================\n"; cout<<"\n\nSTART EXECUTION..."<<endl;

 struct InfPlace pt;
 struct InfTr t;
 Equation eq;
 Elem el;
 SystEqMas se(1,1,places,transitions,itime,seed);
 vector< struct InfPlace> Vpl;

//Transition T0
 t.InPlaces.clear();
 t.InhPlaces.clear();
 t.InOuPlaces.clear();
 t.Places.clear();
 t.discrete= false;
 t.GenFun= "FBA";
 t.FuncT=  &FBA;
 t.rate = 1.0;
 t.InOuPlaces.insert(0);
 pt.Id = 0;
 pt.Card = 1;
 t.Places.push_back(pt);
 se.InsertTran(0,t);

//Place P0
 eq.clear();
 el.setIncDec(1);
 el.setIdTran(0);
 eq.Insert(el);
 se.InsertEq(0,eq,0,0,2147483647);


 if (fbound!="") {
	 if (!(se.readSLUBounds(fbound))) exit(EXIT_FAILURE);;
 }
 if (finit!="") {
	 if (!(se.readInitialMarking(finit))) exit(EXIT_FAILURE);
 }
 if (fparm!="") {
	 if (!(se.readParameter(fparm))) exit(EXIT_FAILURE);
 }
 se.setEpsTAU(epsTAU);

 se.Print();


try	{
	if (SOLVE==-1) {
		 cerr<< "\n\nError: solution methods is not implemented\nYou should use:  ODE-E or ODE-RKF or ODE45 or LSODA or SDE or HODE or HSDE or TAUG or SSA or STEP\n"; 
		 exit(EXIT_FAILURE);
	}

 
	if (SOLVE == 1)
		 se.SolveODEEuler(hini,atolODE,rtolODE,ftime,OUTPUT,stime,argv[1]);
	 else
		 if (SOLVE == 0)
			 se.SolveSDEEuler(hini,atolODE,rtolODE,ftime,runs,OUTPUT,stime,argv[1]);
		 else 
			if (SOLVE == 3)
				 se.SolveSSA(hini,atolODE,rtolODE,ftime,runs,OUTPUT,stime,argv[1]); 
			 else 
				 if (SOLVE == 4)
					  se.HeuristicStep(hini,atolODE,rtolODE,ftime,OUTPUT,stime,argv[1]);   
				 else
					 if (SOLVE == 5)
					  se.SolveODERKF(hini,atolODE,ftime,OUTPUT,stime,argv[1]);   
				 else
					if (SOLVE == 6)
						 se.SolveODE45(hini,atolODE,ftime,OUTPUT,stime,argv[1]);
				 else
					 if (SOLVE == 8)
							 se.SolveHLSODE(hini,atolODE,rtolODE,ftime,runs,OUTPUT,stime,argv[1]);
					 else 
							 if (SOLVE == 7) 
								 se.SolveLSODE(hini,atolODE,rtolODE,ftime,OUTPUT,stime,argv[1]);
							 else  
								 se.SolveTAUG(ftime,runs,OUTPUT,stime,argv[1]);
	se.PrintStatistic(argv[1]);
	}
 catch(Exception obj)
	{
	cerr<<endl<<obj.get()<<endl;
	}

 time(&time_4);

 cout<<"\n\nEND EXECUTION"<<endl;
 cout<<"\nResults are saved in: "<<argv[1]<<endl;
 cout<<"\n=========================== TIME ===========================\n\n\t";
 cout<<"Total time required: "<<(time_4-time_1)<<"s."<<endl;
 cout<<"\n=========================== TIME ===========================\n\n";
 cout<<"\n=========================== MEM. ===========================\n\n\t";
 getrusage(who,&usage);
 cout<<"Total memory used: "<<usage.ru_maxrss<<"KB"<<endl;
 cout<<"\n=========================== TIME ===========================\n\n";

}