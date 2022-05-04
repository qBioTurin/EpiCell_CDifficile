
#include <iostream>
#include "class.hpp"


#include <iostream>
#include "/home/EpitCellDifficile.hpp"

using namespace SDE;
extern double epsilon;

 string places[]={"EpithelialCells","pheme_e","BiomassCD","pro_L_e","leu_L_e","ile_L_e","val_L_e","trp_e","cys_e","Drug","OxiStress","Damage","pro_L_v","leu_L_v","ile_L_v","val_L_v","trp_L_v","cys_L_v"};
 string transitions[]={"EX_pheme_e_in","EX_pro_L_e","EX_biomass_e","Death","DeathCD","EX_val_L_e","EX_ile_L_e","EX_leu_L_e","EX_cys_L_e","EX_trp_L_e","DrugAction","T_cys_L_e","EX_pheme_e_out","T_trp_L_e","T_val_L_e","T_ile_L_e","T_leu_L_e","T_pro_L_e"};

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
 	 std::cerr<<"\n\nUSE:EpitCellDifficile_solve <out_file> [OPTIONS]";
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
 SystEqMas se(18,18,places,transitions,itime,seed);
 vector< struct InfPlace> Vpl;

//Transition EX_pheme_e_in
 t.InPlaces.clear();
 t.InhPlaces.clear();
 t.InOuPlaces.clear();
 t.Places.clear();
 t.discrete= false;
 t.GenFun= "FBA";
 t.FuncT=  &FBA;
 t.rate = 1.0;
 pt.Id = 1;
 pt.Card = 1;
 t.InPlaces.push_back(pt);
 t.InOuPlaces.insert(1);
 pt.Id = 1;
 pt.Card = -1;
 t.Places.push_back(pt);
 se.InsertTran(0,t);

//Transition EX_pro_L_e
 t.InPlaces.clear();
 t.InhPlaces.clear();
 t.InOuPlaces.clear();
 t.Places.clear();
 t.discrete= false;
 t.GenFun= "FBA";
 t.FuncT=  &FBA;
 t.rate = 1.0;
 pt.Id = 3;
 pt.Card = 1;
 t.InPlaces.push_back(pt);
 t.InOuPlaces.insert(3);
 pt.Id = 3;
 pt.Card = -1;
 t.Places.push_back(pt);
 se.InsertTran(1,t);

//Transition EX_biomass_e
 t.InPlaces.clear();
 t.InhPlaces.clear();
 t.InOuPlaces.clear();
 t.Places.clear();
 t.discrete= false;
 t.GenFun= "FBA";
 t.FuncT=  &FBA;
 t.rate = 1.0;
 pt.Id = 2;
 pt.Card = 1;
 t.InPlaces.push_back(pt);
 t.InOuPlaces.insert(2);
 t.InOuPlaces.insert(2);
 se.InsertTran(2,t);

//Transition Death
 t.InPlaces.clear();
 t.InhPlaces.clear();
 t.InOuPlaces.clear();
 t.Places.clear();
 t.discrete= false;
 t.GenFun= "";
 t.FuncT=  nullptr;
 t.rate = 1.000000;
 pt.Id = 0;
 pt.Card = 1;
 t.InPlaces.push_back(pt);
 t.InOuPlaces.insert(0);
 pt.Id = 0;
 pt.Card = -1;
 t.Places.push_back(pt);
 se.InsertTran(3,t);

//Transition DeathCD
 t.InPlaces.clear();
 t.InhPlaces.clear();
 t.InOuPlaces.clear();
 t.Places.clear();
 t.discrete= false;
 t.GenFun= "";
 t.FuncT=  nullptr;
 t.rate = 1.000000;
 pt.Id = 10;
 pt.Card = 1;
 t.InPlaces.push_back(pt);
 t.InOuPlaces.insert(10);
 pt.Id = 10;
 pt.Card = -1;
 t.Places.push_back(pt);
 se.InsertTran(4,t);

//Transition EX_val_L_e
 t.InPlaces.clear();
 t.InhPlaces.clear();
 t.InOuPlaces.clear();
 t.Places.clear();
 t.discrete= false;
 t.GenFun= "FBA";
 t.FuncT=  &FBA;
 t.rate = 1.0;
 pt.Id = 6;
 pt.Card = 1;
 t.InPlaces.push_back(pt);
 t.InOuPlaces.insert(6);
 pt.Id = 6;
 pt.Card = -1;
 t.Places.push_back(pt);
 se.InsertTran(5,t);

//Transition EX_ile_L_e
 t.InPlaces.clear();
 t.InhPlaces.clear();
 t.InOuPlaces.clear();
 t.Places.clear();
 t.discrete= false;
 t.GenFun= "FBA";
 t.FuncT=  &FBA;
 t.rate = 1.0;
 pt.Id = 5;
 pt.Card = 1;
 t.InPlaces.push_back(pt);
 t.InOuPlaces.insert(5);
 pt.Id = 5;
 pt.Card = -1;
 t.Places.push_back(pt);
 se.InsertTran(6,t);

//Transition EX_leu_L_e
 t.InPlaces.clear();
 t.InhPlaces.clear();
 t.InOuPlaces.clear();
 t.Places.clear();
 t.discrete= false;
 t.GenFun= "FBA";
 t.FuncT=  &FBA;
 t.rate = 1.0;
 pt.Id = 4;
 pt.Card = 1;
 t.InPlaces.push_back(pt);
 t.InOuPlaces.insert(4);
 pt.Id = 4;
 pt.Card = -1;
 t.Places.push_back(pt);
 se.InsertTran(7,t);

//Transition EX_cys_L_e
 t.InPlaces.clear();
 t.InhPlaces.clear();
 t.InOuPlaces.clear();
 t.Places.clear();
 t.discrete= false;
 t.GenFun= "FBA";
 t.FuncT=  &FBA;
 t.rate = 1.0;
 pt.Id = 8;
 pt.Card = 1;
 t.InPlaces.push_back(pt);
 t.InOuPlaces.insert(8);
 pt.Id = 8;
 pt.Card = -1;
 t.Places.push_back(pt);
 se.InsertTran(8,t);

//Transition EX_trp_L_e
 t.InPlaces.clear();
 t.InhPlaces.clear();
 t.InOuPlaces.clear();
 t.Places.clear();
 t.discrete= false;
 t.GenFun= "FBA";
 t.FuncT=  &FBA;
 t.rate = 1.0;
 pt.Id = 7;
 pt.Card = 1;
 t.InPlaces.push_back(pt);
 t.InOuPlaces.insert(7);
 pt.Id = 7;
 pt.Card = -1;
 t.Places.push_back(pt);
 se.InsertTran(9,t);

//Transition DrugAction
 t.InPlaces.clear();
 t.InhPlaces.clear();
 t.InOuPlaces.clear();
 t.Places.clear();
 t.discrete= false;
 t.GenFun= "";
 t.FuncT=  nullptr;
 t.rate = 1.000000;
 pt.Id = 9;
 pt.Card = 1;
 t.InPlaces.push_back(pt);
 t.InOuPlaces.insert(9);
 t.InOuPlaces.insert(10);
 pt.Id = 9;
 pt.Card = -1;
 t.Places.push_back(pt);
 pt.Id = 10;
 pt.Card = 1;
 t.Places.push_back(pt);
 se.InsertTran(10,t);

//Transition T_cys_L_e
 t.InPlaces.clear();
 t.InhPlaces.clear();
 t.InOuPlaces.clear();
 t.Places.clear();
 t.discrete= false;
 t.GenFun= "";
 t.FuncT=  nullptr;
 t.rate = 1.000000;
 pt.Id = 0;
 pt.Card = 1;
 t.InPlaces.push_back(pt);
 t.InOuPlaces.insert(0);
 t.InOuPlaces.insert(0);
 t.InOuPlaces.insert(17);
 pt.Id = 17;
 pt.Card = 1;
 t.Places.push_back(pt);
 se.InsertTran(11,t);

//Transition EX_pheme_e_out
 t.InPlaces.clear();
 t.InhPlaces.clear();
 t.InOuPlaces.clear();
 t.Places.clear();
 t.discrete= false;
 t.GenFun= "FBA";
 t.FuncT=  &FBA;
 t.rate = 1.0;
 t.InOuPlaces.insert(1);
 pt.Id = 1;
 pt.Card = 1;
 t.Places.push_back(pt);
 se.InsertTran(12,t);

//Transition T_trp_L_e
 t.InPlaces.clear();
 t.InhPlaces.clear();
 t.InOuPlaces.clear();
 t.Places.clear();
 t.discrete= false;
 t.GenFun= "";
 t.FuncT=  nullptr;
 t.rate = 1.000000;
 pt.Id = 0;
 pt.Card = 1;
 t.InPlaces.push_back(pt);
 t.InOuPlaces.insert(0);
 t.InOuPlaces.insert(0);
 t.InOuPlaces.insert(16);
 pt.Id = 16;
 pt.Card = 1;
 t.Places.push_back(pt);
 se.InsertTran(13,t);

//Transition T_val_L_e
 t.InPlaces.clear();
 t.InhPlaces.clear();
 t.InOuPlaces.clear();
 t.Places.clear();
 t.discrete= false;
 t.GenFun= "";
 t.FuncT=  nullptr;
 t.rate = 1.000000;
 pt.Id = 0;
 pt.Card = 1;
 t.InPlaces.push_back(pt);
 t.InOuPlaces.insert(0);
 t.InOuPlaces.insert(0);
 t.InOuPlaces.insert(15);
 pt.Id = 15;
 pt.Card = 1;
 t.Places.push_back(pt);
 se.InsertTran(14,t);

//Transition T_ile_L_e
 t.InPlaces.clear();
 t.InhPlaces.clear();
 t.InOuPlaces.clear();
 t.Places.clear();
 t.discrete= false;
 t.GenFun= "";
 t.FuncT=  nullptr;
 t.rate = 1.000000;
 pt.Id = 0;
 pt.Card = 1;
 t.InPlaces.push_back(pt);
 t.InOuPlaces.insert(0);
 t.InOuPlaces.insert(0);
 t.InOuPlaces.insert(14);
 pt.Id = 14;
 pt.Card = 1;
 t.Places.push_back(pt);
 se.InsertTran(15,t);

//Transition T_leu_L_e
 t.InPlaces.clear();
 t.InhPlaces.clear();
 t.InOuPlaces.clear();
 t.Places.clear();
 t.discrete= false;
 t.GenFun= "";
 t.FuncT=  nullptr;
 t.rate = 1.000000;
 pt.Id = 0;
 pt.Card = 1;
 t.InPlaces.push_back(pt);
 t.InOuPlaces.insert(0);
 t.InOuPlaces.insert(0);
 t.InOuPlaces.insert(13);
 pt.Id = 13;
 pt.Card = 1;
 t.Places.push_back(pt);
 se.InsertTran(16,t);

//Transition T_pro_L_e
 t.InPlaces.clear();
 t.InhPlaces.clear();
 t.InOuPlaces.clear();
 t.Places.clear();
 t.discrete= false;
 t.GenFun= "";
 t.FuncT=  nullptr;
 t.rate = 1.000000;
 pt.Id = 0;
 pt.Card = 1;
 t.InPlaces.push_back(pt);
 t.InOuPlaces.insert(0);
 t.InOuPlaces.insert(0);
 t.InOuPlaces.insert(12);
 pt.Id = 12;
 pt.Card = 1;
 t.Places.push_back(pt);
 se.InsertTran(17,t);

//Place EpithelialCells
 eq.clear();
 el.setIncDec(-1);
 el.setIdTran(3);
 eq.Insert(el);
 se.InsertEq(0,eq,0,0,2147483647);

//Place pheme_e
 eq.clear();
 el.setIncDec(-1);
 el.setIdTran(0);
 eq.Insert(el);
 el.setIncDec(1);
 el.setIdTran(12);
 eq.Insert(el);
 se.InsertEq(1,eq,0,0,2147483647);

//Place BiomassCD
 eq.clear();
 se.InsertEq(2,eq,0,0,0);

//Place pro_L_e
 eq.clear();
 el.setIncDec(-1);
 el.setIdTran(1);
 eq.Insert(el);
 se.InsertEq(3,eq,0,0,2147483647);

//Place leu_L_e
 eq.clear();
 el.setIncDec(-1);
 el.setIdTran(7);
 eq.Insert(el);
 se.InsertEq(4,eq,0,0,2147483647);

//Place ile_L_e
 eq.clear();
 el.setIncDec(-1);
 el.setIdTran(6);
 eq.Insert(el);
 se.InsertEq(5,eq,0,0,2147483647);

//Place val_L_e
 eq.clear();
 el.setIncDec(-1);
 el.setIdTran(5);
 eq.Insert(el);
 se.InsertEq(6,eq,0,0,2147483647);

//Place trp_e
 eq.clear();
 el.setIncDec(-1);
 el.setIdTran(9);
 eq.Insert(el);
 se.InsertEq(7,eq,0,0,2147483647);

//Place cys_e
 eq.clear();
 el.setIncDec(-1);
 el.setIdTran(8);
 eq.Insert(el);
 se.InsertEq(8,eq,0,0,2147483647);

//Place Drug
 eq.clear();
 el.setIncDec(-1);
 el.setIdTran(10);
 eq.Insert(el);
 se.InsertEq(9,eq,0,0,2147483647);

//Place OxiStress
 eq.clear();
 el.setIncDec(-1);
 el.setIdTran(4);
 eq.Insert(el);
 el.setIncDec(1);
 el.setIdTran(10);
 eq.Insert(el);
 se.InsertEq(10,eq,0,0,2147483647);

//Place Damage
 eq.clear();
 se.InsertEq(11,eq,0,0,0);

//Place pro_L_v
 eq.clear();
 el.setIncDec(1);
 el.setIdTran(17);
 eq.Insert(el);
 se.InsertEq(12,eq,0,0,2147483647);

//Place leu_L_v
 eq.clear();
 el.setIncDec(1);
 el.setIdTran(16);
 eq.Insert(el);
 se.InsertEq(13,eq,0,0,2147483647);

//Place ile_L_v
 eq.clear();
 el.setIncDec(1);
 el.setIdTran(15);
 eq.Insert(el);
 se.InsertEq(14,eq,0,0,2147483647);

//Place val_L_v
 eq.clear();
 el.setIncDec(1);
 el.setIdTran(14);
 eq.Insert(el);
 se.InsertEq(15,eq,0,0,2147483647);

//Place trp_L_v
 eq.clear();
 el.setIncDec(1);
 el.setIdTran(13);
 eq.Insert(el);
 se.InsertEq(16,eq,0,0,2147483647);

//Place cys_L_v
 eq.clear();
 el.setIncDec(1);
 el.setIdTran(11);
 eq.Insert(el);
 se.InsertEq(17,eq,0,0,2147483647);


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