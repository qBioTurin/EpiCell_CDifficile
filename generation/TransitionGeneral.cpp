#include <string.h>
#include <sys/resource.h>
#include <map>
#include <string>
using namespace FBGLPK;

/*
 static double Flag = -1;
 static double Infection_rate;

 void read_constant(string fname, double& Infection_rate)
 {
 ifstream f (fname);
 string line;
 if(f.is_open())
 {
 int i = 1;
 while (getline(f,line))
 {
 switch(i)
 {
 case 1:
 Infection_rate = stod(line);
 //cout << "p" << i << ": " << line << "\t" << p1 << endl;
 break;
 }
 ++i;
 }
 f.close();
 }
 else
 {
 std::cerr<<"\nUnable to open " << fname << ": file do not exists\": file do not exists\n";
 exit(EXIT_FAILURE);
 }
 }
 */

//static double FlagFBA = -1; 1 when the FBA was already run at that time!
static double* Vars;
static double FBAtime = -1;
static unordered_map <string, int> ReactionsNames;
static double Flag = -1;
double rate = 0;
static double EXphemeCost=1;
static map <string, string> FBAmet;

/* Read data from file and fill a map<string,int> */
void read_map_string_int(string fname, unordered_map<string,int>& m)
{
	ifstream f (fname);
	string line;
	if(f.is_open())
	{
		cout << "#### " << fname << "####" << endl;
		int j = 0;
		while (getline(f,line))
		{
			line.erase(remove( line.begin(), line.end(), '\"' ),line.end());
			m.insert(pair<string,int>(line,j));
			cout << line << ";" << j << endl;
			++j;
		}
		f.close();
	}
	else
	{
		std::cerr<<"\nUnable to open " << fname << ": file do not exists\n";
		exit(EXIT_FAILURE);
	}
}

void read_constant(string fname, double& Infection_rate)
{
	ifstream f (fname);
	string line;
	if(f.is_open())
	{
		int i = 1;
		while (getline(f,line))
		{
			switch(i)
			{
			case 1:
				Infection_rate = stod(line);
				//cout << "p" << i << ": " << line << "\t" << p1 << endl;
				break;
			}
			++i;
		}
		f.close();
	}
	else
	{
		std::cerr<<"\nUnable to open " << fname << ": file do not exists\": file do not exists\n";
		exit(EXIT_FAILURE);
	}
}

string str = "FBAfile";
LPprob l(str.c_str());

void init_data_structures()
{
	read_map_string_int("./ReactNames", ReactionsNames);
	read_constant("./EX_phemeSens", EXphemeCost);

	FBAmet["EX_biomass_e"] = "EX_biomass(e)";
	FBAmet["EX_pheme_e_in"] = "EX_pheme(e)";
	FBAmet["EX_cys_L_e"] = "EX_cys_L(e)";
	FBAmet["EX_trp_L_e"] = "EX_trp_L(e)";
	FBAmet["EX_val_L_e"] = "EX_val_L(e)";
	FBAmet["EX_ile_L_e"] = "EX_ile_L(e)";
	FBAmet["EX_leu_L_e"] = "EX_leu_L(e)";
	FBAmet["EX_pro_L_e"] = "EX_pro_L(e)";

	Flag = 1;
}

double FBA(double *Value,
           map <string,int>& NumTrans,
           map <string,int>& NumPlaces,
           const vector<string> & NameTrans,
           const struct InfTr* Trans,
           const int T,
           const double& time)
{

	if( Flag == -1)   init_data_structures();

	// Definition of the function exploited to calculate the rate,
	// in this case for semplicity we define it throught the Mass Action law

	if( FBAtime != time){

		for (map<string, string>::iterator p = FBAmet.begin();
       p != FBAmet.end(); ++p ) {
			// Saving the reactions indexes from the map
			int index = ReactionsNames.find(p->second) -> second ;
			string TypeBound = "GLP_DB";

			double Lb = l.getLwBounds(index);
			double Ub = l.getUpBounds(index);
			// updating the buonds for the sensitivity!!

			Lb = Lb * EXphemeCost;
			Ub = Ub * EXphemeCost;

			l.update_bound(index, TypeBound, Lb, Ub);

		}

		l.solve();
		Vars = l.getVariables();

		FBAtime = time;
	}

	int indexR = 0;
	if(NameTrans[T] == "EX_pheme_e_out")
		indexR = ReactionsNames.find(	FBAmet.find("EX_pheme_e_in") -> second ) -> second ;
	else
		indexR = ReactionsNames.find(	FBAmet.find(NameTrans[T]) -> second ) -> second ;

	cout<<"\nSolution:\n\n";
	//l.print();
	rate=Vars[indexR];

	// trans_in when is neg, otherwise trans_out
	if((NameTrans[T] == "EX_pheme_e_out") & (rate > 0) )
		rate = rate;
	else if((NameTrans[T] == "EX_pheme_e_out") & (rate < 0 ))
		rate = 0 ;
	else if((NameTrans[T] == "EX_pheme_e_in") & (rate > 0) )
		rate = 0 ;
	else if((NameTrans[T] == "EX_pheme_e_in") & (rate < 0) )
		rate = -rate ;

	if(NameTrans[T] == "EX_pheme_e_out")
		cout << NameTrans[T] << " transition with " <<
			FBAmet.find("EX_pheme_e_in") -> second << " flux: " << Vars[indexR] <<
				"(rate: "<< rate << ")" << endl;
	else
		cout << NameTrans[T] << " transition with " <<
			FBAmet.find(NameTrans[T]) -> second << " flux: " << Vars[indexR] <<
				"(rate: "<< rate << ")" << endl;



	if(rate < 0)
		rate = -rate;

	return(1);
}
