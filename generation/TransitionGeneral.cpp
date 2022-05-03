#include <string.h>
#include <sys/resource.h>

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

string str = "FBAfile";
LPprob l(str.c_str());

void init_data_structures()
{
	read_map_string_int("./ReactNames", ReactionsNames);
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
	double rate = 0;
	// Saving the reactions indexes from the map
	int indexR = ReactionsNames.find("EX_pheme(e)") -> second ;
	//char cstr[str.size()+1];
	//strcpy(cstr,str.c_str());

	// Definition of the function exploited to calculate the rate,
	// in this case for semplicity we define it throught the Mass Action law

	if( FBAtime != 0){

		//string TypeBound = "GLP_DB";
		//double Lb = 0;
		//double Ub = 99;

		//l.update_bound(indexR, TypeBound, Lb, Ub);

		l.solve();

		cout<<"\nSolution:\n\n";
		//l.print();
		Vars = l.getVariables();
		cout << "EX_pheme(e) flux: " << Vars[indexR] << endl;
		FBAtime = 0;
		rate=Vars[indexR];
	}


	return(rate);
}
