#include <string.h>
#include <sys/resource.h>
#include <map>
#include <regex>
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
static double gDW_CDmax=0;
static double gDW_IEC=0;
static map <string, string> FBAmet;
static unordered_map <string, double> Vmax;
static unordered_map <string, double> KM;
static double Ktoxin = 0;
static double drIECs = 0;

/* Read data from file and fill a map<string,int> */
void read_map_string_int(string fname, unordered_map<string,int>& m)
{
	ifstream f (fname);
	string line;
	if(f.is_open())
	{
		cout << "#### " << fname << "####" << endl;
		int j = 1;
		while (getline(f,line))
		{
			line.erase(remove( line.begin(), line.end(), '\"' ),line.end());
			m.insert(pair<string,int>(line,j));
			//cout << line << ";" << j << endl;
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

/* Read data from file and fill a map<string,double> */
void read_map_string_double(string fname, unordered_map<string,double>& m)
{
	ifstream f (fname);
	string line;
	if(f.is_open())
	{
		size_t pos = 0, c_pos = 0, length = 0;
		cout << "#### " << fname << "####" << endl;
		int j = 1;
		while (getline(f,line))
		{
			line.erase(remove( line.begin(), line.end(), '\"' ),line.end());
			pos = 0;
			c_pos = 0;
			// read rates
			length = line.length();

				pos = line.find(',');
				if( pos == string::npos)
					pos = length;
				m.insert(pair<string,double>(line.substr(0,pos) , stod(line.substr(pos+1,length))) );
				cout <<line.substr(0,pos) << ": " << stod(line.substr(pos+1,length)) << " ";
			cout << endl;
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
	read_constant("./gDW_CDmax", gDW_CDmax);
	read_constant("./gDW_IEC", gDW_IEC);

	read_map_string_double("./VmaxValues", Vmax);
	read_map_string_double("./KMValues", KM);
	read_map_string_double("./drIECs", drIECs);

	FBAmet["EX_biomass_e_in"] = "EX_biomass(e)";
	FBAmet["EX_pheme_e_in"] = "EX_pheme(e)";
	FBAmet["EX_cys_L_e_in"] = "EX_cys_L(e)";
	FBAmet["EX_trp_L_e_in"] = "EX_trp_L(e)";
	FBAmet["EX_val_L_e_in"] = "EX_val_L(e)";
	FBAmet["EX_ile_L_e_in"] = "EX_ile_L(e)";
	FBAmet["EX_leu_L_e_in"] = "EX_leu_L(e)";
	FBAmet["EX_pro_L_e_in"] = "EX_pro_L(e)";

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

			double Ub = l.getUpBounds(index);
			double Lb = l.getLwBounds(index);

			// if it is in -> updating the buonds!!
			if(p->first == "EX_biomass_e_in"){
				int Biom = Value[NumPlaces.find("BiomassCD") -> second];
				Ub = (gDW_CDmax - Biom);
				if( Ub <= 0 ) Ub = 0.000001;
			}else{
				double Met = Value[Trans[NumTrans.find(p->first) -> second].InPlaces[0].Id];
				cout<< "Trans: " <<  p->first << ", MEt input: " << Met <<";" << endl;
				Lb = (- (Vmax.find(p->first) -> second) * Met) / ((KM.find(p->first) -> second) + Met);
				}
			l.update_bound(index, TypeBound, Lb, Ub);
		}

		l.solve();
		Vars = l.getVariables();

		FBAtime = time;
	}

	int indexR = 0;
	bool Out = 1; // 0 if it is "_in" or in not reversible
	bool In = 0;

	string str = NameTrans[T];
	if(str.find("_out") != string::npos){
		str = std::regex_replace(str, std::regex("_out"), "_in");// replace '_out' -> '_in'
	}
	else{
		Out = 0;
		if(NameTrans[T].find("_in") != string::npos){
			In = 1;
		}
	}

	indexR = ReactionsNames.find(	FBAmet.find(str) -> second ) -> second ;
	cout<<"\nSolution:\n\n";
	//l.print();
	rate=Vars[indexR];

	// trans_in when is neg, otherwise trans_out
	if( (Out) & (rate > 0) )
		rate = rate;
	else if( (Out) & (rate < 0) )
		rate = 0 ;
	else if( (In) & (rate > 0) )
		rate = 0 ;
	else if( (In) & (rate < 0) )
		rate = -rate ;

	cout << NameTrans[T] << " transition with " <<
			FBAmet.find(str) -> second << " flux: " << Vars[indexR] <<
				"(rate: "<< rate << ")" << endl;

	if(rate < 0){
		cout << "WARNING: the rate is negative!!!! see transition: " << NameTrans[T] << endl;
		rate = -rate;
	}


	return(rate);
}

double IECDeath(double *Value,
           map <string,int>& NumTrans,
           map <string,int>& NumPlaces,
           const vector<string> & NameTrans,
           const struct InfTr* Trans,
           const int T,
           const double& time)
{

	if( Flag == -1)   init_data_structures();

	double rate = 0;
	double Den = BiomassCD + Ktoxin;
	double IECs = Value[NumPlaces.find("IECs") -> second];
	double BiomassCD = Value[NumPlaces.find("BiomassCD") -> second];

	(BiomassCD + drIECs)/Den * IECs * gDW_IEC
	(BiomassCD)/Den * IECs * gDW_IEC

}

