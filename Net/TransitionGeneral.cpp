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
static double gDW_CDmax = 0;
static double gDW_IEC = 0;
static map <string, string> FBAmet;
static unordered_map <string, double> Vmax;
static unordered_map <string, double> KM;
// Constants for Inflammation transition
static double Inflammation = 0;
static double DAMAGEmax = 0;
// Constants for Death4Treat transition
static double Emax = 0;
static double Imax = 0;
static double Mmtz = 0;
static double Na = 0;
static double Drug50 = 0;
static double K50 = 0;
// Constants for DeathBac transition
static double gDW_CDmin = 0;
static double h = 0;
// Constants for Dup transition
static double rCDdup = 0;
// Constants for Starv transition
static double RCD = 0;

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
		size_t pos = 0, length = 0;
		cout << "#### " << fname << "####" << endl;
		int j = 1;
		while (getline(f,line))
		{
			line.erase(remove( line.begin(), line.end(), '\"' ),line.end());
			pos = 0;
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
	read_constant("./gDW_CDmin", gDW_CDmin);

	read_constant("./gDW_IEC", gDW_IEC);

	read_constant("./Inflammation_rate", Inflammation);
	read_constant("./DAMAGEmax_rate", DAMAGEmax);

	read_constant("./Emax", Emax);
	read_constant("./Imax", Imax);
	read_constant("./Mmtz",Mmtz);
	read_constant("./Na", Na);
	read_constant("./Drug50", Drug50);
	read_constant("./K50", K50);

	read_constant("./h", h);

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
				Lb = (- (Vmax.find(p->first) -> second) * (Met/Na)) / ((KM.find(p->first) -> second) + (Met/Na));
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

// Inflam transition
double Heam(double *Value,
            map <string,int>& NumTrans,
            map <string,int>& NumPlaces,
            const vector<string> & NameTrans,
            const struct InfTr* Trans,
            const int T,
            const double& time)
{

	if( Flag == -1)   init_data_structures();

	double rate = 0;

	int DamagePlace = Value[NumPlaces.find("Damage") -> second];

	double PercDamage = DamagePlace/DAMAGEmax;
	double g = 0;

	if((PercDamage < 0) & (PercDamage > 1)){

	}
	else if((PercDamage>0.1) & (PercDamage <= 0.7)){
		g = 1/3;
	}
	else if(PercDamage > 0.7){
		g = 1;
	}

	rate = g * Inflammation;
	return(rate);
}

// Death4Treat transition
double Therapy(double *Value,
               map <string,int>& NumTrans,
               map <string,int>& NumPlaces,
               const vector<string> & NameTrans,
               const struct InfTr* Trans,
               const int T,
               const double& time)
{

	if( Flag == -1)   init_data_structures();

	double rate = 0;

	int DrugPlace = Value[NumPlaces.find("Drug") -> second];
	int CDPlace = Value[NumPlaces.find("CD") -> second];
	int HemePlace = Value[NumPlaces.find("pheme_e") -> second];

	rate = (Emax*DrugPlace*CDPlace)/(((Mmtz*Na)/(Drug50*((Imax*HemePlace)/(K50 + HemePlace))))+DrugPlace);

	return(rate);

}

// DeathBac transition
double DeathCD(double *Value,
               map <string,int>& NumTrans,
               map <string,int>& NumPlaces,
               const vector<string> & NameTrans,
               const struct InfTr* Trans,
               const int T,
               const double& time)
{

	if( Flag == -1)   init_data_structures();

	double rate = 0;

	int BiomassCDPlace = Value[NumPlaces.find("BiomassCD") -> second];
	int CDPlace = Value[NumPlaces.find("CD") -> second];

	if((BiomassCDPlace - gDW_CDmin) > 0){
		rate = 0;
	}

	rate = h*CDPlace;

	return(rate);

}

// Dup transition
double Duplication(double *Value,
                   map <string,int>& NumTrans,
                   map <string,int>& NumPlaces,
                   const vector<string> & NameTrans,
                   const struct InfTr* Trans,
                   const int T,
                   const double& time)
{

	if( Flag == -1)   init_data_structures();

	double rate = 0;

	int BiomassCDPlace = Value[NumPlaces.find("BiomassCD") -> second];

	rate = rCDdup*(1 - exp(BiomassCDPlace - gDW_CDmax));

	return(rate);

}

// Starv transition
double Starvation(double *Value,
                  map <string,int>& NumTrans,
                  map <string,int>& NumPlaces,
                  const vector<string> & NameTrans,
                  const struct InfTr* Trans,
                  const int T,
                  const double& time)
{

	if( Flag == -1)   init_data_structures();

	double rate = 0;

	int BiomassCDPlace = Value[NumPlaces.find("BiomassCD") -> second];

	double fbaEXBiomassOut = 1.507725e-06;

	if(fbaEXBiomassOut >= 0){
		rate = RCD*BiomassCDPlace;
	}
		rate = 0;

	return(rate);

}
