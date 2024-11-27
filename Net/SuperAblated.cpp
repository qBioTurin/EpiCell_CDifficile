
#include <string.h>
#include <sys/resource.h>
#include <map>
#include <regex>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <string>
#include<sstream>


static map <string, string> FBAmet;
static map <string, string> FBAplace;
static map <string, double> FBAvars;

static double Flag = -1;

static double FlagDebug = 1;

static double gDW_CDmean = 0;
static double gDW_CDmin = 0;
static double gDW_CDmax = 0;

// https://doi.org/10.1016/j.ymben.2021.10.012
// biomass molecular weigth [g/mmol]
// biomass defined to have a molecular weight (MW) of 1 (g/mmol)
static double MWbio = 0;

static double nBacMax = 0;

// chemical proportionality factors
static double Na = 0;
static double c = 0;

// multiplicative constant for sensitivity analysis in FBA 
static double P = 0;

// Maximum biomass growth tollerance
static double tB = 0;

// Constants for D4T transition
static double Death4Treat = 0;

// Constants for Inflam transition
static double Inflammation = 0;
static double DAMAGEmax = 0;

// Constants for DeathBac transition
static double half_life = 0;

// Constants for Dup transition
static double rCDdup = 0;

// Constants for Efflux transition
static double efflux = 0;

// Constants for Starv transition
static double RCD = 0;

// error difference variation FBAplaces
double eps = 1e-06;

double EX_B_starv = 0;

static map <string, double> ValuePrev{{"BiomassCD", -1}, {"pheme_c", -1}, 
                                      {"cys_L_e", -1}, {"trp_L_e", -1},
                                      {"val_L_e", -1}, {"ile_L_e", -1},
                                      {"leu_L_e", -1}, {"pro_L_e", -1}};

void read_map_string_int(string fname, unordered_map<string,int>& m) {
  ifstream f (fname);
  string line;
  if(f.is_open()) {
    cout << "#### " << fname << "####" << endl;
    int j = 1;
    while (getline(f,line))
    {
      line.erase(remove( line.begin(), line.end(), '\"' ),line.end());
      m.insert(pair<string,int>(line,j));
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

void read_map_string_double(string fname, unordered_map<string,double>& m) {
  ifstream f (fname);
  string line;
  if(f.is_open()) {
    size_t pos = 0, length = 0;
    cout << "#### " << fname << "####" << endl;
    int j = 1;
    while (getline(f,line)) {
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
  else {
    std::cerr<<"\nUnable to open " << fname << ": file do not exists\n";
    exit(EXIT_FAILURE);
  }
}

void read_constant(string fname, double& Infection_rate) {
  ifstream f (fname);
  string line;
  if(f.is_open()) {
    int i = 1;
    while (getline(f,line)) {
      switch(i) {
      case 1:
        Infection_rate = stod(line);
        //cout << "p" << i << ": " << line << "\t" << p1 << endl;
        break;
      }
      ++i;
    }
    f.close();
  }
  else {
    std::cerr<<"\nUnable to open " << fname << ": file do not exists\": file do not exists\n";
    exit(EXIT_FAILURE);
  }
}

void init_data_structures(const struct InfTr* Trans, map <string,int>& NumTrans) {
  
  // constants used:
  
  read_constant("./gDW_CDmin", gDW_CDmin);
  read_constant("./gDW_CDmean", gDW_CDmean);
  read_constant("./gDW_CDmax", gDW_CDmax);
  read_constant("./MWbio", MWbio);
  read_constant("./nBacMax", nBacMax);
  read_constant("./Inflammation", Inflammation);
  read_constant("./DAMAGEmax", DAMAGEmax);
  read_constant("./half_life", half_life);
  read_constant("./rCDdup", rCDdup);
  read_constant("./Efflux", efflux);
  read_constant("./RCD", RCD);
  read_constant("./Death4Treat", Death4Treat);
  read_constant("./Na", Na);
  read_constant("./c", c);
  read_constant("./P", P);
  read_constant("./tB", tB);
  read_constant("./EX_B_starv", EX_B_starv);
  
  // maps used:
  FBAmet["EX_biomass_e_in"] = "EX_biomass_e";
  FBAmet["sink_pheme_c_in"] = "sink_pheme_c";
  FBAmet["EX_cys_L_e_in"] = "EX_cys_L_e";
  FBAmet["EX_trp_L_e_in"] = "EX_trp_L_e";
  FBAmet["EX_val_L_e_in"] = "EX_val_L_e";
  FBAmet["EX_ile_L_e_in"] = "EX_ile_L_e";
  FBAmet["EX_leu_L_e_in"] = "EX_leu_L_e";
  FBAmet["EX_pro_L_e_in"] = "EX_pro_L_e";
  
  // Biomass was normalized
  FBAvars["EX_biomass_e"] =  0.2975482;
  FBAvars["EX_cys_L_e"] = - 0.0000000003436271;
  FBAvars["EX_ile_L_e"] =  - 0.000000006328174;
  FBAvars["EX_leu_L_e" ] = 0.0;
  FBAvars["EX_pro_L_e" ] = 0.0;
  FBAvars["EX_trp_L_e" ] = - 0.0000000004240291;
  FBAvars["EX_val_L_e" ] = 0.0;
  FBAvars["sink_pheme_c"]= - 0.000001172481;
  
  FBAplace["EX_biomass_e_in"] = "BiomassCD";
  FBAplace["sink_pheme_c_in"] = "pheme_c";
  FBAplace["EX_cys_L_e_in"] = "cys_L_e";
  FBAplace["EX_trp_L_e_in"] = "trp_L_e";
  FBAplace["EX_val_L_e_in"] = "val_L_e";
  FBAplace["EX_ile_L_e_in"] = "ile_L_e";
  FBAplace["EX_leu_L_e_in"] = "leu_L_e";
  FBAplace["EX_pro_L_e_in"] = "pro_L_e";
  
  Flag = 1;
  
}

double trunc(double value, double decimal) {
  const double multiplier = std::pow(10.0, decimal);
  return std::floor(value * multiplier) / multiplier;
  //return ((unsigned long int)(value * multiplier)) / multiplier;
}

double FBA(double *Value,
           vector<class FBGLPK::LPprob>& vec_fluxb,
           map <string,int>& NumTrans,
           map <string,int>& NumPlaces,
           const vector<string> & NameTrans,
           const struct InfTr* Trans,
           const int T,
           const double& time) {
  
  double decimalTrunc = 12;
  
  // Check if the data structures need initialization
  if (Flag == -1)
    init_data_structures(Trans, NumTrans);
  
  // Retrieve values for specific places
  double nBac = trunc(Value[NumPlaces["CD"]], decimalTrunc);
  
  bool Out = false;
  bool In = false;
  
  string str = NameTrans[T];
  
  // Check if the string contains "_out" substring
  if (str.find("_out") != string::npos) {
    // Replace "_out" with "_in" in the string
    str.replace(str.find("_out"), 4, "_in");
    Out = true;
  } else {
    In = true;
  }
  
  // Retrieve biomass value
  double Biom = trunc(Value[NumPlaces["BiomassCD"]], decimalTrunc);
  
  // Retrieve rateFBA value
  double rateFBA = FBAvars[FBAmet[str]];
  double rate = 0;
  
  if (str == "EX_biomass_e_in") {
    // Calculate the rate for "EX_biomass_e_in" case
    rate = MWbio * rateFBA * Biom * (1 - (Biom / gDW_CDmax));
  } else {
    double met = Value[NumPlaces[FBAplace[str]]];
    
    rate = rateFBA * met;
    
    // rateFBA = solution (mmol/g/h)
    // ratePT = Transition parameters (Cmolecules/h)
    // Na = (molecules/mmol)
    // c = (molecules)
    // nBac = bacterial cells number (#)
    // Biom = PT biomassCD place (pg = g*1e-12)
    
    // Calculate the rate for "_L_e" case
    // ratePT = (rateFBA * Na * (1 / c)) * (nBac * Biom * 1e-12) = (C_molecules) 
    // ((Met * cNa) / (nBac * Biom * 1e-12))
    
    // Calculate the rate for other cases
    double r = (rate * (str.find("_L_e") != string::npos ? Na * (1 / met) : 1e+09)) * (nBac * Biom * 1e-12);
    
    // Set the rate based on the conditions
    if ((Out && rateFBA > 0) || (In && rateFBA < 0))
      rate = r;
    else if ((Out && rateFBA < 0) || (In && rateFBA > 0))
      rate = 0;
  }
  
  return rate;
}


// Inflam transition
double Heam(double *Value,
            vector<class FBGLPK::LPprob>& vec_fluxb,
            map<string, int>& NumTrans,
            map<string, int>& NumPlaces,
            const vector<string>& NameTrans,
            const struct InfTr* Trans,
            const int T,
            const double& time) {
  
  // Check if the data structures need initialization
  if (Flag == -1)
    init_data_structures(Trans, NumTrans);
  
  // Retrieve the value of the Damage place and calculate the percentage of damage
  double DamagePlace = Value[NumPlaces["Damage"]];
  double PercDamage = DamagePlace / DAMAGEmax;
  
  // Calculate the rate based on the percentage of damage
  double rate = PercDamage * Inflammation;
  
  // Print the rate if FlagDebug is set to 1
  if (FlagDebug == 1) {
    cout << "Inflammation (pmol): " << rate << endl;
  }
  
  return rate;
}

// DeathBac transition
double DeathCD(double *Value,
               vector<class FBGLPK::LPprob>& vec_fluxb,
               map<string, int>& NumTrans,
               map<string, int>& NumPlaces,
               const vector<string>& NameTrans,
               const struct InfTr* Trans,
               const int T,
               const double& time) {
  
  if(Flag == -1) init_data_structures(Trans, NumTrans);
  
  double rate = 0;
  
  double Biom = Value[NumPlaces.find("BiomassCD") -> second];
  double nBac = Value[NumPlaces.find("CD") -> second];
  
  double k = 5;
  
  rate = half_life*nBac*(1/(2 + exp(k*(Biom - gDW_CDmean))));
  
  if(FlagDebug == 1) {
    cout<< "half_life: " << half_life << endl;
  }
  
  return(rate);
  
}

// Dup transition
double Duplication(double *Value,
                   vector<class FBGLPK::LPprob>& vec_fluxb,
                   map<string, int>& NumTrans,
                   map<string, int>& NumPlaces,
                   const vector<string>& NameTrans,
                   const struct InfTr* Trans,
                   const int T,
                   const double& time) {
  
  if(Flag == -1) init_data_structures(Trans, NumTrans);
  
  double rate = 0.0;
  
  double Biom = Value[NumPlaces.find("BiomassCD") -> second];
  double nBac = Value[NumPlaces.find("CD") -> second];
  
  if(FlagDebug == 1) {
    cout<< "Biom (pg): " << Biom << endl;
    cout<< "nBac (cell): " << nBac << endl;
  }
  
  // if((Biom - gDW_CDmin) > tB){
  //   rate = Biom*nBac*rCDdup*(1 - (nBac/nBacMax));
  // }
  rate = ((Biom - gDW_CDmin)/(gDW_CDmax - gDW_CDmin)) *nBac* rCDdup* (1 - (nBac/nBacMax));
  
  
  if(FlagDebug == 1) {
    cout << "Biom (pg) - gDW_CDmin (pg): " << Biom - gDW_CDmin << endl;
    cout << "Dup rate (cell): " << rate << endl;
  }
  
  return(rate);
}

// Starv transition
double Starvation(double *Value,
                  vector<class FBGLPK::LPprob>& vec_fluxb,
                  map<string, int>& NumTrans,
                  map<string, int>& NumPlaces,
                  const vector<string>& NameTrans,
                  const struct InfTr* Trans,
                  const int T,
                  const double& time) {
  
  double rate = 0;
  
  if(Flag == -1) {
    init_data_structures(Trans, NumTrans);
  }
  
  double Biom = Value[NumPlaces.find("BiomassCD") -> second];
  
  rate = (EX_B_starv)*(Biom)*MWbio;
  
  if(FlagDebug == 1) {
    cout << "EX_B_starv (mmol/gDW*h): " << EX_B_starv << ";" << endl;
    cout << "starv rate (pg): " << rate << ";" << endl;
  }
  
  return(rate);
  
}

// Death4Treat
double DfourT(double *Value,
              vector<class FBGLPK::LPprob>& vec_fluxb,
              map<string, int>& NumTrans,
              map<string, int>& NumPlaces,
              const vector<string>& NameTrans,
              const struct InfTr* Trans,
              const int T,
              const double& time) {
  
  double rate = 0;
  
  double nBac = Value[NumPlaces.find("CD") -> second];
  double DrugPlace = Value[NumPlaces.find("Drug") -> second];
  
  //   MM = 171.16*1e-12 # (g/pmol)
  //
  // // ref (1) (doi:_?)
  //   MIC.conc = 1e-06 # (g/mL)
  //   MIC.conc = (MIC.conc/MM)*1e-09 # (mmol)
  //   MIC.conc = MIC.conc*1e+09 # (pmol)
  //   MIC.conc = MIC.conc*602214150000 # (molecule)
  //
  // // dose = 0.5 # (g/day) - prescribed dose -
  //   dose = 0.05 # (g/day) - MIC dose -
  // // dose = 0.01 # (g/day) - SubMIC dose -
  
  double DrugPlaceMIC = 5e+03;
  
  // testing.function = function (DrugPlace, nBac, Death4Treat, DrugPlaceMIC) {
  //       DrugPlace*nBac*Death4Treat*(exp((DrugPlace/DrugPlaceMIC)))
  // }
  //
  // plot(testing.function(c(0, 1e+3, 1e+5, 3e+15, 6e+15, 1e+20, 1e+50), 1e+09, 1e-12, 3e+15))
  
  // rate = DrugPlace*nBac*Death4Treat*(exp(((DrugPlace/DrugPlaceMIC) - 1)));
  rate = DrugPlace*nBac*Death4Treat*(exp((DrugPlace/DrugPlaceMIC)));
  
  if(FlagDebug == 1) {
    cout << "Death4Treat: " << Death4Treat << endl;
    cout << "nBac = " << nBac << " (cell)" << endl;
    cout << "DrugPlace = " << DrugPlace << " (pmol)" << endl;
    cout << "rate Death4Treat = " << rate << " (cell)" << endl;
  }
  
  return(rate);
  
}

// Efflux transition
double Efflux(double *Value,
              vector<class FBGLPK::LPprob>& vec_fluxb,
              map<string, int>& NumTrans,
              map<string, int>& NumPlaces,
              const vector<string>& NameTrans,
              const struct InfTr* Trans,
              const int T,
              const double& time) {

  double nBac = Value[NumPlaces.find("CD") -> second];
  double DrugPlace = Value[NumPlaces.find("Drug") -> second];
  
  double rate = DrugPlace*nBac*efflux;

  
  return(rate);
  
}
