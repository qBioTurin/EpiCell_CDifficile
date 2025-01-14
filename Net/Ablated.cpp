
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
static double MWbio = 0;
static double nBacMax = 0;
static double Na = 0;
static double c = 0;
static double P = 0;
static double tB = 0;
static double Death4Treat = 0;

// Constants for Inflam transition
static double Inflammation = 0;
static double DAMAGEmax = 0;
static double half_life = 0;
static double rCDdup = 0;
static double efflux = 0;
static double RCD = 0;

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
  
  FBAvars["EX_biomass_e"] = 0.297548222125;
  FBAvars["EX_cys_L_e"] = 0.0;
  FBAvars["EX_ile_L_e"] = -4.4089e-08;
  FBAvars["EX_leu_L_e" ] = 0.0;
  FBAvars["EX_pro_L_e" ] = -4.4089e-08;
  FBAvars["EX_trp_L_e" ] = -4.25e-10;
  FBAvars["EX_val_L_e" ] = -2.2879e-08;
  FBAvars["sink_pheme_c"]= -1.172481e-06;
  
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
  
  if (Flag == -1)
    init_data_structures(Trans, NumTrans);
  
  bool Out = 0;
  bool In = 0;
  
  string str = NameTrans[T];
  
  if(str.find("_out") != string::npos){
    str = regex_replace(str, regex("_out"), "_in");
    Out = 1;
  } else {
    In = 1;
  }
  
  double nBac = trunc(Value[NumPlaces["CD"]], decimalTrunc);
  double Biom = trunc(Value[NumPlaces["BiomassCD"]], decimalTrunc);
  
  cout << "str = " << str << endl;
  
  double rateFBA = trunc(FBAvars[FBAmet[str]], decimalTrunc);
  
  double r = 0;
  double rate = 0;
  
  if (str == "EX_biomass_e_in") {
    
    r = MWbio*rateFBA*Biom*(1 - (Biom/gDW_CDmax));
    
    cout << "Biom = " << Biom << endl;
    cout << "rateFBA = " << rateFBA << endl;
    cout << "----------------------------------------" << endl;
    
  } else if(str.find("_L_e") != string::npos) {
    
    double met = trunc(Value[NumPlaces[FBAplace[str]]], decimalTrunc);
    
    r = rateFBA*met*1e+09*(nBac*Biom * 1e-12);
    
    cout << "met = " << met << endl;
    cout << "rateFBA = " << rateFBA << endl;
    cout << "rateFBA*met = " << (rateFBA*met) << endl;
    cout << "----------------------------------------" << endl;
    
  } else{
    
    double met = trunc(Value[NumPlaces[FBAplace[str]]], decimalTrunc );
    
    r = (rateFBA*met*Na*(1/c))*(nBac*Biom*1e-12); 
    
    cout << "met = " << met << endl;
    cout << "rateFBA = " << rateFBA << endl;
    cout << "rateFBA*met = " << (rateFBA*met) << endl;
    cout << "----------------------------------------" << endl;
    
  }
  
  r = trunc(r, decimalTrunc);
  
  cout << "r = " << r << " | " << str << " | " << "time = "<< time << endl;
  cout << "----------------------------------------" << endl;
  
  if((Out) && (rateFBA > 0))
    rate = r;
  else if((Out) && (rateFBA < 0))
    rate = 0;
  else if((In) && (rateFBA > 0))
    rate = 0;
  else if((In) && (rateFBA < 0))
    rate = -r;
  
  return rate;
  
}

// Inflam transition
double Heam(double *Value,
            vector<class FBGLPK::LPprob>& vec_fluxb,
            map <string,int>& NumTrans,
            map <string,int>& NumPlaces,
            const vector<string> & NameTrans,
            const struct InfTr* Trans,
            const int T,
            const double& time) {
  
  double decimalTrunc = 12;
  
  if(Flag == -1) init_data_structures(Trans, NumTrans);
  
  double rate = 0;
  
  double DamagePlace = Value[NumPlaces.find("Damage") -> second];
  double PercDamage = DamagePlace/DAMAGEmax;
  
  rate = PercDamage*Inflammation;
  rate = trunc(rate, decimalTrunc);
  
  return(rate);
}

// DeathBac transition
double DeathCD(double *Value,
               vector<class FBGLPK::LPprob>& vec_fluxb,
               map <string,int>& NumTrans,
               map <string,int>& NumPlaces,
               const vector<string> & NameTrans,
               const struct InfTr* Trans,
               const int T,
               const double& time) {
  
  double decimalTrunc = 12;
  
  if(Flag == -1) init_data_structures(Trans, NumTrans);
  
  double rate = 0;
  
  double Biom = Value[NumPlaces.find("BiomassCD") -> second];
  double nBac = Value[NumPlaces.find("CD") -> second];
  
  double k = 5;
  
  rate = half_life*nBac*(1/(2 + exp(k*(Biom - gDW_CDmean))));
  rate = trunc(rate, decimalTrunc);
  
  return(rate);
  
}

// Dup transition
double Duplication(double *Value,
                   vector<class FBGLPK::LPprob>& vec_fluxb,
                   map <string,int>& NumTrans,
                   map <string,int>& NumPlaces,
                   const vector<string> & NameTrans,
                   const struct InfTr* Trans,
                   const int T,
                   const double& time) {
  
  double decimalTrunc = 12;
  
  if(Flag == -1) init_data_structures(Trans, NumTrans);
  
  double rate = 0.0;
  
  double Biom = Value[NumPlaces.find("BiomassCD") -> second];
  double nBac = Value[NumPlaces.find("CD") -> second];
  
  if(FlagDebug == 1) {
  }
  
  rate = ((Biom - gDW_CDmin)/(gDW_CDmax - gDW_CDmin)) *nBac* rCDdup* (1 - (nBac/nBacMax));
  rate = trunc(rate, decimalTrunc);
  
  return(rate);
}

// Starv transition
double Starvation(double *Value,
                  vector<class FBGLPK::LPprob>& vec_fluxb,
                  map <string,int>& NumTrans,
                  map <string,int>& NumPlaces,
                  const vector<string> & NameTrans,
                  const struct InfTr* Trans,
                  const int T,
                  const double & time) {
  
  double decimalTrunc = 12;
  
  if(Flag == -1) {
    init_data_structures(Trans, NumTrans);
  }
  
  double Biom = Value[NumPlaces.find("BiomassCD") -> second];
  
  double rate = (EX_B_starv)*(Biom)*MWbio;
  rate = trunc(rate, decimalTrunc);
  
  return(rate);
  
}

// Death4Treat
double DfourT(double *Value,
              vector<class FBGLPK::LPprob>& vec_fluxb,
              map <string,int>& NumTrans,
              map <string,int>& NumPlaces,
              const vector<string> & NameTrans,
              const struct InfTr* Trans,
              const int T,
              const double & time) {
  
  double decimalTrunc = 12;
  
  double rate = 0;
  
  double nBac = Value[NumPlaces.find("CD") -> second];
  double DrugPlace = Value[NumPlaces.find("Drug") -> second];
  
  double DrugPlaceMIC = 5e+03;
  
  rate = DrugPlace*nBac*Death4Treat*(exp((DrugPlace/DrugPlaceMIC)));
  rate = trunc(rate, decimalTrunc);
  
  return(rate);
  
}

// Efflux transition
double Efflux(double *Value,
              vector<class FBGLPK::LPprob>& vec_fluxb,
              map <string,int>& NumTrans,
              map <string,int>& NumPlaces,
              const vector<string> & NameTrans,
              const struct InfTr* Trans,
              const int T,
              const double& time) {
  
  double decimalTrunc = 12;
  
  double nBac = Value[NumPlaces.find("CD") -> second];
  double DrugPlace = Value[NumPlaces.find("Drug") -> second];
  
  double rate = DrugPlace*nBac*efflux;
  rate = trunc(rate, decimalTrunc);
  
  return(rate);
  
}
