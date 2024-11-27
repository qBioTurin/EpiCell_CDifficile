
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

// param toxin production
static double kT = 0;

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
  read_constant("./kT", kT);
  read_constant("./EX_B_starv", EX_B_starv);
  
  // maps used:
  
  FBAmet["EX_biomass_e_in"] = "EX_biomass(e)";
  FBAmet["sink_pheme_c_in"] = "sink_pheme(c)";
  FBAmet["EX_cys_L_e_in"] = "EX_cys_L(e)";
  FBAmet["EX_trp_L_e_in"] = "EX_trp_L(e)";
  FBAmet["EX_val_L_e_in"] = "EX_val_L(e)";
  FBAmet["EX_ile_L_e_in"] = "EX_ile_L(e)";
  FBAmet["EX_leu_L_e_in"] = "EX_leu_L(e)";
  FBAmet["EX_pro_L_e_in"] = "EX_pro_L(e)";
  
  FBAvars["EX_biomass(e)"] =  0.2975482; // ho normalizzato anche la biomassa
  FBAvars["EX_cys_L(e)"] = - 0.0000000003436271;
  FBAvars["EX_ile_L(e)"] =  - 0.000000006328174;
  FBAvars["EX_leu_L(e)" ] = 0.0;
  FBAvars["EX_pro_L(e)" ] = 0.0;
  FBAvars["EX_trp_L(e)" ] = - 0.0000000004240291;
  FBAvars["EX_val_L(e)" ] = 0.0;
  FBAvars["sink_pheme(c)"]= - 0.000001172481;
  
  
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
           map <string,int>& NumTrans,
           map <string,int>& NumPlaces,
           const vector<string> & NameTrans,
           const struct InfTr* Trans,
           const int T,
           const double& time) {
  
  double rate = 0;
  double rateFBA = 0;
  bool FBAmarking = 0;
  
  // truncation floating
  double decimalTrunc = 12;
  
  if(Flag == -1) init_data_structures(Trans, NumTrans);
  
  double nBac = trunc(Value[NumPlaces.find("CD") -> second], decimalTrunc);
  double DrugPlace = trunc(Value[NumPlaces.find("Drug") -> second], decimalTrunc);
  
  bool Out = 0;
  bool In = 0;

  string str = NameTrans[T];
  
  if(str.find("_out") != string::npos){
    str = std::regex_replace(str, std::regex("_out"), "_in");
    Out = 1;
  } else {
    In = 1;
  }
  
  double Biom = trunc(Value[NumPlaces.find("BiomassCD") -> second], decimalTrunc );
  
  rateFBA = FBAvars.find(FBAmet.find(str) -> second) -> second ;
  double r = 0;
  
  if (str == "EX_biomass_e_in") {
  
    r = MWbio*rateFBA*Biom*(1 - (Biom/gDW_CDmax));

  } else if(str.find("_L_e") != string::npos) {
    
    double met = Value[NumPlaces.find(FBAplace.find(str) -> second) -> second];
    
    rate = rateFBA*met;
    // rateFBA = solution (mmol/g/h)
    // ratePT = Transition parameters (Cmolecules/h)
    // Na = (molecules/mmol)
    // c = (molecules)
    // nBac = bacterial cells number (#)
    // Biom = PT biomassCD place (pg = g*1e-12)
    
    // trans_in when is neg, otherwise trans_out
    // ratePT = (rateFBA*Na*(1/c))*(nBac*Biom*1e-12) = (C_molecules) 
    // ((Met*cNa)/(nBac*Biom*1e-12))
    
    r = (rate*Na*(1/c))*(nBac*Biom*1e-12);

  } else {
    
    double met = Value[NumPlaces.find(FBAplace.find(str) -> second) -> second];
    
    rate = rateFBA*met;
    r = (rate*1e+09)*(nBac*Biom*1e-12);
   
  }
  
  if((Out) && (rateFBA > 0))
    rate = r ;
  else if((Out) && (rateFBA < 0))
    rate = 0;
  else if((In) && (rateFBA > 0))
    rate = 0;
  else if((In) && (rateFBA < 0))
    rate = -r;
  
  return(rate);
  
}


// Inflam transition
double Heam(double *Value,
            map <string,int>& NumTrans,
            map <string,int>& NumPlaces,
            const vector<string> & NameTrans,
            const struct InfTr* Trans,
            const int T,
            const double& time) {
  
  if(Flag == -1) init_data_structures(Trans, NumTrans);
  
  double rate = 0;
  double g = 0;
  
  double DamagePlace = Value[NumPlaces.find("Damage") -> second];
  double PercDamage = DamagePlace/DAMAGEmax;
  
  // if(PercDamage <= 0.1){
  //   g = 0;
  // } else if((PercDamage > 0.1) && (PercDamage <= 0.7)){
  //   g = 1/3;
  // } else if(PercDamage > 0.7){
  //   g = 1;
  // }
  
  rate = PercDamage*Inflammation;
  
  if(FlagDebug == 1) {
    cout << "Inflammation (pmol): " << rate << endl;
  }
  
  return(rate);
}

// DeathBac transition
double DeathCD(double *Value,
               map <string,int>& NumTrans,
               map <string,int>& NumPlaces,
               const vector<string> & NameTrans,
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
                   map <string,int>& NumTrans,
                   map <string,int>& NumPlaces,
                   const vector<string> & NameTrans,
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
                  map <string,int>& NumTrans,
                  map <string,int>& NumPlaces,
                  const vector<string> & NameTrans,
                  const struct InfTr* Trans,
                  const int T,
                  const double & time) {
  
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
              map <string,int>& NumTrans,
              map <string,int>& NumPlaces,
              const vector<string> & NameTrans,
              const struct InfTr* Trans,
              const int T,
              const double & time) {
  
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
  double DrugDose = 1e+04;
  
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
              map <string,int>& NumTrans,
              map <string,int>& NumPlaces,
              const vector<string> & NameTrans,
              const struct InfTr* Trans,
              const int T,
              const double& time) {

  double nBac = Value[NumPlaces.find("CD") -> second];
  double DrugPlace = Value[NumPlaces.find("Drug") -> second];
  
  double rate = DrugPlace*nBac*efflux;

  
  return(rate);
  
}

// ToxinP
double ToxinP(double *Value,
              vector<class FBGLPK::LPprob>& vec_fluxb,
              map <string,int>& NumTrans,
              map <string,int>& NumPlaces,
              const vector<string> & NameTrans,
              const struct InfTr* Trans,
              const int T,
              const double & time) {
  
  double rate = 0;
  
  double nBac = Value[NumPlaces.find("CD") -> second];
  double CysPlace = Value[NumPlaces.find("cys_L_e") -> second];
  
  // cys_L_e = 0.25 (μmol/mL)
  // cys_L_e = (6.0221415e+17*(1/6.022e8)) (C_molecules/mL)
  
  // Tmax = 0.001 (pM) -> (mmol/ml) -> (pmol/ml)
  
  double Tmax = 0.001;
  
  // from: https://doi.org/10.3892/ol.2018.7921
  // 300 (ng/ml) of toxin increases human cells apoptotic rates
  // TcdA_MM = 308 (g/mmol) -> 1e+09*308 (ng/mmol)
  // 1/((1/100)*(1e+09*308)) = 3.25e-10 (mmol) -> 3.25e-10*1e+9 = 0.325 (pmol)
  
  // testing.function = function(CD, cys_L_e, Tmax, kT) {
  //   ((CD*Tmax)/(kT + exp(cys_L_e)))
  // }
  
  // plot(testing.function(CD = 45700000,
  //                   cys_L_e = c(1000023497, 100002349, 10000234,
  //                               1000023, 10000, 1000, 100, 50, 25, 12, 10, 1, 0),
  //                               Tmax = 0.001, kT = 1e+05))
  
  rate = ((nBac*Tmax)/(kT + exp(CysPlace)));
  
  if(FlagDebug == 1) {
    cout << "ToxinP rate: " << rate <<  "(pmol/h)" << endl;
    cout << "CysPlace = " << (CysPlace*6.022e+8)/6.0221415e+17 << " (mmol/mL)" << endl;
  }
  
  return(rate);
  
}
