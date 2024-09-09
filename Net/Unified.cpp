
#include <string.h>
#include <sys/resource.h>
#include <map>
#include <regex>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <string>
#include<sstream>

using namespace FBGLPK;

static double* Vars;

static double FBAtime = -1;

static map <string, string> FBAmet;
static map <string, string> FBAplace;

static double Flag = -1;
static bool FlagPrev = 0;

static bool FlagDebug = 1;

double rate = 0;

static double gDW_CDmean = 0;
static double gDW_CDmin = 0;
static double gDW_CDmax = 0;

// https://doi.org/10.1016/j.ymben.2021.10.012
// biomass molecular weigth (g/mmol)
// biomass defined to have a molecular weight (MW) of 1 (g/mmol)
static double MWbio = 0;

static double nBacMax = 0;

// chemical proportionality factors
static double Na = 0;
static double c = 0;
static double pack = 0;
// constant for sensitivity analysis in FBA
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

// error difference variation FBAplaces
double eps = 1e-06;

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
      cout << line.substr(0,pos) << ": " << stod(line.substr(pos+1,length)) << " ";
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
  read_constant("./Death4Treat", Death4Treat);
  read_constant("./Na", Na);
  read_constant("./c", c);
  read_constant("./P", P);
  read_constant("./tB", tB);
  
  pack = 1*(Na*(1/c));
  
  FBAmet["EX_biomass_e_in"] = "EX_biomass_e";
  FBAmet["sink_pheme_c_in"] = "sink_pheme_c";
  FBAmet["EX_cys_L_e_in"] = "EX_cys_L_e";
  FBAmet["EX_trp_L_e_in"] = "EX_trp_L_e";
  FBAmet["EX_val_L_e_in"] = "EX_val_L_e";
  FBAmet["EX_ile_L_e_in"] = "EX_ile_L_e";
  FBAmet["EX_leu_L_e_in"] = "EX_leu_L_e";
  FBAmet["EX_pro_L_e_in"] = "EX_pro_L_e";
  
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
  
  bool FBAmarking = 0;
  // truncation floating
  double decimalTrunc = 12;
  
  if(Flag == -1) init_data_structures(Trans, NumTrans);
  
  double nBac = trunc(Value[NumPlaces.find("CD") -> second], decimalTrunc );

  if (ValuePrev["BiomassCD"] != -1) {
    FlagPrev = 1;
  } else {
    FBAmarking = 1;
  }
  
  if (FlagPrev) {
    
    // cout << "err evaluating ..." << endl;
    
    map <string, double>::iterator p = ValuePrev.begin();
    
    while(p != ValuePrev.end() && !FBAmarking) {
      
      double MPrev = p -> second;
      double M = trunc(Value[NumPlaces.find(p -> first) -> second], decimalTrunc );
      
      double err = 0;
      
      if (MPrev == 0.0) {
        err = abs(M - MPrev);
      } else {
        err = (abs(M - MPrev)/MPrev);
      }
      
      // cout << "err: " << err << endl;
      
      if (err > eps) {
        // cout << "FBA must be call" << endl;
        FBAmarking = 1;
      }
      ++p;
    }
  }
  
  if(FBAmarking && FBAtime != time) {
    
    // cout << "FBA calling ..." << endl;
    
    for(map <string, double>::iterator p = ValuePrev.begin(); p != ValuePrev.end(); ++p) {
      p -> second = trunc(Value[NumPlaces.find(p -> first) -> second], decimalTrunc );
    }
    
    for (map<string, string>::iterator p = FBAmet.begin(); p != FBAmet.end(); ++p) {
      
      int index = vec_fluxb[0].fromNametoid(p->second);
      string TypeBound = "GLP_DB";
      
      double Ub = vec_fluxb[0].getUpBounds(index);
      double Lb = vec_fluxb[0].getLwBounds(index);
      
      double Biom = trunc(Value[NumPlaces.find("BiomassCD") -> second], decimalTrunc );
      
      if(p -> first == "EX_biomass_e_in") {
        
        double required = 0.14;
        double max = 0.2;
        double stiff = 150;

        if (nBac < 1) {
          
          Ub = Biom*nBac;
          
        } else {
          
          if ((gDW_CDmax - Biom) > tB) {
            Ub = (gDW_CDmax - Biom);
          } else {
            Ub = tB;
          }
          
          Ub = Ub > max ? max : Ub;
          
          if (Biom < gDW_CDmin) {
            Ub = 0.0;
          } else if (Biom == gDW_CDmin) {
            Ub = max;
          }
          
          Ub = Ub / (1 + exp(-stiff * (Biom - gDW_CDmin)));
          Ub = Ub * (1 - 1 / (1 + exp(-stiff * (Biom - required))));
        }
        
      }
      
      else if (p -> first == "sink_pheme_c_in") {
        
        double Met = trunc(Value[NumPlaces[FBAplace[p -> first]]], decimalTrunc );
        double E = 1e-24; // (pg)
        
        if (nBac < 1) {
          Lb = - (Met*1e-09)/((Biom*1e-12) + E);
        } else {
          Lb = - (Met*1e-09)/((nBac*(Biom*1e-12)) + E);
        }
        
        if (Biom < gDW_CDmin) {
          Lb = Lb * (Biom/gDW_CDmin);
        }
        
        if (Lb < P) {
          Lb = P;
        }
      }
      
      else {
        
        double Met = trunc(Value[NumPlaces[FBAplace[p->first]]], decimalTrunc );
        double E = 1e-24; // (pg)
        
        if (nBac < 1) {
          Lb = - (Met/pack)/((Biom*1e-12) + E);
        } else {
          Lb = - (Met/pack)/((nBac*(Biom*1e-12)) + E);
        }
        
        if (Biom < gDW_CDmin) {
          Lb = Lb * (Biom/gDW_CDmin);
        }
        
        if (Lb < P) {
          Lb = P;
        }
        
        Ub = 0.0;
      }
      
      double Lbtr = trunc( Lb, decimalTrunc );
      double Ubtr = trunc( Ub, decimalTrunc );
      
      vec_fluxb[0].update_bound(index, TypeBound, Lbtr, Ubtr);
      
      if (FlagDebug == 1) {
        cout << "Transition:" << p -> first << endl;
        cout << "bounds: [" << Lbtr << " , " << Ubtr << "]" << endl;
      }
      
    }
    
    vec_fluxb[0].solve();
    Vars = vec_fluxb[0].getVariables();
    
    FBAtime = time;
    
  }
  
  bool Out = 0;
  bool In = 0;
  
  string str = NameTrans[T];
  
  if(str.find("_out") != string::npos){
    str = std::regex_replace(str, std::regex("_out"), "_in");
    Out = 1;
  } else {
    In = 1;
  }
  
  int index = vec_fluxb[0].fromNametoid(FBAmet.find(str) -> second);

  double Biom = trunc(Value[NumPlaces.find("BiomassCD") -> second], decimalTrunc );
  
  rate = trunc(Vars[index], decimalTrunc );
  double r = 0;
  
  if (str == "EX_biomass_e_in") {
    r = Biom*MWbio*rate;

  } else if(str.find("_L_e") != string::npos) {
    r = (rate*pack)*(nBac*Biom*1e-12);
    
  } else {
    r = ((rate*1e+09)*(nBac*Biom*1e-12));
  }
  
  if((Out) && (rate > 0))
    rate = r;
  else if((Out) && (rate < 0))
    rate = 0;
  else if((In) && (rate > 0))
    rate = 0;
  else if((In) && (rate < 0))
    rate = -r;
  
  return(rate);
  
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
    cout << "Inflammation: " << rate << endl;
  }
  
  return rate;
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
  
  if(Flag == -1) init_data_structures(Trans, NumTrans);
  
  double rate = 0;
  
  double Biom = Value[NumPlaces.find("BiomassCD") -> second];
  double nBac = Value[NumPlaces.find("CD") -> second];
  
  double k = 100;
  
  rate = half_life*nBac*(1/(2 + exp(k*(Biom - gDW_CDmean))));
  
  if(FlagDebug == 1) {
    cout << "Biom (pg): " << Biom << endl;
    cout << "nBac (cell): " << nBac << endl;
    cout<< "DeathCD: " << rate << endl;
  }
  
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
  
  if(Flag == -1) init_data_structures(Trans, NumTrans);
  
  double Biom = Value[NumPlaces.find("BiomassCD") -> second];
  double nBac = Value[NumPlaces.find("CD") -> second];
  
  double d = (Biom - gDW_CDmin) / (gDW_CDmax - gDW_CDmin);
  double cap = 1 - (nBac / nBacMax);
  double exp_factor = exp(d * cap);
  double rate = d * nBac * rCDdup * cap * exp_factor;
  
  if(FlagDebug == 1) {
    cout << "Biom (pg): " << Biom << endl;
    cout << "nBac (cell): " << nBac << endl;
    cout << "Duplication rate: " << rate << endl;
  }
  
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
  
  if(Flag == -1) {
    init_data_structures(Trans, NumTrans);
  }
  
  double Biom = Value[NumPlaces.find("BiomassCD") -> second];
  
  double EX_B_starv = 0.0080178; // (pg/h)
  double max_rate = 0.011454; // (pg/h)
  
  double rate;
  
  double biom_gDW_ratio = Biom / gDW_CDmax;
  rate = max_rate * pow(biom_gDW_ratio, 3) / (pow(biom_gDW_ratio, 3) + max_rate);
  
  if(FlagDebug == 1) {
    cout << "Biom (pg): " << Biom << endl;
    cout << "Starvation rate: " << rate << endl;
  }
  
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
  
  double rate = 0;
  
  double nBac = Value[NumPlaces.find("CD") -> second];
  double DrugPlace = Value[NumPlaces.find("Drug") -> second];
  
  double DrugPlaceMIC = 5842.487;

  rate = DrugPlace*nBac*Death4Treat*(exp((DrugPlace/DrugPlaceMIC)));
  
  if(FlagDebug == 1) {
    cout << "nBac (cell): " << nBac << endl;
    cout << "DrugPlace (pmol): " << DrugPlace << endl;
    cout << "DfourT rate: " << rate << endl;
  }
  
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
  
  double nBac = Value[NumPlaces.find("CD") -> second];
  double DrugPlace = Value[NumPlaces.find("Drug") -> second];
  
  double rate = DrugPlace*nBac*efflux;
  
  if(FlagDebug == 1) {
    cout << "nBac (cell): " << nBac << endl;
    cout << "DrugPlace (pmol): " << DrugPlace << endl;
    cout << "Efflux rate: " << rate << endl;
  }
  
  return(rate);
  
}
