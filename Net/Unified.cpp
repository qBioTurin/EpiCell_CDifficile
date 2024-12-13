
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
static double cNa = 0;
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
  read_constant("./RCD", RCD);
  read_constant("./Death4Treat", Death4Treat);
  read_constant("./Na", Na);
  read_constant("./c", c);
  read_constant("./P", P);
  read_constant("./tB", tB);
  read_constant("./EX_B_starv", EX_B_starv);
  read_constant("./eps", eps);
  
  cNa = c/Na;
  
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
  
  if(FBAmarking) {// && FBAtime != time) {
    
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

        if ((gDW_CDmax - Biom) > tB) {
          Ub = (gDW_CDmax - Biom);
        } else {
          Ub = tB;
        }
        
      }
      else if (p -> first == "sink_pheme_c_in") {
        
        double Met = trunc(Value[NumPlaces[FBAplace[p -> first]]], decimalTrunc );
        
        if(nBac < 1){
          Lb = -(Met*1e-09)/(Biom*1e-12);
        } else {
          Lb = (-(Met*1e-09)/(nBac*Biom*1e-12) );
        }
        
        Ub = 10;
      }
      else {
        
        double Met = trunc(Value[NumPlaces[FBAplace[p->first]]], decimalTrunc );
        
        // in the model AAs are measured as (C_molecules)
        // in FBA bounds are in (mmol)
        // quick conversion: 1 mmol = 6.02214154+20 molecules
        // 6.02214154+20 molecules = Na
        
        if(nBac < 1){
          //Lb = -( Met*cNa );
          Lb = -( (Met*cNa)/(Biom*1e-12)  );
        } else {
          Lb = -( (Met*cNa)/(nBac*Biom*1e-12) );
        }
        //}
      
        Ub = 0.0;
        
      }
      
      double Lbtr = trunc(Lb, decimalTrunc );
      double Ubtr = trunc(Ub, decimalTrunc );
      
      // cout << "Transition:" << p -> first << endl;
      // cout << "bounds: [" << Lb <<" , " << Ub << "]" << endl;
      // cout << "Truncated bounds: [" << Lbtr <<" , " << Ubtr << "]" << endl;
      
      vec_fluxb[0].update_bound(index, TypeBound, Lbtr, Ubtr);
      
    }
    
    vec_fluxb[0].solve();
    Vars = vec_fluxb[0].getVariables();
  
    //double marking = 0;
    //double flux = 0;
    
   // for(map <string, string>::iterator p = FBAmet.begin(); p != FBAmet.end(); ++p) {
     //int indexFBA = vec_fluxb[0].fromNametoid(p -> second);
      //flux = trunc(Vars.at(indexFBA), decimalTrunc );
      
      //int indexPN = NumPlaces.find(FBAplace.find(p -> first) -> second) -> second;
      //marking = Value[indexPN];
      
      //if (marking == 0) {
      //  Vars.at(indexFBA) = flux;
      //} else {
      //  Vars.at(indexFBA) = (flux/marking);
      //}
      //}
    
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
  double rate = trunc(Vars[index], decimalTrunc ) ;
  double r = 0;

  if (str == "EX_biomass_e_in") {
      r = MWbio*rate*Biom*(1 - (Biom/gDW_CDmax));
    } else if(str.find("_L_e") != string::npos) {
      r = rate * 1e+09 * (nBac * Biom * 1e-12);
    } else{
      r = (rate*Na*(1/c))*(nBac*Biom*1e-12); 
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
            map <string,int>& NumTrans,
            map <string,int>& NumPlaces,
            const vector<string> & NameTrans,
            const struct InfTr* Trans,
            const int T,
            const double& time) {
  
  if(Flag == -1) init_data_structures(Trans, NumTrans);
  
  double rate = 0;
  
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
  
  double rate = (EX_B_starv)*(Biom)*MWbio;
  
  if(FlagDebug == 1) {
    cout << "EX_B_starv (mmol/gDW*h): " << EX_B_starv << ";" << endl;
    cout << "starv rate (pg): " << rate << ";" << endl;
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
  //double DrugDose = 1e+04;
  
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
    cout << "nBac = " << nBac << " (cell)" << endl;
    cout << "DrugPlace = " << DrugPlace << " (pmol)" << endl;
    cout << "efflux = " << efflux << " (1/pmol*h*cell)" << endl;
    cout << "Drug Efflux = " << rate << "(pmol)" << endl;
  }
  
  return(rate);
  
}
