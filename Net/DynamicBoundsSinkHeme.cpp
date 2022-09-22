
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

//static double FlagFBA = -1; 1 when the FBA was already run at that time!

static double* Vars;

static double FBAtime = -1;

static map <string, string> FBAmet;
static map <string, string> FBAplace;

static double Flag = -1;

double rate = 0;

static double gDW_CDmean = 0;
static double gDW_CDmin = 0;
static double gDW_CDmax = 0;

// https://doi.org/10.1016/j.ymben.2021.10.012
// biomass molecular weigth [g/mmol]
static double MWbio = 0;

static double nBacMax = 0;

// chemical proportionality factors
static double Na = 0;
static double c = 0;

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
  
  read_constant("./Na", Na);
  read_constant("./c", c);
  
  read_constant("./Death4Treat", Death4Treat);
  
  FBAmet["EX_biomass_e_in"] = "EX_biomass(e)";
  FBAmet["sink_pheme_c_in"] = "sink_pheme(c)";
  FBAmet["EX_cys_L_e_in"] = "EX_cys_L(e)";
  FBAmet["EX_trp_L_e_in"] = "EX_trp_L(e)";
  FBAmet["EX_val_L_e_in"] = "EX_val_L(e)";
  FBAmet["EX_ile_L_e_in"] = "EX_ile_L(e)";
  FBAmet["EX_leu_L_e_in"] = "EX_leu_L(e)";
  FBAmet["EX_pro_L_e_in"] = "EX_pro_L(e)";
  
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

double FBA(double *Value,
           vector<class FBGLPK::LPprob>& vec_fluxb,
           map <string,int>& NumTrans,
           map <string,int>& NumPlaces,
           const vector<string> & NameTrans,
           const struct InfTr* Trans,
           const int T,
           const double& time) {
  
  if(Flag == -1) 
    
    init_data_structures(Trans, NumTrans);
  
  // Definition of the function exploited to calculate the rate,
  // in this case for semplicity we define it throught the Mass Action law
  
  cout << "FBAtime: " << FBAtime << endl;
  
  double nBac = Value[NumPlaces.find("CD") -> second];
  cout<< "nBac [cell]: " << nBac << endl;
  
  double DrugPlace = Value[NumPlaces.find("Drug") -> second];
  cout<< "DrugPlace [pmol]: " << DrugPlace << endl;
  
  if(FBAtime != time) {
    
    for (map<string, string>::iterator p = FBAmet.begin();
         p != FBAmet.end(); ++p) {
      
      // Saving the reactions indexes from the map
      int index = vec_fluxb[0].fromNametoid(p->second);
      string TypeBound = "GLP_DB";
      
      double Ub = vec_fluxb[0].getUpBounds(index);
      double Lb = vec_fluxb[0].getLwBounds(index);
      
      double Biom = Value[NumPlaces.find("BiomassCD") -> second];
      
      if(p->first == "EX_biomass_e_in") {
        
        cout << "Biom: " << Biom << endl;
        cout << "gDW_CDmax - Biom: " << (gDW_CDmax - Biom) << endl;
        
        if ((gDW_CDmax - Biom) > 1e-12) {
          Ub = (gDW_CDmax - Biom);
        } else {
          Ub = 1e-12;
        }
        
        cout << "Trans: " <<  p->first << ", Met [pg]: " << Biom << ";" << endl;
        cout << "Trans: " <<  p->first << ", Lb [mmol]: " << Lb << ";" << endl;
        cout << "Trans: " <<  p->first << ", Ub [mmol]: " << Ub << ";" << endl;
        
      } else {
        
        if (p->first == "sink_pheme_c_in") {
          
          double Met = Value[NumPlaces[FBAplace[p->first]]];
          
          if(nBac < 1) {
            Lb = -(Met*1e+09)*nBac;
          } else {
            Lb = -((Met*1e+09)/nBac*Biom*1e-12);
            }
          
          if(abs(Lb) > (Met*1e+09)){
            if(-(Met*1e+09) >= vec_fluxb[0].getLwBounds(index)){
              Lb = -(Met*1e+09);
            } else {
              Lb = vec_fluxb[0].getLwBounds(index);
            }
          }
          
          Ub = 1e-04;
          
          cout << "Trans: " <<  p->first << ", Met [mmol]: " << Met << ";" << endl;
          cout << "Trans: " <<  p->first << ", Met [pmol]: " << (Met*1e+09) << ";" << endl;
          cout << "Trans: " <<  p->first << ", Lb [mmol/gDW]: " << Lb << ";" << endl;
          cout << "Trans: " <<  p->first << ", Ub [mmol/gDW]: " << Ub << ";" << endl;
          
        } else {
          
          double Met = Value[NumPlaces[FBAplace[p->first]]];
          
          if(nBac < 1) {
            Lb = -((Met*c)/Na)*nBac;
          } else {
            Lb = -((Met*c)/Na)/(nBac*Biom*1e-12);
          }
          
          if(abs(Lb) > ((Met*c)/Na)){
            if(-((Met*c)/Na) >= vec_fluxb[0].getLwBounds(index)){
              Lb = -((Met*c)/Na);
            } else {
              Lb = vec_fluxb[0].getLwBounds(index);
            }
          }
          
          // essentail aa could be only produced and not consumed
          Ub = 1e-12;
          
          cout<< "Trans: " <<  p->first << ", Met [molecule]: " << Met*c <<";" << endl;
          cout<< "Trans: " <<  p->first << ", Met [mmol]: " << ((Met*c)/Na) <<";" << endl;
          cout << "Trans: " << p->first << ", Lb [mmol/gDW]: " << Lb << ";" << endl;
          cout << "Trans: " << p->first << ", Ub [mmol/gDW]: " << Ub << ";" << endl;
          
        }
      }
      
      vec_fluxb[0].update_bound(index, TypeBound, Lb, Ub);
      
    }
    
    vec_fluxb[0].solve();
    Vars = vec_fluxb[0].getVariables();
    
    FBAtime = time;
    
  }
  bool Out = 0;
  bool In = 0;
  
  string str = NameTrans[T];
  
  if(str.find("_out") != string::npos){
    // replace '_out' -> '_in
    str = std::regex_replace(str, std::regex("_out"), "_in");
    Out = 1;
  } else {
    In = 1;
  }
  
  int index = vec_fluxb[0].fromNametoid(FBAmet.find(str) -> second);
  
  cout<<"\nSolution:\n\n";
  
  double Biom = Value[NumPlaces.find("BiomassCD") -> second];
  
  if (str == "EX_biomass_e_in") {
    
    rate = Vars[index];
    
    // to convert moles to grams you need to multiply the MW by concentration:
    
    if((Out) && (rate > 0))
      rate = (MWbio*1e-03*(rate*Biom*1e-12))*1e+12;
    else if((Out) && (rate < 0))
      rate = 0;
    else if((In) && (rate > 0))
      rate = 0;
    else if((In) && (rate < 0))
      rate = -((MWbio*1e-03*(rate*Biom*1e-12))*1e+12);
    
    cout << "Firing transition: " << NameTrans[T] << endl;
    cout << "transition associated to: " << FBAmet.find(str) -> second << endl;
    cout << "flux estimated from fba [mmol/gDW*h]:" << Vars[index] << endl;
    cout << "rate [pg]: " << rate << endl;
    
  } else if(str.find("_L_e") != string::npos) {
    
    rate = Vars[index];
    
    // convert mmol to pmol
    // trans_in when is neg, otherwise trans_out
    
    if((Out) && (rate > 0))
      rate = rate*nBac*Biom*1e-12*(Na/c);
    else if((Out) && (rate < 0))
      rate = 0;
    else if((In) && (rate > 0))
      rate = 0;
    else if((In) && (rate < 0))
      rate = -(rate*nBac*Biom*1e-12*(Na/c));
    
    cout << "Firing transition: " << NameTrans[T] << endl;
    cout << "transition associated to: " << FBAmet.find(str) -> second << endl;
    cout << "flux estimated from fba [mmol/gDW*h]:" << Vars[index] << endl;
    cout << "rate [molecule]: " << (rate*c) << endl;
    
  } else {
    
    rate = Vars[index];
    
    if((Out) && (rate > 0))
      rate = rate*1e+09*nBac*Biom*1e-12;
    else if((Out) && (rate < 0))
      rate = 0;
    else if((In) && (rate > 0))
      rate = 0;
    else if((In) && (rate < 0))
      rate = -(rate*1e+09*nBac*Biom*1e-12);
    
    cout << "Firing transition: " << NameTrans[T] << endl;
    cout << "transition associated to: " << FBAmet.find(str) -> second << endl;
    cout << "flux estimated from fba [mmol/gDW*h]:" << Vars[index] << endl;
    cout << "rate [pmol]: " << rate << endl;
    
  }
  
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
    double g = 0;
    
    double DamagePlace = Value[NumPlaces.find("Damage") -> second];
    double PercDamage = DamagePlace/DAMAGEmax;
    
    if((PercDamage < 0) && (PercDamage < 0.1)){
      g = 0;
    } else if((PercDamage > 0.1) && (PercDamage <= 0.7)){
      g = 1/3;
    } else if(PercDamage > 0.7){
      g = 1;
    }
    
    rate = g*Inflammation;
    
    cout << "Inflammation [pmol]: " << rate << endl;
    
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
    
    if((Biom - gDW_CDmin) < 1e-12){
      if(nBac < 1) {
        half_life = 0.0;
      }
      
      cout<< "half_life: " <<  half_life << endl;
      rate = half_life*nBac;
      
    }
    
    cout<< "DeathBac [cells]: " << rate << endl;
    
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
    
    double rate = 0;
    
    double Biom = Value[NumPlaces.find("BiomassCD") -> second];
    double nBac = Value[NumPlaces.find("CD") -> second];
    
    cout<< "Biom [pg]: " << Biom << endl;
    cout<< "nBac [cell]: " << nBac << endl;
    
    if((Biom - gDW_CDmin) > 1e-12){
      
      rate = nBac*rCDdup*(1 - (nBac/nBacMax));
    } 
    
    cout << "Biom - gDW_CDmin: " << Biom - gDW_CDmin << endl;
    cout << "Dup rate [cell]: " << rate << endl;
    
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
    
    if(Flag == -1) init_data_structures(Trans, NumTrans);
    
    double rate = 0;
    
    double fbaEXBiomassOut = 0;
    double BiomassCDPlace = Value[NumPlaces.find("BiomassCD") -> second];
    
    if(FBAtime != time) {
      fbaEXBiomassOut = FBA(Value, 
                            vec_fluxb,
                            NumTrans, 
                            NumPlaces, NameTrans, Trans, 
                            NumTrans.find("EX_biomass_e_in") -> second, time);
    } else {
      int index = vec_fluxb[0].fromNametoid(string("EX_biomass(e)"));
      fbaEXBiomassOut = Vars[index];
    }
    
    cout<< "fbaEXBiomassOut [mmol/gDW*h]: " << fbaEXBiomassOut <<";" << endl;
    
    if(fbaEXBiomassOut <= 1e-04){
      // RCD [pg/h]
      rate = RCD*BiomassCDPlace;
      }
    
    cout<< "Starv rate [pg]: " << rate <<";" << endl;
    
    return(rate);
    
    }

// Starv Death4Treat

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
  double drug = Value[NumPlaces.find("Drug") -> second];
  
  rate = Death4Treat*nBac*drug;
  
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
  double Drug = Value[NumPlaces.find("Drug") -> second];
  
  double rate = Drug*nBac*efflux;
  
  return(rate);
}
  