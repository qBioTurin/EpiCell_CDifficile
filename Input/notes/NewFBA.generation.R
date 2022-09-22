
# The aim of this method is to generate a suitable ".Rdata" file model that 
# can be compiled into *FBAmodel* file subjected to epimod framework through model.generation()

# Here, model template lower bounds will be adjusted for the exchange reactions 
# of each metabolite with respect to the type of diet simulated.

# This methos also thought to implemented diet formulation obtained from VMH in models

FBA.generation = function(RDataPresent, NutritionModeling, ShowExBounds, write) {
  
  #'     Inputs: 
  #'             
  #'             RDataPresent: a logical value indicating whether FBAmodel should be generated starting from 
  #'             .mat file or the .RData is already available.
  #'             
  #'             NutritionModeling: character string specifying the diet implemetnaiton method used. 
  #'             NutritionModeling = "Macronutriets" (uptake rates were calculated by multiplying percentages 
  #'             of each macro-nutrient by a basal uptake rate B (data contained in "basal_uptake_rates.xlsx")).
  #'             Finally, uptake rates were varied by adjusting lower bounds for the exchange reactions of 
  #'             each metabolite with respect to the type of diet simulated.
  #'             NutritionModeling = "VMH" (uptake rates were calculated by replace template lower bounds for the exchange 
  #'             reactions of each metabolite with respect to the type of diet simulated)
  #'             
  #'             BasalUptakeRatesPresent: a logical value indicating whetherbasal uptake reta data are available,
  #'             if not they are randomly gerated based of the exchane reactions in model. 
  #'             
  #'             ShowExBounds: a logical value indicating whether model's bound were printed to check
  
  setwd("~/Documents/ClostridiumDiff_FBAandPN")
  
  if (RDataPresent == F) {
    
    # the generic function readMATmod imports matlab cobra models into sybil model files
    source("./code/readMATmod.R")
    
    model.name = readline(prompt = "Enter GSMM name (e.g. CD196HemeSink): ")
    path = readline(prompt = "Enter full directory path pointing the target 'model.mat': ")
    
    model.mat = readMATmod(path)
    R.matlab::Matlab$startServer()
    matlab <- R.matlab::Matlab()
    isOpen <- open(matlab)
    
    # Setting COBRA toolbox wd and starting COBRA
    R.matlab::evaluate(matlab, "cd /home/riccardo/cobratoolbox")
    R.matlab::evaluate(matlab, "initCobraToolbox")
    R.matlab::evaluate(matlab, paste("cd /home/riccardo/Documents/ClostridiumDiff_FBAandPN/EpiCell_CDifficile/Input/CDmodels/",
                             model.name, sep = ""))
    
    R.matlab::evaluate(matlab, paste("model = importdata('", model.name, ".mat');", sep = ""))
    # Building the rxnGeneMat reaction-to-gene mapping in sparse matrix form.
    R.matlab::evaluate(matlab,"model = buildRxnGeneMat(model);")
    R.matlab::evaluate(matlab, "rxnGeneMat = model.rxnGeneMat;")
    R.matlab::evaluate(matlab, "save('rxnGeneMat.mat', 'rxnGeneMat');")
    # Generating model's reaction names 
    R.matlab::evaluate(matlab, "model.rxnNames(cellfun('isempty', model.rxnNames)) = {'NA'};")
    R.matlab::evaluate(matlab, "writecell(model.rxnNames,'rxnNames.xlsx','Sheet',1);")
    
    # R.matlab::evaluate(matlab, paste("writeCbModel(model, 'format', 'mat', 'fileName','", 
    #                                  model.name, ".mat');", sep = ""))
    
    R.matlab::evaluate(matlab, paste("save('", model.name, ".mat', 'model');", sep = ""))
    # close the MATLAB server
    close(matlab)
    file.remove("InputStreamByteWrapper.class")
    file.remove("MatlabServer.m")
    # Saving on object in RData format
    model.mat = readMATmod(path)
    save(model.mat, 
         file = paste("./EpiCell_CDifficile/Input/CDmodels/", 
                      model.name, "/", model.name, ".RData", sep = ""))
    rm(model.mat) 
    
    # To load the data again
    load(file = paste("./EpiCell_CDifficile/Input/CDmodels/", model.name, "/", model.name, ".RData", sep = ""))
    rxnGeneMat <- R.matlab::readMat(paste("./EpiCell_CDifficile/Input/CDmodels/", model.name, "/", 
                                          "rxnGeneMat.mat", sep = ""))
    rxnGeneMat = rxnGeneMat[[1]]
    # changing model identification name
    model.mat@mod_desc = model.name
    # Searching for metabolites, which appear only once in the stoichiometric matrix.
    model.mat@met_single = sybil::singletonMetabolites(model.mat, retIds = FALSE)[["smet"]]
    model.mat@react_single = sybil::singletonMetabolites(model.mat, retIds = FALSE)[["sreact"]]
    model.mat@met_de = sybil::deadEndMetabolites(model.mat, retIds = FALSE)[["dem"]]
    model.mat@react_de = sybil::deadEndMetabolites(model.mat, retIds = FALSE)[["der"]]
    # Replacing empty slot
    model.mat@rxnGeneMat = rxnGeneMat
    
    # print reaction equation from reaction name
    # The function modelorg2tsv writes the content of an instance of class modelorg to text files in a
    # character-separated value format adopted from the BiGG database output.
    sybil::modelorg2tsv(model.mat, paste("BIGGdata", model.mat@mod_desc, sep = "_"), 
                        extMetFlag = "b", fielddelim = "\t", entrydelim = ", ", makeClosedNetwork = FALSE, 
                        onlyReactionList = FALSE, minimalSet = TRUE, 
                        fpath = paste("./EpiCell_CDifficile/Input/CDmodels/", model.name, sep = ""))
    
    # loading model's reaction metadata and default constraint values
    assign("BiGGdata", read.table(file = paste("./EpiCell_CDifficile/Input/CDmodels/",
                                               model.mat@mod_desc, "/BIGGdata_",
                                               model.mat@mod_desc, "_react.tsv", 
                                               sep = ""), sep = '\t', header = TRUE))
    
    # reading full reaction names and complete metadata
    rxnNames <- readxl::read_excel(paste("./EpiCell_CDifficile/Input/CDmodels/", model.mat@mod_desc,
                                         "/rxnNames.xlsx", sep = ""), col_names = FALSE)
    
    BiGGdata = cbind(rxnNames, BiGGdata)
    colnames(BiGGdata)[1] <- "rxnNames"
    
    # converting and saving BIGG data to a data frame
    write.table(BiGGdata, file = paste("./EpiCell_CDifficile/Input/CDmodels/", 
                                       model.mat@mod_desc, "/BIGGdata_", 
                                       model.mat@mod_desc, "_react.tsv", 
                                       sep = ""), 
                quote=FALSE, sep='\t', col.names = NA)
    
    # Saving the modified model
    save(model.mat, file = paste("./EpiCell_CDifficile/Input/CDmodels/", 
                                 model.name, "/", model.name, ".RData", sep = ""))
    
  } else {
    
    path_rdata = readline(prompt = "Enter full directory path pointing the target 'model.RData': ")
    cat("RData loading ...")
    load(path_rdata)
    
  }
  
  switch(NutritionModeling,
         
         Macronutriets = {
           cat("Macronutriets method selected ... ")
           
           cat("Example of percentages of macro-nutrients:")
           paste(c("Western C50,F35,P15", "High Carb C75,F5,P20", "High Fat C5,F75,P20", "High Prot C5,F20,P75"))
           cat("                                            ", "\n")
           
           Carb = readline(prompt = "Enter percentages of macro-nutrients (C): ")
           Carb = as.integer(Carb)
           Fat = readline(prompt = "Enter percentages of macro-nutrients (F): ")
           Fat = as.integer(Fat)
           Prot = readline(prompt = "Enter percentages of macro-nutrients (P): ")
           Prot = as.integer(Prot)
           
           sumok = Carb + Fat + Prot
           
           if(sumok != 100) {
             stop("percentages NO add up to 100")
           }
           
           # path_basal: full directory path pointing the target basal_uptake_rates.xlsx
           path_basal_rate = readline(prompt = "Enter basal uptake rate file full path: ")
           # EpiCell_CDifficile/Input/Diets/macronutrients/basal_rate_uptakes.csv
           
           # e.g https://doi.org/10.1371/journal.pcbi.1006714
           
           basal_rate_uptakes <- read_csv(path_basal_rate)
           
           ex <- sybil::findExchReact(model.mat)
           ex = ex@react_id
           
           value = basal_rate_uptakes$uptakes[which(basal_rate_uptakes$rnx_name %in% ex)]
           name = basal_rate_uptakes$rnx_name[which(basal_rate_uptakes$rnx_name %in% ex)]
           type = basal_rate_uptakes$NutrientType[which(basal_rate_uptakes$rnx_name %in% ex)]
           
           not_name = ex[which(!(ex %in% basal_rate_uptakes$rnx_name))][which(stringr::str_detect(a, "EX"))]
           value_not_name = rep(0.01, length(not_name))
           type_not_name = rep("other", length(not_name))
           
           df = data.frame(name = c(name, not_name),
                           value = c(value, value_not_name),
                           type = c(type, type_not_name))
           
           carbuptake = ((Carb/100)*df$value[which(df$type == "carb")])/-0.1 # multiply basal rates for relevant macronutrients by their respective random percentages
           fatuptake = ((Fat/100)*df$value[which(df$type == "fat")])/-0.1 # fats and protein rates are multiplied by 2 again to scale up to 100%
           protuptake = ((Prot/100)*df$value[which(df$type == "prot")])/-0.1 # divide by -0.1 to rescale uptake values(since metabolites are taken up, values should be negative)
           
           otheruptake = ((1/100)*df$value[which(df$type == "other")])/-0.1
           
           # recombine into single vector of uptakes
           df$value = c(carbuptake, fatuptake, protuptake, otheruptake)
           
           index = c()
           
           uppbndBiomass = readline(prompt = "Enter a Upper Bound value (ub >= 0) for EX_biomass(e):")
           # e.g. 1.082392115563062
           lowbndBiomass = readline(prompt = "Enter a Lower Bound value (lb <= 0) for EX_biomass(e):")
           # e.g lbNOTdiet
           
           for (i in 1:length(df$name)) {
             index[i] = which(model.mat@react_id == df$name[i])
           }
           
           model.mat@lowbnd[index] = df$value
           
           model.mat@uppbnd[which(model.mat@react_id == "EX_biomass(e)")] = uppbndBiomass
           model.mat@lowbnd[which(model.mat@react_id == "EX_biomass(e)")] = lowbndBiomass
           
           ans = readline(prompt = "Want to save .RData copy (Y or N) ?")
           
           if (ans == "Y") {
             file_name = readline(prompt = "enter file name (e.g. mod_descWithEnv.RData): ")
             save(model.mat, file = file_name)
           }
           
           },
         
         VMH = {
           
           cat("VMH method selected ...")
           
           path_diet = readline(prompt = "Enter full directory path pointing the diet data: ")
           # e.g. ~/Documents/ClostridiumDiff_FBAandPN/EpiCell_CDifficile/Input/Diets/diets.RData
           load(path_diet)
           cat("Diet formulation provided:", names(diets), "\n")
           cat("                                            ", "\n")
           diet_name = readline(prompt = "Enter diet name: ")
           
           d = diets[[paste(diet_name)]]
           
           cat("                                            ", "\n")
           cat("Diet is envisaged. Diet formulation: ", diet_name, "\n")
           
           ex <- sybil::findExchReact(model.mat)
           ex = ex@react_id
           
           value = d$`Flux Value`[which(d$Reaction %in% ex)]
           name = d$Reaction[which(d$Reaction %in% ex)]
           index = c()
           
           for (i in 1:length(name)) {
             index[i] = which(model.mat@react_id == name[i])
           }
           
           react_id_ex = which(model.mat@react_id %in% unlist(ex))
           
           cat("It could be that not all the Exchanges included in the model are covered by dietary data ... \n")
           
           lbNOTdiet = readline(prompt = "Enter Lower Bound value (lb <= 0) for those metabolites not envisaged in the diet:")
           # e.g. lbNOTdiet = -1e-05
           lbNOTdiet = as.double(lbNOTdiet)
           
           uppbndBiomass = readline(prompt = "Enter Upper Bound value (ub >= 0) for EX_biomass(e):")
           # e.g. 1.082392115563062
           lowbndBiomass = readline(prompt = "Enter Lower Bound value (lb <= 0) for EX_biomass(e):")
           # e.g lbNOTdiet
           
           uppbndBiomass = as.double(uppbndBiomass)
           lowbndBiomass = as.double(lowbndBiomass)
           
           model.mat@lowbnd[react_id_ex[which(model.mat@lowbnd[react_id_ex] != 0)]] =
             rep(lbNOTdiet, length(model.mat@lowbnd[react_id_ex[which(model.mat@lowbnd[react_id_ex] != 0)]]))
           
           # total initial gDW = 4.57e07*0.5*1e-12 = 2.285e-05 [gDW]
           nBac = 4.57e07; Biom = 0.5; 
           # 1e-04 is a rescaling factor 
           model.mat@lowbnd[index] = (value*1e-04)/(nBac*Biom*1e-12) # [mmol/per community]
           
           Ub = readline(prompt = "Enter Upper bounds values (ub >= 0):")
           # e.g. Ub = 0.1
           Ub = as.double(Ub)
           
           model.mat@uppbnd[react_id_ex[which(model.mat@uppbnd[react_id_ex] != 0)]] =
             rep(Ub, length(model.mat@uppbnd[react_id_ex[which(model.mat@uppbnd[react_id_ex] != 0)]]))
           
           aaUb = readline(prompt = "Enter Upper bounds values for AAs (ub >= 0):")
           aaUb = as.double(aaUb)
           
           model.mat@uppbnd[which(model.mat@react_id %in% 
                                    c("EX_pro_L(e)", "EX_leu_L(e)",
                                      "EX_ile_L(e)", "EX_val_L(e)", 
                                      "EX_trp_L(e)", "EX_cys_L(e)"))] = rep(aaUb, 6)
           
           model.mat@uppbnd[which(model.mat@react_id == "EX_biomass(e)")] = uppbndBiomass
           model.mat@lowbnd[which(model.mat@react_id == "EX_biomass(e)")] = lowbndBiomass
           
           ans = readline(prompt = "Want to save .RData copy (Y or N) ?")
           
           model.mat <<- model.mat
           
           if (ans == "Y") {
             file_name = readline(prompt = "enter file name (e.g. mod_descWithEnv.RData): ")
             save(model.mat, file = file_name)
             }
           
           },
         
         Template = {
           cat("Diet is not envisaged. Bounds of EX reactions are kept as model template ...")
           cat("                                                                            ")
         },
         
         {
           stop("incorrect diet method")
           }
         
         )
  
  if(ShowExBounds == T) {
    
    source("./EpiCell_CDifficile/Rfunction/flux_balance_fba.R")
    
    sol = flux_balance(model = model.mat, ShowExBounds = F,
                       typename = model.mat@mod_desc, diet = F, 
                       lbNOTdiet = lbNOTdiet, symmetric = F, 
                       uppbndBiomass = diet_name)
    
    v = sol[[1]][["fluxes"]]
    
    EX_biomass = c()
    EX_ile = c(); EX_pro = c(); EX_leu = c()
    EX_val = c(); EX_cys = c(); EX_trp = c()
    Sink_heme = c()
    
    EX_biomass = v[which(model.mat@react_id == 'EX_biomass(e)')]
    EX_leu = v[which(model.mat@react_id == 'EX_leu_L(e)')]
    EX_pro = v[which(model.mat@react_id == 'EX_pro_L(e)')]
    EX_ile = v[which(model.mat@react_id == 'EX_ile_L(e)')]
    EX_val = v[which(model.mat@react_id == 'EX_val_L(e)')]
    EX_cys = v[which(model.mat@react_id == 'EX_cys_L(e)')]
    EX_trp = v[which(model.mat@react_id == 'EX_trp_L(e)')]
    Sink_heme = v[which(model.mat@react_id == 'sink_pheme(c)')]
    
    DietData = data.frame(EX_biomass,
                          EX_pro, EX_leu, EX_val,
                          EX_ile, EX_cys, EX_trp,
                          Sink_heme, DietSet = diet_name)
    
    DietData = tidyr::gather(DietData, key, Value, -DietSet)
    View(DietData)
    
    ggplot(DietData, aes(x = DietSet, y = Value, fill = key)) + 
      geom_bar(aes(x = DietSet, y = Value, fill = key), 
               stat = "identity", position = "dodge") +
      geom_text(aes(label = round(Value, 9)), size = 1.8, 
                position=position_dodge(width = 0.9), vjust = -2.5) +
      scale_y_continuous("Estimated Flux [mmol/gDW*h", expand = c(0, 0)) + 
      scale_x_discrete("Reaction") +
      scale_fill_manual(values = c("#999999", "#E69F00", "#56B4E9", "#009E73",
                                   "#F0E442", "#0072B2", "#D55E00", "#CC79A7")) +
      hrbrthemes::theme_ipsum() +
      theme(axis.text.x = element_text(angle = 0, hjust = 1, vjust = 0),
            axis.line = element_blank(), axis.ticks.x = element_blank())
    
    }
  
  # FBAfile generation phase
  
  # Problem: S*x = b -. size x = n and size of S = m X n
  # With this function the file to pass to the GLPKsolver is generated:
  
  # 1 row) n_row n_col GLP_MAX (or GLP_MIN if the obj has to max or min)
  # 2 row) the coeff for defining the obj function
  
  #           (the num of coeff has to be == length of x)
  
  # m rows)   definition of the S row bounds (i.e. b)
  # n rows)   definition of the x bounds
  # m*n rows) definition of the S coeffs: row_index col_index value
  
  #### Start writing the csv file to pass to GLPK solver
  
  Matlab.file = model.mat
  
  ReactionsNames <- unlist(Matlab.file@react_id)
  # [e] = extracellular metabolites
  # [c] = cytosolic metabolites
  ReagentsNames <- unlist(Matlab.file@met_id)
  
  as.matrix(Matlab.file@S) -> S
  
  ncol=length(S[1,])
  nrow=length(S[,1])
  
  matrix(0,nrow = ncol, ncol = 1) -> b
  as.matrix(Matlab.file@lowbnd) -> lb
  as.matrix(Matlab.file@uppbnd) -> ub
  c(Matlab.file@obj_coef) -> obj
  
  rb <- cbind(b,b)
  cb <- cbind(lb,ub)
  
  if(write) {
    fileName = "FBAmodel"
    write(paste(ReactionsNames, collapse = " "), file = fileName)
    write(paste0(nrow," ; ",ncol," ; ","GLP_MAX" ),file = fileName, append = T)
    write(paste(obj, collapse = " "),file = fileName,append = T )
    
    for(i in 1:nrow) {
      write(paste0("GLP_F"," ; ",paste0(rb[i,],collapse = " ; ") ),file = fileName ,append = T)
    }
    
    for(j in 1:ncol){
      write(paste0("GLP_DB"," ; ",paste0(cb[j,],collapse = " ; ")),file = fileName ,append = T)
    }
    
    for(i in 1:nrow) {
      for(j in 1:ncol){
        write(paste0(i," ; ", j, " ; ", S[i,j]),file = fileName ,append = T)
      }
    }
    
    FinalMatrix = readr::read_csv("FBAmodel", col_names = FALSE)
    
    # save(FinalMatrix,file = "Input/FinalMatrix.RData")
    
    return(0)
    
  } else {
    
    FinalMatrix <- paste0(nrow," ; ",ncol," ; ","GLP_MAX" )
    FinalMatrix <- c(FinalMatrix, paste(obj,collapse = " "))
    for(i in 1:nrow){
      FinalMatrix <- c(FinalMatrix, paste0("GLP_F"," ; ",paste0(rb[i,],collapse = " ; ") ))
    }
    for(j in 1:ncol){
      FinalMatrix <- c(FinalMatrix,paste0("GLP_DB"," ; ",paste0(cb[j,],collapse = " ; ")) )
    }
    for(i in 1:nrow){
      for(j in 1:ncol){
        FinalMatrix <- c(FinalMatrix,paste0(i," ; ", j, " ; ", S[i,j]))
      }
    }
  }
  
  return(matrix(FinalMatrix,ncol=1))
  
}

