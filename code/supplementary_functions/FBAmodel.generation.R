
devtools::install_github("https://github.com/qBioTurin/epimod_FBAfunctions")

# The aim of this method is to generate a suitable ".Rdata" file model that 
# can be compiled into *FBAmodel* file subjected to epimod framework through model.generation()

# Here, model template lower bounds will be adjusted for the exchange reactions 
# of each metabolite with respect to the type of diet simulated.

# This method also thought to implemented diet formulation obtained from VMH in models

FBAmodel.gen = function(model.type,
                        model.name,
                        NutritionModeling,
                        write,
                        saveRData,
                        diet_name,
                        change_manual,
                        lb_manual,
                        gene.exp,
                        gene.exp_sample,
                        wd) {
  
  source(paste0(wd, "/code/processing_functions/ExtractEx.R"))
  source(paste0(wd, "/code/supplementary_functions/Functions.R"))
  
  source(paste0(wd, "/code/processing_functions/FBAmodel_Metadata.R"))
  source(paste0(wd, "/code/processing_functions/createReactionString.R"))
  source(paste0(wd, "/code/processing_functions/modelorg2tsv.R"))
  source(paste0(wd, "/code/processing_functions/Eval_ReactionExpression.R"))
  source(paste0(wd, "/code/processing_functions/findReactEq.R"))

  mat.file = paste0(wd, "/input/", model.type, "/", model.name, "/", model.name, ".mat")
  model.file = paste0(wd, "/input/", model.type, "/", model.name, "/", model.name, ".RData")
  all.react.file = paste0(wd, "/input/", model.type, "/", model.name, "/all.react.rds")
  diet_file = paste0(wd, "/input/diets/diets.RData")
  
  # single cell C.diff minimum observed biomass
  gDW_CD_min = ComputeBiomassBac(diameter = 0.3, len = 2) # (pg)
  # # single cell C.diff mean observed biomass
  gDW_CDmean = ComputeBiomassBac(diameter = 0.5, len = 5.5) # (pg)
  # # single cell C.diff maximun observed biomass
  gDW_CDmax = ComputeBiomassBac(diameter = 0.7, len = 9) # (pg)
  
  Vcd = ((pi*(0.5)^2)/4)*(((4*0.5)/6) + 5.5) # (μm^3)
  Biom = (1/3)*Vcd*1.3 # (pg)
  nBac = 45700000 # (bac)
  Delta = 1e-04 # is a rescaling factor (unit)
  
  # generate S3 R object FBA_model
  model = FBA4Greatmod.generation(fba_mat = mat.file)
  
  # moving files
  system(paste0("cd ", wd))
  
  files = c("geneAssociation.rds", 
            "geni.rds", "officialName.rds", 
            "subsystem.rds", "genesFromGeneAss.rds",
            "met_KEGGID.rds", "rxn_KEGGID.rds")
  
  for(f in files) {
    system(paste0("mv ", f, " ", wd, "/input/", model.type, "/", model.name, "/", f))
  } 
  
  # saving RData FBA model
  save(model, file = model.file)
  
  # generate FBA model metadata and saving FBA model reaction data
  bigg.path = paste0(wd, "/input/", model.type, "/", model.name, 
                     "/BIGGdata_", model.name, "_react.tsv")
  
  metadataEQ(model = model, 
             model.name = model.name,
             prefix = paste("BIGGdata", model.name, sep = "_"),
             suffix = "tsv",
             extMetFlag = "b", 
             fielddelim = "\t", 
             entrydelim = ", ", 
             makeClosedNetwork = FALSE, 
             onlyReactionList = FALSE, 
             minimalSet = TRUE, 
             fpath = paste0(wd, "/input/", model.type, "/", model.name))
  
  all.react = FBAmodel.metadata(model = model, 
                                model.name = model.name, 
                                wd = wd, 
                                bigg.path = bigg.path)
  
  saveRDS(all.react, all.react.file)
  
  # saving biochemical reaction equations
  BIGGdata = read.delim2(bigg.path)
  
  saveRDS(BIGGdata$equation, file = paste0(wd, "/input/", model.type, "/", model.name, "/equation.rds"))
  
  # model annotation
  all.react = readRDS(all.react.file)
  br = all.react[[2]]
  all.react = all.react[[1]]
  
  switch(NutritionModeling,
         
         VMH = {
           
           cat("VMH method selected ...")
           
           load(diet_file)
           d = diets[[diet_name]]
           
           cat("                                                             ", "\n")
           cat("VMH diet envisaged: diet formulation ...", "'", diet_name, "'", "\n")
           
           value = d$`Flux Value`[which(d$Reaction %in% br$react.id)]
           react_id_ex = d$Reaction[which(d$Reaction %in% br$react.id)]
           
           ReactionsNames <<- unlist(model@react_id)
           
           index = c()
           
           for (i in 1:length(react_id_ex)) {
             index[i] = which(model@react_id == react_id_ex[i])
           }
           
           lbNOTdiet = -(Delta)/(nBac*Biom*1e-12) # (mmol/per community)
           
           t_lb = dplyr::filter(all.react, React_Lb != 0 & ReactionDistrict == "boundary")
           
           model@lowbnd[t_lb$ReactionPos] = rep(lbNOTdiet, length(model@lowbnd[t_lb$ReactionPos]))
           model@lowbnd[index] = -(value*Delta)/(nBac*Biom*1e-12) # (mmol/per community)
           
           },
         
         Template = {
           cat("Diet is not envisaged. Bounds of EX reactions are kept as model template ...")
           cat("                                                                            ")
         },
         
         {
           stop("incorrect diet method")
         }
         
         )
  
  if (gene.exp == T) {
    
    # from An intact S-layer is advantageous to Clostridioides difficile within the host.
    # GSE205747
    # Bulk RNA-seq and DNA-Seq of Clostridioides difficile strains FM2-5, R20291, Rv9, Rv189 and RW2-5.
    # importing gene expression
    
    # genesFromGeneAss = readRDS(paste0(wd, "/input/", model.type, "/", model.name, "/genesFromGeneAss.rds"))
    # 
    # GE = as.data.frame(matrix(round(runif(length(genesFromGeneAss)*12, min = 1, max = 100), 3), ncol = 12))
    # GE = cbind(gene_id_model = genesFromGeneAss,
    #            GE)
    # 
    # write.table(GE,
    #             paste0(wd, "/input/", model.type, "/", model.name, "/gene_expression/Table_expression_data.txt"),
    #             row.names = F, col.names = T, quote = F)
    
    GE = read_table(paste0(wd, "/input/", model.type, "/", model.name, 
                           "/gene_expression/Table_expression_data.txt"),
                    col_types = cols(gene_id_model = col_character()))
    
    # gene expression file for the network parametrization
    mod_gprRules.file = paste0(wd, "/input/", model.type, "/", model.name, "/geneAssociation.rds")
    
    GSEfile = paste0(wd, "/input/", model.type, "/", model.name, "/gene_expression/GSEfile.RData")
    
    if (!file.exists(GSEfile)) {
      
      cl = makeCluster(getOption("cl.cores", detectCores()))
      clusterExport(cl, ls())
      
      GSE.allScenarios = parLapply(
        
        cl = cl, 
        X = colnames(GE)[2:(ncol(GE))],
        
        fun = function(patient){
          
          map_geID = GE[, c("gene_id_model", patient)]
          map_geID[which(is.na(as.double(map_geID[, 2]))), 2] = 1
          geni = map_geID[, 1]
          gene_expr = as.double(map_geID[, 2])
          
          ## gene expression influence
          colnames(map_geID) = c("gene", "val")
          
          # gene-reaction association rule for each reaction
          mod_gprRules = readRDS(mod_gprRules.file)
          
          # loop over the array of the gene expressions, in order to replace the names
          # of genes in geni_reazioni.mat with their values
          
          GSE = rep(0, length(mod_gprRules))
          
          for (i in 1:length(mod_gprRules)) {
            expression = mod_gprRules[i]
            if(expression != ""){
              # evaluation of reaction expression
              GSE[i] =  Eval.ReactionExpression(txt = expression,
                                                GeneExpr = as.data.frame(map_geID)
              )
            } else {
              # replaces all the empty cells of gene expression (e.g. exchange reactions) 
              # with 1, i.e. gene expressed normally
              GSE[i] = 1
            }
          }
          return(list(GSE = GSE))
        }
        
      )
      
      stopCluster(cl)
      
      names(GSE.allScenarios) = colnames(GE)[2:(ncol(GE))]
      
      save(GSE.allScenarios, file = GSEfile)
      
    }
    
    load(GSEfile)
    
    GSE = GSE.allScenarios[[gene.exp_sample]]$GSE
    
    for (i in 1:length(GSE)) {
      
      if (GSE[i] >= 1) {
        
        model@lowbnd[i] = model@lowbnd[i] * (1 + log(GSE[i]))
        model@uppbnd[i] = model@uppbnd[i] * (1 + log(GSE[i]))
        
      } else {
        
        model@lowbnd[i] = model@lowbnd[i] / (1 + abs(log(GSE[i])))
        model@uppbnd[i] = model@uppbnd[i] / (1 + abs(log(GSE[i])))
        
      }
    }
  }
  
  model@lowbnd[which(model@react_id %in% change_manual)] = lb_manual
  model@uppbnd[which(model@react_id == "EX_biomass_e")] = (gDW_CDmax - gDW_CDmean)
  
  model <<- model
  
  save(model, file = model.file)
  
  all.react = FBAmodel.metadata(model = model, 
                                model.name = model.name, 
                                wd = wd, 
                                bigg.path = bigg.path)
  
  saveRDS(all.react, paste0(wd, "/input/", model.type, "/", model.name, "/all.react_diet.rds"))
  
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
  
  ReagentsNames = unlist(model@met_id)
  # [e] = extracellular metabolites
  # [c] = cytosolic metabolites
  
  S = as.matrix(model@S)
  
  ncol = length(S[1, ])
  nrow = length(S[, 1])
  
  matrix(0, nrow = ncol, ncol = 1) -> b
  as.matrix(model@lowbnd) -> lb
  as.matrix(model@uppbnd) -> ub
  c(model@obj_coef) -> obj
  
  rb = cbind(b, b)
  cb = cbind(lb, ub)
  
  if(write) {
  
    fileName = paste0(model.name, ".txt")
    
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
    
    system(paste0("mv ", fileName, " ", wd, "/input/CompiledModels"))
    
    return(0)
    
  }
}

