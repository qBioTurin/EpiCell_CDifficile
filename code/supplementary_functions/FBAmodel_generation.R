
# The aim of this method is to generate a suitable ".Rdata" file model that 
# can be compiled into *FBAmodel* file subjected to epimod framework through model.generation()

# Here, model template lower bounds will be adjusted for the exchange reactions 
# of each metabolite with respect to the type of diet simulated.

# This method also thought to implemented diet formulation obtained from VMH in models

FBAmodel.gen = function(model_type,
                        model_name,
                        NutritionModeling,
                        write,
                        saveRData,
                        diet_name,
                        change_manual,
                        lb_manual,
                        gene.exp,
                        gene.exp_sample,
                        lb_off,
                        ub_manual,
                        GeneExpress_file,
                        wd) {
  
  source(paste0(wd, "/code/processing_functions/ExtractEx.R"))
  source(paste0(wd, "/code/processing_functions/FBAmodel_Metadata.R"))
  source(paste0(wd, "/code/processing_functions/createReactionString.R"))
  source(paste0(wd, "/code/processing_functions/modelorg2tsv.R"))
  source(paste0(wd, "/code/processing_functions/Eval_ReactionExpression.R"))
  source(paste0(wd, "/code/processing_functions/findReactEq.R"))
  
  source(paste0(wd, "/code/epimod_FBAfunctions/R/class_generation.R"))
  source(paste0(wd, "/code/epimod_FBAfunctions/R/FBAgreatmodeClass.R"))
  source(paste0(wd, "/code/epimod_FBAfunctions/R/readMat.R"))
  
  # setting files paths
  model_wd = paste0(wd, "/input/", model.type, "/", model.name, "/")
  mat.file = paste0(model_wd, "/", model_name, ".mat")
  model.file = paste0(model_wd, "/", model_name, ".RData")
  all.react.file = paste0(model_wd, "/all_react.rds")
  all.react.file_param = paste0(model_wd, "/param_all_react.rds")
  diet_file = paste0(wd, "/input/diets/diets.RData")
  bigg.path = paste0(model_wd, "/BIGGdata_", model_name, "_react.tsv")
  mod_gprRules.file = paste0(model_wd, "/GRrules.rds")
  mod_genes.file = paste0(model_wd, "/genesfromGRrules.rds")
  
  # generate S3 R object FBA_model
  model = FBA4Greatmod.generation(fba_mat = mat.file)
  
  # moving files
  system(paste0("cd ", wd))
  
  # save and move metabolic model fundamental fields
  files = c("GRrules.rds", 
            "geni.rds", "officialName.rds", 
            "subsystem.rds", "genesfromGRrules.rds",
            "met_KEGGID.rds", "rxn_KEGGID.rds")
  
  for(f in files) {
    system(paste0("mv ", f, " ", wd, "/input/", model_type, "/", model_name, "/", f))
  }
  
  metadataEQ(model = model, 
             model.name = model_name,
             prefix = paste("BIGGdata", model_name, sep = "_"),
             suffix = "tsv",
             extMetFlag = "b", 
             fielddelim = "\t", 
             entrydelim = ", ", 
             makeClosedNetwork = FALSE, 
             onlyReactionList = FALSE, 
             minimalSet = TRUE, 
             fpath = model_wd)
  
  all_react = FBAmodel.metadata(model = model, 
                                model.name = model_name, 
                                wd = wd, 
                                bigg.path = bigg.path)
  
  saveRDS(all_react, all.react.file)
  
  # saving biochemical reaction equations
  BIGGdata = read.delim2(bigg.path)
  saveRDS(BIGGdata$equation, file = paste0(model_wd, "/equation.rds"))
  
  br = all_react[[2]]
  all_react = all_react[[1]]
  
  t_lb = dplyr::filter(all_react, ReactionType == "Exchange" & all_react$React_Lb != 0)
  
  switch(NutritionModeling,
         
         VMH = {
           
           # Load diet file
           load(diet_file)
           d = diets[[diet_name]]
           
           # Print VMH diet envisaged
           cat("VMH diet envisaged: diet formulation ...'", diet_name, "'\n")
           
           # Get flux value and reaction ID
           value = d$`Flux Value`[d$Reaction %in% br$react.id] # (mmol/h*person)
           react_id_ex = d$Reaction[d$Reaction %in% br$react.id]
           
           # Get index of reaction ID
           index = sapply(react_id_ex, function(x) which(model@react_id == x))
           
           # Set lower bounds - new default
           if( !is.na(lb_off) ) {
             model@lowbnd[t_lb$ReactionPos] = rep(lb_off, dim(t_lb)[1]) 
           }
           
           # nBac as the number of bacterial cells in the volume
           nBac = 45700000 # (cell)
           # initial biomass as the mean biomass of a bacterial cell
           gDW_CDmean = 0.05727*1e-12 # (g)
           
           value = value / (nBac*gDW_CDmean) # (mmol/gDW*h)
           
           # Set lower bounds - diet
           model@lowbnd[as.vector(index)] = -value
           
           cat("model ex reacts lb ...'", dim(t_lb)[1], "'\n")
           cat("diet ex reacts lb ...'", length(d$Reaction), "'\n")
           cat("ex reacts lb constrained by diet...'", length(react_id_ex), "'\n")
           
         },
         
         Template = {
           
           if( is.na(lb_off) ) {
             # Diet is not envisaged. Bounds of EX reactions are kept as model template
             cat("Diet is not envisaged. Bounds of ex reacts are kept as model template ...\n")
           } else {
             # Diet is not envisaged. Bounds of EX reactions changing
             cat("Diet is not envisaged. Bounds of ex reacts lb set as ...", lb_off, "\n")
             
             pos = all_react[match(t_lb$React_ID, all_react$React_ID), ]$ReactionPos
             # Find these reactions in the model and set their lower bound to lb_off
             model@lowbnd[pos] = rep(lb_off, length(pos))
             
           }
         },
         
         {
           # Incorrect diet method exception
           stop("incorrect diet method")
         }
  )
  
  if ( gene.exp ) {
    
    # Genome-scale metabolic models (GSMMs) computationally describe 
    # gene-protein-reaction associations for entire metabolic genes in an organism.
    # These associations are formulated on the basis of genome annotation.
    
    # The gene expression data could be integrated with the GSMM. 
    # Each reaction in the model is associated with one or more genes
    # and the expression level of these genes can be used to set the bounds through each reaction.
    
    cat("Gene Expression data envisaged: ... \n")
    cat("total model reactions: ... ", dim(all_react)[1], "\n")
    
    GE = read.csv(GeneExpress_file)
    colnames(GE)[1] = "gene_id_model"
    
    cat("num model genes ...'", length(readRDS(mod_genes.file)), "'\n")
    cat("num model genes matched in expression dataset ...'", dim(GE)[1], "'\n")
    
    GSEfile = paste0(wd, "/input/gene_expression/GSE_", model_name, ".csv")
    
    if ( !file.exists(GSEfile) ) {
      
      library(parallel)
      
      cl = makeCluster(getOption("cl.cores", detectCores()))
      clusterExport(cl, ls())
      
      GSE.allScenarios = parLapply(
        
        cl = cl, 
        X = colnames(GE)[2:(ncol(GE))],
        
        fun = function(p){
          
          map_geID = GE[, c("gene_id_model", p)]
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
            tryCatch({
              if (expression != "") {
                # evaluation of reaction expression
                GSE[i] = Eval.ReactionExpression(
                  txt = expression, GeneExpr = as.data.frame(map_geID))
              } else {
                # replaces all the empty cells of gene expression (e.g. exchange reactions) 
                # with 1, i.e. gene expressed normally
                GSE[i] = 1
              }
            }, error = function(e) {
              cat("Error in expression ", i, ": ", e$message, "\n")
            })
          }
          return(list(GSE = GSE))
        }
      )
      
      stopCluster(cl)
      
      names(GSE.allScenarios) = colnames(GE)[2:(ncol(GE))]
      save(GSE.allScenarios, file = GSEfile)
      
    } else { load(GSEfile) }
    
    GSE = GSE.allScenarios[[gene.exp_sample]]$GSE
    
    # Adjusting the lower and upper bounds of a metabolic model based on gene expression data:
    #
    # (1) If the gene expression level GSE[i] is greater than or equal to 1, it means the gene is expressed. 
    # In this case, it increases the bounds of the ith reaction flux by (1 + log(GSE[i])). 
    # The use of the logarithm dampens the effect of very high expression levels.
    #
    # (2) If the gene expression level GSE[i] is less than 1, it means the gene is lowly expressed or not expressed. 
    # In this case, the code decreases the bounds of the ith reaction flux by dividing them by (1 + abs(log(GSE[i]))). 
    # The absolute value of the logarithm is used to avoid negative values, 
    # and adding 1 ensures that the bounds do not become infinite when the expression level is 0.
    
    for ( i in 1:length(GSE) ) {
      if ( GSE[i] >= 1 ) {
        
        model@lowbnd[i] = model@lowbnd[i] * (1 + log(GSE[i]))
        model@uppbnd[i] = model@uppbnd[i] * (1 + log(GSE[i]))
        
      } else {
        
        model@lowbnd[i] = model@lowbnd[i] / (1 + abs(log(GSE[i])))
        model@uppbnd[i] = model@uppbnd[i] / (1 + abs(log(GSE[i])))
        
      }
    }
  }
  
  pos_manual = all_react[match(change_manual, all_react$React_ID), ]$ReactionPos
  
  model@lowbnd[pos_manual] = lb_manual
  model@uppbnd[pos_manual] = ub_manual
  
  cat("ex reacts 'lb' and 'ub' manually constrained ...'", length(pos_manual), "'\n")
  
  all_react = FBAmodel.metadata(model = model, 
                                model.name = model_name, 
                                wd = wd, 
                                bigg.path = bigg.path)
  
  saveRDS(all_react, all.react.file_param)
  save(model, file = model.file)
  
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
  
  S = as.matrix(model@S)
  
  ncol = length(S[1, ])
  nrow = length(S[, 1])
  
  matrix(0, nrow = ncol, ncol = 1) -> b
  as.matrix(model@lowbnd) -> lb
  as.matrix(model@uppbnd) -> ub
  c(model@obj_coef) -> obj
  
  rb = cbind(b, b)
  cb = cbind(lb, ub)
  
  if ( write ) {
    # Set file name
    fileName = paste0(model_name, ".txt")
    
    # Write reactions names, number of rows and columns, and GLP_MAX
    write(paste(unlist(model@react_id), collapse = " "), file = fileName)
    write(paste0(nrow," ; ",ncol," ; ","GLP_MAX" ),file = fileName, append = T)
    write(paste(obj, collapse = " "),file = fileName,append = T )
    
    # Write GLP_F
    for (i in 1:nrow) {
      write(paste0("GLP_F"," ; ",paste0(rb[i,],collapse = " ; ") ),file = fileName ,append = T)
    }
    
    # Write GLP_DB
    for (j in 1:ncol){
      write(paste0("GLP_DB"," ; ",paste0(cb[j,],collapse = " ; ")),file = fileName ,append = T)
    }
    
    # Write S
    for (i in 1:nrow) {
      for (j in 1:ncol){
        write(paste0(i," ; ", j, " ; ", S[i,j]),file = fileName ,append = T)
      }
    }
    
    # Move file to input/compiled_models directory
    system(paste0("mv ", fileName, " ", wd, "/input/compiled_models"))
    
  }
}
