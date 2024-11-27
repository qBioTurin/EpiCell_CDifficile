
modelProcessing <- function(path, model.name) {
  
  #'     Inputs: 
  #'             path <- full directory path pointing the target model.mat
  #'             model.name <- the name to associated at the model
  #'    
  #'     Output <-  
  
  setwd("~/Documents/ClostridiumDiff_FBAandPN")
  
  # the generic function readMATmod imports matlab cobra models into sybil model files
  source("./code/readMATmod.R")
  
  # path del model.mat
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
  save(model.mat, file = paste("./EpiCell_CDifficile/Input/CDmodels/", model.name, "/", model.name, ".RData", sep = ""))
  
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
  save(model.mat, file = paste("./EpiCell_CDifficile/Input/CDmodels/", model.name, "/", model.name, ".RData", sep = ""))
  
  ################################################################################
  ############################# Nomenclature problem #############################
  ################################################################################
  
  # link for gene mapping, also old identifiers
  # http://bigg.ucsd.edu/api/v2/models/Recon3D/genes
  
  # Give the input file name to the function.
  # result <- rjson::fromJSON(file = "./data/Recon/Recon3D_301/Recon3D_genes.json")
  # result <- result$results
  # genes = as.matrix(do.call("rbind", result))
  
  # saving the two different nomenclatures (old vs new)
  # write.table(genes, file = "./data/gene_expression/genesID.txt", quote = F, sep = ",")
  
  # reading data for nomenclature mapping
  # genesID <- readr::read_csv("./data/gene_expression/genesID.txt")
  
}
