
################################################
# Function: metadataEQ
#
# 
# 

metadataEQ = function(model, model.name, 
                      prefix, suffix, extMetFlag = "b",
                      fielddelim = "\t", entrydelim = ", ",
                      makeClosedNetwork = FALSE,
                      onlyReactionList,
                      minimalSet,
                      fpath = SYBIL_SETTINGS("PATH_TO_MODEL")) {
  
  # validate model structure before writing
  validObject(model)
  
  # filenames
  if (missing(prefix)) {
    prefix = gsub("\\s+", "_", model.name)
  }
  
  if (missing(suffix)) {
    suffix = switch(fielddelim,
                     "\t" = { "tsv" },
                     ";"  = { "csv" },
                     ","  = { "csv" },
                     "|"  = { "dsv" },
                     { "dsv" }
    )
  }
  
  fnameR = paste(paste(prefix, "react", sep = "_"), suffix, sep = ".")
  fnameM = paste(paste(prefix, "met",   sep = "_"), suffix, sep = ".")
  fnameD = paste(paste(prefix, "desc",  sep = "_"), suffix, sep = ".")
  
  # path to output file
  tsvfileR = file.path(fpath, fnameR)
  tsvfileM = file.path(fpath, fnameM)
  tsvfileD = file.path(fpath, fnameD)
  
  #--------------------------------------------------------------------------#
  # reactions list
  #--------------------------------------------------------------------------#
  
  # create reaction strings
  rstr = GenerateReactionEquations(model,
                                   makeClosedNetwork,
                                   entrydelim,
                                   extMetFlag)
  
  met.id = model@met_id
  
  if( sum(unique(stringr::str_sub(met.id, -1)) == "]") > 1) {
    met.id = gsub("\\[",replacement = "_", met.id)
    met.id = gsub("\\]",replacement = "", met.id)
  }
  
  mod_compart = unique(stringr::str_sub(met.id, -1))
  
  met_comp = c()
  
  for (i in 1:length(met.id)) {
    comp = stringr::str_sub(met.id[i], -1)
    met_comp[i] = which(comp == mod_compart)
    }
  
  if (isTRUE(onlyReactionList)) {
    write.table(x = data.frame(equation = rstr$equat),
                row.names = FALSE, 
                file = tsvfileR, 
                sep = fielddelim)
  } else if (isTRUE(minimalSet)) {
    write.table(x = data.frame(abbreviation = model@react_id,
                               equation = rstr$equat,
                               lowbnd = model@lowbnd,
                               uppbnd =  model@uppbnd,
                               obj_coef = model@obj_coef),
                row.names = FALSE, 
                file = tsvfileR, 
                sep = fielddelim)
    } 
  else {
    
    write.table(x = data.frame(
      abbreviation = model@react_id,
      equation     = rstr$equat,
      reversible   = rstr$revers,
      compartment  = rstr$compat,
      lowbnd       = model@lowbnd,
      uppbnd       = model@uppbnd,
      obj_coef     = model@obj_coef),
      row.names = FALSE, 
      file = tsvfileR, 
      sep = fielddelim)
  }
  
  #--------------------------------------------------------------------------#
  # metabolites list
  #--------------------------------------------------------------------------#
  
  if ( (!isTRUE(onlyReactionList)) && (!isTRUE(minimalSet)) ) {
    
    metunq <- sort(unique(rstr$metab))
    mpos <- lapply(metunq, function(x) rstr$metab %in% x)
    
    # metabolite names
    metNames <- lapply(mpos, function(x) unique(model@met_id[x]))
    metNames <- unlist(lapply(metNames, paste, collapse = entrydelim))
    
    # metabolite compartments
    metCompart = lapply(mpos, function(x) mod_compart[met_comp[x]])
    metCompart = unlist(lapply(metCompart, paste, collapse = entrydelim))
    
    write.table(x = data.frame(
      abbreviation = metunq,
      name         = metNames,
      compartment  = metCompart),
      row.names = FALSE, file = tsvfileM, sep = fielddelim)
    
  }
  
  #--------------------------------------------------------------------------#
  # model description
  #--------------------------------------------------------------------------#
  
  if ((!isTRUE(onlyReactionList)) && (!isTRUE(minimalSet))) {
    
    # get id's of metabolites in different compartments
    # (one per compartment)
    metDiffComp = match(mod_compart, mod_compart[met_comp])
    metAbbrevComp = character(length(metDiffComp))
    
    # get the compartment abbreviations
    metALl = grepl("^.+(\\[\\w+\\])$", model@met_id[metDiffComp])
    metAbbrevComp[metALl] = sub("^.+(\\[\\w+\\])$", "\\1", model@met_id[metDiffComp[metALl]])
    metAbbrevComp[!metALl] = mod_compart[!metALl]
    
    # generate output format
    ma = paste(metAbbrevComp, collapse = entrydelim)
    mc = paste(mod_compart, collapse = entrydelim)
    
    write.table(x = data.frame(
      name         = model.name,
      id           = "FBAmodel",
      description  = model.name,
      compartment  = mc,
      abbreviation = ma,
      Nmetabolites = dim(model@S)[1],
      Nreactions   = dim(model@S)[2]),
      row.names = FALSE, file = tsvfileD, sep = fielddelim)
    
    }
  
  #--------------------------------------------------------------------------#
  # end
  #--------------------------------------------------------------------------#
  
  return(TRUE)
  
}