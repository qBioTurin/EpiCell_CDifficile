
################################################
# Function: GenerateReactionEquations
#
#
#

GenerateReactionEquations <- function(model,
                                      makeClosedNetwork,
                                      entrydelim = ", ",
                                      extMetFlag = "b") {
  
  # required data structures
  equat  <- vector(mode = "character", length = dim(model@S)[2])
  compat <- vector(mode = "character", length = dim(model@S)[2])
  
  react_rev = c()
  
  b = cbind(lb = model@lowbnd, ub = model@uppbnd)
  
  for(i in 1:nrow(b)) {
    react_rev[i] = !any(b[i, ] == 0)
  }
  
  revers = ifelse(react_rev, "Reversible", "Irreversible")
  arrow = ifelse(react_rev, "<==>", "-->")
  
  # remove compartment flag if existing
  metab = sub("\\[\\w+\\]$", "", model@met_id)
  
  met.id = model@met_id
  
  if( sum(unique(stringr::str_sub(met.id, -1)) == "]") >= 1 ) {
    met.id = gsub("\\[",replacement = "_", met.id)
    met.id = gsub("\\]",replacement = "", met.id)
  }
  
  metcp = stringr::str_sub(met.id, -1)
  mod_compart = unique(stringr::str_sub(met.id, -1))
  
  met_comp = c()
  
  for (i in 1:length(met.id)) {
    comp = stringr::str_sub(met.id[i], -1)
    met_comp[i] = which(comp == mod_compart)
  }
  
  for (j in 1:dim(model@S)[2]) {
    
    column = model@S[, j]
    # row indices
    constr_ind = which(column != 0)
    # stoichiometric coefficients
    stcoef = column[constr_ind]
    
    # check if reaction is empty
    if (length(constr_ind) > 0) {
      
      comp = unique(mod_compart[met_comp[constr_ind]])
      
      # reaction involves more than one compartment
      if (length(comp) > 1) {
        if (metcp[constr_ind[1]] == metab[constr_ind[1]]) {
          metr = paste(model@met_id[constr_ind], 
                       "_", mod_compart[met_comp[constr_ind]], sep = "")
        } 
        else {
          metr = met.id[constr_ind]
        }
        compat[j] = paste(comp, collapse = entrydelim)
        compflag  = ""
      } 
      else {
        
        # Check if the current reaction is an exchange reaction.
        # In order to build a closed network, we need to add a 'boundary'
        # metabolite [b].
        
        if ((isTRUE(makeClosedNetwork)) && (length(constr_ind) == 1)) {
          if (metcp[constr_ind[1]] == metab[constr_ind[1]]) {
            metIN = paste(met.id[constr_ind], "_",
                          mod_compart[met_comp[constr_ind[1]]], sep = "")
          }
          else {
            metIN = met.id[constr_ind]
          }
          metr = c(metIN, paste(metab[constr_ind], "_", extMetFlag, sep = ""))
          constr_ind = c(constr_ind, constr_ind)
          stcoef =  c(stcoef, (stcoef * -1))
          compat[j] = comp
          compflag = ""
        }
        else {
          metr = metab[constr_ind]
          compat[j] <- comp
          # if yes, the metabolite id does not contain the metabolite compartment
          if (metcp[constr_ind[1]] == metab[constr_ind[1]]) {
            compflag  <- paste0("_", mod_compart[met_comp[constr_ind[1]]])
          }
          else {
            compflag  <- metcp[constr_ind[1]]
          }
          metr = paste0(metab[constr_ind], "_", comp)
          metr = strex::str_singleize(metr, paste0( "_", comp))
        }
      }
      
      educt <- vector(mode = "list")
      product <- vector(mode = "list")
      
      for (i in seq(along = constr_ind)) {
        if (stcoef[i] > 0) {
          stoich <- ifelse(stcoef[i] != 1,
                           paste("(", stcoef[i], ") ", sep = ""),
                           "")
          product[metr[i]] <- paste(stoich, metr[i], sep = "")
        }
        else {
          stoich <- ifelse(stcoef[i] != -1,
                           paste("(", (stcoef[i] * -1), ") ", sep = ""),
                           "")
          educt[metr[i]] <- paste(stoich, metr[i], sep = "")
        }
      }
      
      equattmp = paste(paste(educt, collapse = " + "),
                       arrow[j],
                       paste(product, collapse = " + "))
      
      equat[j] = equattmp
      
    }
  }
  
  return(list(equat  = equat,
              compat = compat,
              revers = revers,
              metab  = metab))
  
}
