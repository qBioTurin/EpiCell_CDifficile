
ExtractEx <- function(model) {
  
  StoichM <- model@S
  
  # identifying columns with only one entry
  oneEntry <- Matrix::colSums(StoichM != 0) == 1
  
  if (sum(oneEntry) > 0) {
    
    # exchange reactions can be with a -1 or 1
    ExReact = (Matrix::colSums(StoichM[, oneEntry, 
                                       drop = FALSE] == 1) == 1 | 
                 Matrix::colSums(StoichM[, oneEntry, 
                                          drop = FALSE] == -1) == 1)
    
    # finding Ex_reaction's IDs
    ex = c(1:dim(StoichM)[2])[oneEntry[ExReact]]
    
    # extracting metabolites involved in exchange reactions
    Met = which(Matrix::rowSums(
      abs(StoichM[, ex, drop = FALSE])) > 0)[StoichM[which(
        Matrix::rowSums(abs(StoichM[, ex, drop = FALSE])) > 0), 
        ex, drop = FALSE]@i+1]
    
    reactions = new("reactId_Exch", mod_id  = model@mod_id,
                    mod_key = model@mod_key, rpnt = ex, rid = model@react_id[ex],
                    # finding uptake reactions, uptake reactions are those reactions
                    # will proceed to the left (see ReadMe.R)
                    upt = model@lowbnd[ex] < 0, mpnt = Met, mid = model@met_id[Met],
                    lb = model@lowbnd[ex], ub = model@uppbnd[ex])
  } else {
    
    warning("model without exchage ractions")
    reactions <- NULL
    
  }
  
  reactions@react_id = reactions@react_id[which(grepl("EX_", reactions@react_id, fixed = TRUE))]
  
  reactions@react_id = gsub("\\(", replacement = "_", reactions@react_id)
  reactions@react_id = gsub("\\)", replacement = "", reactions@react_id)
    
  return(reactions)
}
