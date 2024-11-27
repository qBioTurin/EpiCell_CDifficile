
ExtractEx <- function(model, model.name, bigg.path) {
  
  StoichM = model@S
  
  # identifying columns with only one entry
  oneEntry = Matrix::colSums(StoichM != 0) == 1
  
  if (sum(oneEntry) > 0) {
    
    # exchange reactions can be with a -1 or 1
    ExReact = (Matrix::colSums(StoichM[, oneEntry, drop = FALSE] == 1) == 1 | 
                 Matrix::colSums(StoichM[, oneEntry, drop = FALSE] == -1) == 1)
    
    # finding Ex_reaction's IDs
    ex = c(1:dim(StoichM)[2])[oneEntry[ExReact]]
    
    # extracting metabolites involved in exchange reactions
    Met = which(Matrix::rowSums(
      abs(StoichM[, ex, drop = FALSE])) > 0)[StoichM[which(
        Matrix::rowSums(abs(StoichM[, ex, drop = FALSE])) > 0), 
        ex, drop = FALSE]+1]
    
    BIGGdata = read.delim2(bigg.path)
    
    reactions = data.frame(index = ex,
                           react.id = model@react_id[ex],
                           lb = model@lowbnd[ex], 
                           ub = model@uppbnd[ex],
                           react.eq = BIGGdata[ex, ]$equation
    )
    
  } else {
    
    warning("model without exchage ractions")
    reactions <- NULL
    
  }
  
  reactions$react.id = gsub("\\(", replacement = "_", reactions$react.id)
  reactions$react.id = gsub("\\[", replacement = "_", reactions$react.id)
  reactions$react.id = gsub("\\)", replacement = "", reactions$react.id)
  reactions$react.id = gsub("\\]", replacement = "", reactions$react.id)
  
  reactions$react.id = gsub("_c_", replacement = "_c", reactions$react.id)
  
  return(reactions)
  
}
