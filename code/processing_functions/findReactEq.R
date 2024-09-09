findReactEq <- function(model, reagent, metadata) {
  
  # Returns a list of equations (formulas) of reactions in which the 
  # metabolite (reagent) is involved
  #
  # USAGE:
  #
  #   findReactEq(model, reagent)
  #
  # INPUTS:
  #    model:             Model structure
  #    reagent:           Metabolite
  #    metadata:          Reaction/Metabolite and equations data
  #
  # OUTPUTS:
  #     List of reactions + Reaction formulas corresponding and indexing

  # reactions nomenclature in accordance with VMH
  ReactionsNames = unlist(model@react_id)
  # reactant's nomenclature in accordance with VMH
  ReagentsNames = unlist(model@met_id)
  
  S = as.matrix(model@S)
  ReagentsIndex = which(ReagentsNames == reagent)
  ReactionIndex = which(S[ReagentsIndex, ] != 0)
  reactList = ReactionsNames[ReactionIndex]
  
  ReactEq = metadata[which(metadata$abbreviation %in% reactList), c(1, 2, 3, 4)]
  # sybil::printReaction(model, react = reactList, printOut = TRUE)
  
  return(ReactEq)
}
