
flux_balance <- function(model, typename, diet) {

  ## Inputs:   model <- An object of class modelorg
  ##           typename <- A single character string indicating organism type
  ##           diet <-  A single character string indicating the the diet formulation
  ##                    "FALSE if diet is not envisaged"
  ##
  ##
  ## Outputs : An object of class optsol_optimizeProb
  ##
  
  if (diet != F) {
    
    cat("                                            ", "\n")
    cat("Diet is envisaged. Diet formulation = ", diet, "\n")
    
    d = diets[[paste(diet)]]
    
    model@react_id = gsub("\\(", replacement = "_", model@react_id)
    model@react_id = gsub("\\)", replacement = "", model@react_id)
    
    source("./Rfunction/ExtractEx.R")
    
    reactions = ExtractEx(model)
    EX = reactions@react_id
    
    if (diet == "Diet_default") {
      
      name = d[which(d[, 1] %in% EX), ][, 1]
      reacts = which(model@react_id %in% name)
      
      # conserving original bounds
      lobnd <- model@lowbnd
      upbnd <- model@uppbnd
      
      model@lowbnd[reacts] = -d[which(d[, 1] %in% EX), ][, 3]
      
    } else {
      
      name = d[which(d$Reaction %in% EX), ][, 1]
      reacts = which(model@react_id %in% unlist(name))
      
      # conserving original bounds
      lobnd <- model@lowbnd
      upbnd <- model@uppbnd
      
      model@lowbnd[reacts] =
        unlist(-(d[which(unlist(d[, 1]) %in% EX), ][, 2])*
                 ((3.14*291*2.5*1.57*6.5*13*1e-04)/1.13e-04)^-1)
    }
    
    name = d[which(d$Reaction %in% EX), ][, 1]
    reacts = which(model@react_id %in% unlist(name))
  } 
  else {
    cat("DIET is not envisaged. Bounds of EX reactions are kept as model template ...")
    cat("                                                                            ")
    cat("Initializing structure of the S4 class 'Organism' ...")
    cat("                                                                            ") 
  }
  
  #' Structure of the S4 class "Organism"
  #'
  #' Structure of the S4 class \code{Organism}
  #' representing the organisms present in the environment.
  #'
  #' @export Organism
  #' @exportClass Organism
  #' @import sybil
  #' @importFrom stats na.omit
  #' @rdname Organism
  #'
  #' @slot lbnd A numeric vector containing the lower bounds of the model structure.
  #' @slot ubnd A numeric vector containing the upper bounds of the model structure.
  #' @slot type A character vector containing the description of the organism.
  #' @slot medium A character vector containing all exchange reactions of the organism.
  #' @slot lpobj A sybil optimization object containing the linear programing problem.
  #' @slot fbasol A list with the solutions of the flux balance analysis.
  #' @slot kinetics A List containing Km and v_max values for each reactions.
  #' @slot cellarea A numeric value indicating the surface that one organism occupies (default (E.coli): 4.42 mu_m^2)
  #' @slot cellweight_mean A numeric giving the mean of starting biomass (default (E.coli): 0.489 pg)
  #' @slot model Object of class sybil::modelorg containging the genome scale metabolic model
  #' @slot algo Algorithm to be used during optimization (default fba)
  #' @slot rbiomass Name of biomass reactions which is used for growth model (set automatically but needs input if objective is not biomass optimization)
  #' @slot coupling_constraints List with coupling parameters.

  setClass("Organism",
           representation(
             lbnd="numeric",
             ubnd="numeric",
             type="character",
             medium="character",
             lpobj="sysBiolAlg",
             fbasol="list",
             kinetics="list",
             cellarea="numeric",
             cellweight_mean = "numeric",
             model="modelorg",
             algo="character",
             rbiomass="character",
             coupling_constraints="list"),
  				 prototype(kinetics = list(),
  				 					cellarea = 4.42,
  				 					cellweight_mean = 0.496))

  #' Find biomass reaction in model
  #' @description Helper function to search for biomass reaction in available reactions of a model
  #'
  #' @param model Object of class sybil::modelorg containging the genome sclae metabolic model
  #' @param keys Vector with strings which are used to find biomass reaction in model
  #' @return Vector with reaction ids for biomass reaction(s)

  findrBiomass <- function(model, keys = c("biom", "cpd11416")){

    ex <- sybil::findExchReact(model)
    ex_pos <- ex@react_pos
    ex_biom <- c(grep(paste0(keys, collapse = "|"),
    									ex@met_id, ignore.case = TRUE),
                grep(paste0(keys, collapse = "|"),
                		 sybil::met_name(model)[ex@met_pos], ignore.case = TRUE))

    rbio <- vector()

    for(k in keys){
      idx <- grep(k, sybil::react_id(model), ignore.case = TRUE)
      if(length(idx)==0) idx <- grep(k, sybil::react_name(model), ignore.case = TRUE)
      if(length(idx)>0) rbio <- c(rbio, idx)
    }

    if(length(ex_biom) > 0 ){ # take care of biomass metabolite
      rbio   <- c(rbio, ex@react_pos[ex_biom])
      ex_pos <- setdiff(ex_pos, ex@react_pos[ex_biom])
    }

    if(length(rbio)==0) return(NULL)
    rbio <- setdiff(rbio, ex_pos) # exclude exchange reactions
    return(sybil::react_id(model)[rbio])
  }

  #' @export
  #' @name Organism-constructor
  #'
  #' @param model Object of class sybil::modelorg containing the genome scale metabolic model
  #' @param algo A single character string giving the name of the algorithm to use. See \link[sybil]{SYBIL_SETTINGS}
  #' @param ex Identifier for exchange reactions
  #' @param ex_comp Defining exchange reactions whose compounds should be added to the medium of the arena (default: all)
  #' @param typename A string defining the name (set to model name in default case)
  #' @param setExInf Enable if all lower bounds of exchange reaction which are set to zero (i.e. no uptake possible!) should be set to -infitity (default: true)
  #' @param setAllExInf Enable if all lower bounds of exchange reaction should be set to -infitity (default: false)
  #' @param coupling_constraints List with coupling parameters.
  #' @return Object of class Organism

  Organism <- function(model, algo="fba", ex="EX_",
                       ex_comp=NA,
                       typename=NA, setExInf=TRUE,
                       setAllExInf=FALSE,
                       coupling_constraints=list()) {
    
    pot_biomass <- findrBiomass(model)
    
    cat("FBA model biomass-related reactions: ", pot_biomass, "\n")

    if(all(model@obj_coef==0)){
      if(length(pot_biomass)==0) stop("No objection function set in model")
      print("No objective function set, try set automatically...")
      print(paste("found:", pot_biomass))
      print(paste("set new biomass function:", pot_biomass[1]))
      sybil::obj_coef(model)[which(sybil::react_id(model)==pot_biomass[1])] <- 1
      rbiomass <- pot_biomass[1]

    } else {

      idx_obj <- which(sybil::obj_coef(model)!=0)
      idx_bio <- which(sybil::react_id(model) %in% pot_biomass)

      if(any(idx_obj %in% idx_bio)) rbiomass <- sybil::react_id(model)[idx_obj[which(idx_obj %in% idx_bio)]] # easy case: objective is (also) optimiziation of biomass

      else if(length(idx_bio)>0){ # if objective is not optimization of biomass
        cat("Optimization of biomass seems to be not the objective. Even if not optimized, a biomass reactions is needed for the growth model.\n")
        cat("Objective functions:", sybil::react_id(model)[idx_obj], "\n")
        cat("Available biomass reactions:", sybil::react_id(model)[idx_bio], "\n")
        cat("Biomass reaction used for growth model:", sybil::react_id(model)[idx_bio][1], "\n")
        rbiomass <- sybil::react_id(model)[idx_bio][1]
      } else {

        # if no biomass reaction is found
        cat("Optimization of biomass seems to be not the objective. Even if not optimized, a biomass reactions is needed for the growth model.\n")
        cat("Objective functions ID:", sybil::react_id(model)[idx_obj], "\t", "Objective functions name:", sybil::react_name(model)[idx_obj], "\n")

        if(length(idx_obj) == 0){
          stop("No objective found for the model")
        }

        rbiomass <- sybil::react_id(model)[idx_obj]
        warning("No biomass objective set. Please check that current objective makes sense.")
      }
    }

    if(algo=="coupling" & length(coupling_constraints)==0) stop("Please provide couling constraints for coupling!")

    if(is.na(ex)){

      medc <- sybil::react_id(sybil::findExchReact(model))
      names(medc) <- model@met_name[sybil::findExchReact(model)@met_pos]

    } else {

      exf <- sybil::findExchReact(model)
      medc <- exf@react_id[grep(ex, exf@react_id)]
      names(medc) <- model@met_name[exf@met_pos[grep(ex, exf@react_id)]]

    }

    if (!all(duplicated(medc)==FALSE)){

      warning("Model file contains duplicated reaction IDs. Attention, duplicated exchange reaction will be removed.")
      print(medc[which(medc==medc[duplicated(medc)])])
      medc <- unique(medc)

    }

    if(!is.na(ex_comp)){
      medc <- medc[grep(ex_comp, medc)]
    }

    rxname<- sybil::react_id(model)
    lobnd <- sybil::lowbnd(model)
    upbnd <- sybil::uppbnd(model)
    names(lobnd) = rxname

    if (setAllExInf){
      lobnd[which(names(lobnd) %in% medc ) ] <- -1000
    }

    if(setExInf){
    	# if setExInf is true then set lower bound of all
    	# exchange reactions which have zero values to -INF
      lobnd[ which(names(lobnd) %in% medc & lobnd==0) ] <- -1000
    }

    lobnd.ex <- lobnd[match(medc, rxname)]
    lobnd.ex.med <- stats::median(lobnd.ex[ which( lobnd.ex < 0 & lobnd.ex > -1000 ) ])
    cat("lobnd.ex.med", lobnd.ex.med, "\n")

    if( !is.na(lobnd.ex.med) ){
      print(paste0("Median lower bound for non-zero and non-Inf exchanges is:",
      						 round(lobnd.ex.med), 6))
      if( lobnd.ex.med > -10 ){
        warning("Many lower bounds of the model seems to be set to non -infinity. Please be aware that they will be used as maximal uptake rates even when the available medium is more abundant! (set setAllExInf=TRUE to reset all exchanges to -INF)")
        #print( lobnd.ex[ which( lobnd.ex < 0 & lobnd.ex >  lobnd.ex.med ) ] )
      }
    }

    if(is.na(typename)) typename <- ifelse(length(sybil::mod_desc(model)) > 0,
    																			 sybil::mod_desc(model),
    																			 ifelse( length(sybil::mod_name(model)) > 0,
    																			 				sybil::mod_name(model),
    																			 				ifelse( length(sybil::mod_id(model)) > 0,
    																			 								sybil::mod_id(model), "test" )))
    if(algo=="coupling"){
      lpobject <- sybil::sysBiolAlg(model,
      															algorithm="mtfEasyConstraint",
      															easyConstraint=coupling_constraints)
    } else lpobject <- sybil::sysBiolAlg(model, algorithm=algo)

    fbasol <- sybil::optimizeProb(lpobject, react=1:length(lobnd), ub=upbnd, lb=lobnd)
    names(fbasol$fluxes) = rxname
    upbnd = sybil::uppbnd(model)
    names(upbnd) = rxname

    new("Organism", lbnd=lobnd, ubnd=upbnd, type=typename, medium=medc, lpobj=lpobject,
        fbasol=fbasol, model=model, algo=algo,rbiomass=rbiomass,
        coupling_constraints=coupling_constraints)
  }

  setGeneric("lbnd", function(object){standardGeneric("lbnd")})
  setMethod("lbnd", "Organism", function(object){return(object@lbnd)})
  
  setGeneric("ubnd", function(object){standardGeneric("ubnd")})
  setMethod("ubnd", "Organism", function(object){return(object@ubnd)})
  
  setGeneric("type", function(object){standardGeneric("type")})
  setMethod("type", "Organism", function(object){return(object@type)})
  
  setGeneric("medium", function(object){standardGeneric("medium")})
  setMethod("medium", "Organism", function(object){return(object@medium)})
  
  setGeneric("lpobj", function(object){standardGeneric("lpobj")})
  setMethod("lpobj", "Organism", function(object){return(object@lpobj)})
  
  setGeneric("fbasol", function(object){standardGeneric("fbasol")})
  setMethod("fbasol", "Organism", function(object){return(object@fbasol)})
  
  setGeneric("kinetics", function(object){standardGeneric("kinetics")})
  setMethod("kinetics", "Organism", function(object){return(object@kinetics)})
  
  setGeneric("model", function(object){standardGeneric("model")})
  setMethod("model", "Organism", function(object){return(object@model)})
  
  setGeneric("cellweight_mean", function(object){standardGeneric("cellweight_mean")})
  setMethod("cellweight_mean", "Organism", function(object){return(object@cellweight_mean)})
  
  setGeneric("cellarea", function(object){standardGeneric("cellarea")})
  setMethod("cellarea", "Organism", function(object){return(object@cellarea)})
  
  setGeneric("algo", function(object){standardGeneric("algo")})
  setMethod("algo", "Organism", function(object){return(object@algo)})
  
  setGeneric("coupling_constraints", function(object){standardGeneric("coupling_constraints")})
  setMethod("coupling_constraints", "Organism", function(object){return(object@coupling_constraints)})
  
  setGeneric("optimizeWithEnvironment", 
             function(model_env, lpob=model_env@lpobj,
                      lb=model_env@lbnd, ub=model_env@ubnd,
                      cutoff=1e-6, j, sec_obj="mtf",
                      with_shadow=FALSE){standardGeneric("optimizeWithEnvironment")})
  
  setMethod("optimizeWithEnvironment", "Organism", 
            function(model_env, lpob=model_env@lpobj,
                     lb=model_env@lbnd, ub=model_env@ubnd,
                     cutoff=1e-6, j, sec_obj="mtf",
                     with_shadow=FALSE){
    
    # the LP sobject has to be updated according to the objective
    
    eval.parent(substitute(model_env@lpobj <- sybil::sysBiolAlg(model, algorithm = "mtf")))

    fbasl <- sybil::optimizeProb(lpob, react=1:length(lb),
                                 ub=ub, lb=lb,
                                 resetChanges = FALSE)

    mtf_sol <- sybil::optimizeProb(model, algorithm="mtf", mtfobj=fbasl$obj,
                                   gene = NULL, react = NULL, lb = NULL, ub = NULL,
                                   obj_coef = NULL, fldind = TRUE,
                                   retOptSol = T);
    return(mtf_sol)

  })

  model_env = Organism(model, algo = "fba", ex = "EX_", ex_comp = "all", typename = typename)

  mtf_sol = optimizeWithEnvironment(model_env)
  mtf_sol = list(mtf_sol, model@lowbnd, model@uppbnd)
  
  # v = sybil::getFluxDist(mtf_sol[[1]])
  
  cat("solver: ", mtf_sol[[1]]@solver, "\n")
  cat("number of variables: ", model.mat@react_num, "\n")
  cat("number of constraints: ", model.mat@met_num, "\n")
  cat("value of objective function (mtf): ", mtf_sol[[1]]@alg_par[["wtobj"]], "\n")
  cat("how to obtain flux distribution: {r} v = sybil::getFluxDist(mtf_sol[[1]])", "\n")
  
  return(mtf_sol)

}
