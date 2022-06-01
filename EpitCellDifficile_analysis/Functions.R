init.gen <- function(n_file)
{
	yini.names <- readRDS(n_file)
	y_ini <- c(
		# IECs
		4.57e05,
		# Damage
		0,
		# BiomassCD
		4.57e05*100*(1/3)*(((pi*(0.5)^2)/4)*(((4*0.5)/6) + 5.5))*1.3e-12,
		# CD
		4.57e05*100,
		# Drug
		0.17116*6.022e20*0.15e-06,
		# pheme_e
		3.5e-04*6.022e20,
		# pro_L_e
		0.36*6.022e20,
		# leu_L_e
		0.57*6.022e20,
		# ile_L_e
		0.32*6.022e20,
		# val_L_e
		0.58*6.022e20,
		# trp_L_e
		0.098*6.022e20,
		# cys_L_e
		0.25*6.022e20,
		0, 0, 0, 0, 0, 0)

	names(y_ini) <- c("IECs", "Damage", "BiomassCD", "CD", "pheme_e", "Drug",
										"pro_L_e", "leu_L_e", "ile_L_e", "val_L_e", "trp_L_e", "cys_L_e",
										"pro_L_v", "leu_L_v", "ile_L_v", "val_L_v", "trp_L_v", "cys_L_v")

	y_ini = y_ini[yini.names]
	return(y_ini)
}

FBA.generation = function(model,write=F,load=F){

	if(load)
	{
		load(model)
		FinalMatrix = FinalMatrix$X1
		return(matrix(FinalMatrix,ncol=1))
	}else{

		# Problem: S*x = b -. size x = n and size of S = m X n
		# With this function the file to pass to the GLPKsolver is generated:
		# 1 row)    n_row n_col GLP_MAX (or GLP_MIN if the obj has to max or min)
		# 2 row)    the coeff for defining the obj function
		#           (the num of coeff has to be == length of x)
		# m rows)   definition of the S row bounds (i.e. b)
		# n rows)   definition of the x bounds
		# m*n rows) definition of the S coeffs: row_index col_index value
		x = load(model)
		Matlab.file = get(x)

		#### Start writing the csv file to pass to GLPK solver
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

		if(write){
			fileName = "FBAmodel"
			write(paste0(nrow," ; ",ncol," ; ","GLP_MAX" ),file = fileName )
			write(paste(obj,collapse = " "),file = fileName,append = T )
			for(i in 1:nrow){
				write(paste0("GLP_F"," ; ",paste0(rb[i,],collapse = " ; ") ),file = fileName ,append = T)
			}
			for(j in 1:ncol){
				write(paste0("GLP_DB"," ; ",paste0(cb[j,],collapse = " ; ")),file = fileName ,append = T)
			}
			for(i in 1:nrow){
				for(j in 1:ncol){
					write(paste0(i," ; ", j, " ; ", S[i,j]),file = fileName ,append = T)
				}
			}
			FinalMatrix = readr::read_csv("FBAmodel", col_names = FALSE)
			save(FinalMatrix,file = "Input/FinalMatrix.RData")

			return(0)
		}else{
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
}

saveNames= function(model){

	x = load(model)
	Matlab.file = get(x)

	ReactionsNames <- unlist(Matlab.file@react_id)

	return(matrix(ReactionsNames,ncol = 1))
}

saveMMConstant = function(type) {

	met.name = c("EX_pheme_e_in", "EX_pro_L_e_in", "EX_leu_L_e_in", "EX_ile_L_e_in",
							 "EX_val_L_e_in", "EX_trp_L_e_in", "EX_cys_L_e_in")

	Vmax = c(0.01,
					 0.01, 0.01, 0.01, 0.01, 0.01, 0.01)
	KM = c(0.005,
				 0.005, 0.005, 0.005, 0.005, 0.005, 0.005)

	if(type == "KM"){
		y = data.frame(met.name, KM)
	}else{
		y = data.frame(met.name, Vmax)
	}

	return(y)
}

InflammationFunction = function(gamma) {

	Dv = 0.03 #[mm]
	Vrbc = 1.1e04 #[mm/h]
	RHOrbc = 5 # [cell/muL]
	HeamRBC = 5e09 # [molecule]
	Aw = 1.13e-04 # [m^2]
	W = 0.1 # [mm]
	L = 1 # [mm]
	EM = 13 # [unit]

	Inflam = gamma*pi*(Dv/2)^2*(Vrbc*RHOrbc*HeamRBC*Aw/2*pi*(W/2)*L*EM)
	return(Inflam)

}

ComputeBiomassBac = function(diameter, len) {

	mass = (1/3)*(((pi*(diameter)^2)/4)*(((4*diameter)/6) + len))*1.3e-12

	return(mass)
}

EvalDiet = function(diet, ex) {

	Lin = 291 # [cm]
	Din = 2.5 # [cm]
	P = 1.57 # [unit]
	LM = 6.5 # [unit]
	EM = 13 # [unit]

	vheme = 1.14 # [g/(day*person)]
	Na = 6.022e20 # [molecule]
	Aw = 1.13e-04 # [m^2]
	Mpheme = 0.616487 # [g/mmol]

	Ain = pi*Lin*Din*P*LM*EM*1e-04 # [m^2]

	if (is.null(diet)) {

		vD = (vheme*Na)/(Mpheme*24)*((Ain/Aw)^(-1))

	} else {

		# diet = readr::read_delim(paste0("./Input/", diet, ".tsv", sep=""), "\t", escape_double = FALSE, trim_ws = TRUE)

		# diet[["Reaction"]] = gsub("\\(", replacement = "_", diet[["Reaction"]])
		# diet[["Reaction"]] = gsub("\\)", replacement = "", diet[["Reaction"]])

		# diet[["Reaction"]] = gsub("\\[", replacement = "_", diet[["Reaction"]])
		# diet[["Reaction"]] = gsub("\\]", replacement = "", diet[["Reaction"]])

		# diet[["Flux Value"]][c("EX_cys_L_e", "EX_trp_L_e", "EX_val_L_e", "EX_ile_L_e", "EX_leu_L_e", "EX_pro_L_e") %in% ex]

		fluxes = data.frame(Reaction = c("EX_cys_L_e", "EX_trp_L_e", "EX_val_L_e",
																		 "EX_ile_L_e", "EX_leu_L_e", "EX_pro_L_e"),
												Flux = c(7.511832, 28.921675, 46.887810,
																 53.815279, 4.554825, 37.784377))

		vD = ((fluxes$Flux[fluxes$Reaction == ex]*Na)/24)*((Ain/Aw)^(-1))
	}
	return(vD)
}

EvalTransport = function(yield) {

	IECs_t0 = 4.57e05
	vT = yield/IECs_t0
}

EvalTreat = function(dose){

	Mmtz = 0.17116 # [g/mmol]
	Na = 6.022e20 # [molecule]
	Aw = 1.13e-04 # [m^2]

	Lin = 291 # [cm]
	Din = 2.5 # [cm]
	P = 1.57 # [unit]
	LM = 6.5 # [unit]
	EM = 13 # [unit]

	Ain = pi*Lin*Din*P*LM*EM*1e-04 # [m^2]

	ther = (1/2)*((Mmtz*Na)/dose)*((Ain/Aw)^(-1))

	return(ther)

}
