init.gen <- function(n_file)
{
	yini.names <- readRDS(n_file)
	y_ini <- rep(100,length(yini.names))
	names(y_ini) <- yini.names
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
