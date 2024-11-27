
params <- readRDS("/home/docker/data/Results/FVA_Ablated_NoDrug_Model_Analysis/params_Model_sensitivity.RDS")

# setwd("/home/")
# folder_trace = paste0("/home/",basename(params$folder_trace) )

folder_trace = "/home/docker/data/Results/CDiff_Ablated_NoDrug_Model_Analysis/"
folder_sensitivity = "/home/docker/data/Results/FVA_Ablated_NoDrug_Model_Analysis/"

flux_fname_file = paste0(folder_sensitivity,"/flux_fname")
fls = list.files(folder_trace,pattern = ".flux$")
fva_name = gsub(fls,pattern = ".flux$",replacement = "")
fva_name = gsub(fva_name,pattern = "analysis",replacement = "fva")
fva_name = paste0(fva_name,"-")
names(fva_name) = fls

####### ATTENZIONE
## Per piu' problemi di FBA, bisogna leggere il rispettivo file e
## filtrare le colonne in base ai flussi presenti nei fbafile, questo perche'
## VARIABILITY.sh funziona solo su singolo fbafile!!!!

fls; fl = fls

chk_dir <- function(path){
  pwd <- basename(path)
  return(paste0(file.path(dirname(path),pwd, fsep = .Platform$file.sep), .Platform$file.sep))
}

t = "MIN"

fba_files_index = 1 + as.numeric(gsub(fl,pattern = paste0("(",params$out_fname_analysis,"-[0-9]-)|(.flux)"),replacement = ""))
flux <- read.csv(file = paste0(folder_trace,"/",fl), sep = "", header = F)
flID = flux[1,]
flux_fname =  params$flux_fname[params$flux_fname %in% flID]

dfList = lapply(flux_fname, function(f){
  initFile = paste0(folder_sensitivity, fva_name[fl],f,"_",t)
  df = read.table( initFile,sep = " ",header = T )
  finalFile = gsub(pattern = "\\(|\\)|\\[|\\]| ",replacement = "_",initFile)
  file.rename(from = initFile,
              to = finalFile)
  colnames(df) = c("Time",f)
  return(df)
})
