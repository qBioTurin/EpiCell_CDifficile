
# GeneExpr must be a dtaframe with two columns named as "gene" and "val"

Eval.ReactionExpression = function(txt, GeneExpr){
  
  txt1 = gsub(" ", replacement = "", txt )
  txt1 = gsub("\\(",replacement = "[", txt1 )
  txt1 = gsub("\\)",replacement = "]", txt1 )
  
  txtsep <- base::strsplit(x = txt1,
                           split = paste0("(?<=.)(?=", "\\[|\\]", ")"),
                           perl = TRUE)
  
  txtsep = txtsep[[1]]
  txtsep <- base::strsplit(x = txtsep,
                           split = paste0("(?<=", "\\[|\\]", ")"),
                           perl = TRUE)
  txtaux = unlist(txtsep)
  txtaux = c("[",txtaux,"]")
  ANDind = grep("(and)|(AND)",txtaux)
  
  if(length(ANDind) > 0 ){
    # if(length(grep("(\\[)|(\\])",txtaux) )== 0)
    #   txtaux = c("[",txtaux,"]")
    
    start = (grep("\\[",txtaux)) # the first (
    end = (grep("\\]",txtaux))
    startChar = rep(1,length(start)) # open
    endChar = rep(2,length(end)) # closed
    Brackets = data.frame(type = c(startChar,endChar), Index = c(start,end))
    Brackets.tmp =  Brackets[order(Brackets$Index),]
    
    while(length(Brackets.tmp[,1])>0){
      Brackets.tmp$Conseq = c(0,diff(Brackets.tmp$type))
      conseqIndex = which(Brackets.tmp$Conseq == 1)
      ClB = Brackets.tmp$Index[conseqIndex]
      OpB = Brackets.tmp$Index[conseqIndex-1]
      ANDbtwBr =do.call("rbind",lapply(1:length(ClB), function(x){
        ind = which( ANDind < ClB[x] & ANDind > OpB[x])
        if(length(ind) >0)
          return(data.frame(ClB = ClB[x],OpB =OpB[x], ind = ind))
      } ))
      
      if( length( ANDbtwBr )>0 ){
        txtaux[ANDbtwBr$OpB ] = "min(["
        txtaux[ANDbtwBr$ClB ] = "])"
        ANDind = ANDind[-c(ANDbtwBr$ind) ] 
      }
      
      Brackets.tmp = Brackets.tmp[-c(conseqIndex,conseqIndex-1),]
      
    }
  }
  
  txtMin = paste0(txtaux,collapse = "")
  txtFinal = gsub("(or)|(OR)",replacement = " + ", txtMin )  
  txtFinal = gsub("(and)|(AND)",replacement = " , ", txtFinal )  
  txtFinal = gsub("\\[",replacement = "c(", txtFinal )
  txtFinal = gsub("\\]",replacement = ")", txtFinal )
  
  gsub(pattern = "\\.",replacement = "_",x = txtFinal) -> txtNum
  gsub(pattern = "(min)|(c)|(\\()|(\\))|( )",replacement = "",x = txtNum) -> txtNum2
  
  Gend = grep(x = unlist(strsplit(gsub("(,)|(\\+)","~", txtNum2), "~")),
              pattern = "[0-9]+",value = T)
  
  ggExpr = gsub(pattern = "\\.", replacement = "_", x = GeneExpr$gene)
  
  for(g in Gend) {
    
    print(paste0("#### g | ", g))
    print(paste0("#### replacement | ", GeneExpr$val[which(ggExpr == g)]))
    
    txtNum = gsub(paste0("\\b",g,"\\b"), 
                  replacement = GeneExpr$val[which(ggExpr == g)], 
                  txtNum)
    
    }
  
  #GSE = data.frame(Val = eval(parse(text=txtNum)) , ExpNum = txtNum, ExpGenes =txtFinal )
  #return(GSE)
  return( eval(parse(text=txtNum)) )
}


###### Example 
# txt = "( (527_AT1 or 8992_AT1) and 50617_AT1 and 533_AT1 and (527_AT1 or 8992_AT1) and 245972_AT1) OR (50617_AT1 and 533_AT1 and 527_AT1 and 8992_AT1 and 9114_AT1)"
# GeneExpr = data.frame(gene = c("527_AT1","8992_AT1","50617_AT1" , "533_AT1" ,"527_AT1" , "8992_AT1" ,"245972_AT1" , "50617_AT1" , "533_AT1" , "527_AT1", "8992_AT1" , "9114_AT1"),
#                       lev = 4 )
# Eval.ReactionExpression(txt,GeneExpr)
