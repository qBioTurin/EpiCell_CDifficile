
FBAmodel.metadata <- function(model, model.name, wd) {
  
  source(paste0(wd, "/Rfunction/ExtractEx.R"))
  
  br = ExtractEx(model, model.name = model.name)
  
  data = rbind(data.frame(value = br$lb,
                          dir = rep("lowbnd", length(br$lb)),
                          react.id = br$react.id,
                          react.index = br$index),
               data.frame(value = br$ub,
                          dir = rep("uppbnd", length(br$ub)),
                          react.id = br$react.id,
                          react.index = br$index))
  
  all.react = data.frame(React_ID = model@react_id,
                  ReactionType = rep(".", dim(model@S)[2]),
                  React_Lb = model@lowbnd,
                  React_Ub = model@uppbnd,
                  ReactionDistrict = rep(".", dim(model@S)[2]),
                  ReactionPos = 1:dim(model@S)[2])
  
  all.react$ReactionDistrict[data$react.index] = "boundary"
  all.react$ReactionDistrict[which(all.react$ReactionDistrict != "boundary")] = "core"

  a = dplyr::filter(all.react, grepl("EX_", React_ID))
  a$ReactionType = rep("Exchange", length(a$ReactionType))

  b = dplyr::filter(all.react, grepl("DM_", React_ID))
  b$ReactionType = rep("Demand/Sink", length(b$ReactionType))

  c = dplyr::filter(all.react, grepl("sink_", React_ID))
  c$ReactionType = rep("Demand/Sink", length(c$ReactionType))

  d = dplyr::filter(all.react, grepl("tra", React_ID))
  d$ReactionType = rep("Transcription", length(d$ReactionType))

  e = dplyr::filter(all.react, grepl("rep", React_ID))
  e$ReactionType = rep("Replication", length(e$ReactionType))

  f = dplyr::filter(all.react, grepl("pbios", React_ID))
  f$ReactionType = rep("Biosynthesis", length(e$ReactionType))

  g = dplyr::filter(all.react, !grepl("biomass205|rep|tra|EX_|DM_|sink_|pbios", React_ID))
  g$ReactionType = rep("Internals/Transporters", length(f$ReactionType))

  h = dplyr::filter(all.react, grepl("biomass205", React_ID))
  h$ReactionType = rep("Objective", length(h$ReactionType))

  all.react = rbind(a, b, c, d, e, f, g, h)
  
  all.react$React_ID = gsub("\\(", replacement = "_", all.react$React_ID)
  all.react$React_ID = gsub("\\[", replacement = "_", all.react$React_ID)
  all.react$React_ID = gsub("\\)", replacement = "", all.react$React_ID)
  all.react$React_ID = gsub("\\]", replacement = "", all.react$React_ID)
  
  all.react$React_ID = gsub("_c_", replacement = "_c", all.react$React_ID)
  
  saveRDS(all.react, paste(wd, "/Input/CDmodels/", model.name, "/Reaction_metadata.rds", sep = ""))
  
}
