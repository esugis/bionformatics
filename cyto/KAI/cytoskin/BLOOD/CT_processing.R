setwd("~/Documents/KAI/cytoskin/CT_blood2")
files_csv=list.files(path = "~/Documents/KAI/cytoskin/CT_blood2", pattern ="*.csv", all.files = FALSE,
           full.names = FALSE, recursive = FALSE,
           ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)

#substitute "(" and ")" and "-" in the files from /CT_blood2/ to "_"
for (i in 1:length(files_csv)){
  file=files_csv[i]
fname=files_csv[i]
fname=gsub("\\)","_",fname)
fname=gsub("\\(","_",fname)
fname=gsub("-","_",fname)
system(sprintf("mv %s %s ", paste("'",as.character(file),"'",collapse=", ",sep=""), fname))
}

