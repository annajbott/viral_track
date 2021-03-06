## pipeline_viral_track R module file
## Containing R functions etc. then sourced into R scripts

## Auxiliary functions
## Function to extract information from STAR Log output

Extraction_Log_final = function(path_to_Log_file) {
  #Loading the final log file from STAR
  Log_file =read.delim(path_to_Log_file,header = F,sep="\t")
  Name_variable = as.character(Log_file$V1)
  Value_variable = as.character(Log_file$V2)
  
  #Extracting the informations about the mapping quality 
  Uniquely_mapped_percent = Value_variable[Name_variable=="                        Uniquely mapped reads % |"]
  Multiple_mapped_percent = Value_variable[Name_variable=="             % of reads mapped to multiple loci |"]
  Unmapped_too_short_percent = Value_variable[Name_variable=="                 % of reads unmapped: too short |"]
  Unmapped_mismatch_percent = Value_variable[Name_variable=="       % of reads unmapped: too many mismatches |"]
  Unmapped_other_percent = Value_variable[Name_variable=="                     % of reads unmapped: other |"]
  
  remove_percent = function(x) {
    l = nchar(x)
    x = substr(x,1,l-1)
    x=as.numeric(x)
    return(x)
  }  
  
  Uniquely_mapped_percent=remove_percent(Uniquely_mapped_percent)
  Multiple_mapped_percent=remove_percent(Multiple_mapped_percent)
  Unmapped_too_short_percent=remove_percent(Unmapped_too_short_percent)
  Unmapped_mismatch_percent = remove_percent(Unmapped_mismatch_percent)
  Unmapped_other_percent = remove_percent(Unmapped_other_percent)
  Total_unmapped = Unmapped_too_short_percent + Unmapped_mismatch_percent + Unmapped_other_percent
  
  Matrix_mapping = matrix(c(Uniquely_mapped_percent,Multiple_mapped_percent,Total_unmapped),nrow = 3)
  
  #####Looking at length of mapping, deletion and insertion
  
  Mean_mapped_length = as.numeric(Value_variable[Name_variable=="                          Average mapped length |"])
  Mean_deletion_length = as.numeric(Value_variable[Name_variable=="                        Deletion average length |"])
  Mean_insertion_length = as.numeric(Value_variable[Name_variable=="                       Insertion average length |"])
  #####Looking at rates of  of mismatch, insertion and deletion
  
  Mismatch_rate= remove_percent(Value_variable[Name_variable=="                      Mismatch rate per base, % |"])
  Deletion_rate = remove_percent(Value_variable[Name_variable=="                         Deletion rate per base |"])
  Insertion_rate = remove_percent(Value_variable[Name_variable=="                        Insertion rate per base |"])
  
  List_elements = list(Mapping_result= Matrix_mapping , 
                       Length_vector = c(Mean_mapped_length,Mean_insertion_length,Mean_deletion_length),
                       Rate_vector = c(Mismatch_rate,Insertion_rate,Deletion_rate) )
  
  return(List_elements)
}

## ---------------------------------------------------------------
## Function to compute similarity between genome sequences

alignement_score <- function(x) { # x is a vector of sequences in the Biostring format (DNAstring format)
  dist_matrix = matrix(0,nrow = length(x),ncol = length(x))
  comparison_sequence = c()
  for (i in 1:nrow(dist_matrix)){
    for (j in 1:ncol(dist_matrix)) {
      dist_matrix[i,j] = pairwiseAlignment(pattern = x[i], subject = x[j], type = "local-global", substitutionMatrix = nucleotideSubstitutionMatrix(match = 1, mismatch = 0, baseOnly = FALSE, type = "DNA") , gapOpening = 1, gapExtension = 0, scoreOnly=TRUE)
    }
  }
  dist_matrix = -dist_matrix
  rownames(dist_matrix) = names(x)
  colnames(dist_matrix) = names(x)
  return(dist_matrix)
}

string.to.colors = function (string, colors = NULL) 
{
  if (is.factor(string)) {
    string = as.character(string)
  }
  if (!is.null(colors)) {
    if (length(colors) != length(unique(string))) {
      (break)("The number of colors must be equal to the number of unique elements.")
    }
    else {
      conv = cbind(unique(string), colors)
    }
  }
  else {
    conv = cbind(unique(string), rainbow(length(unique(string))))
  }
  unlist(lapply(string, FUN = function(x) {
    conv[which(conv[, 1] == x), 2]
  }))
}
