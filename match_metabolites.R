#rt library matching 
#steps 
library(XML)

#read in the retention time library which has been built by researchers, saved in metabolights 
rt_lib <- read.csv("E:/748/Assessment 3/match_rt/rt_library.csv", header = TRUE, fill = TRUE)
head(rt_lib)

#for matches we need only the metabolite name and the retention time if a match is found we can report the other meta data later on aswell 
match_rt_lib <- rt_lib[, c("Number","metabolite_identification", "retention_time", "CAS_Number", "Formula", "Isomeric_SMILES", "InChI")]
head(match_rt_lib)

#read in the experimental data from a tumourigenic breast epithelial cell sample
exp_data <- read.csv("E:/748/Assessment 3/match_rt/experimental_rt.csv")
head(exp_data)

#make a new dataframe with only the id and the retention time 
match_exp_data <- exp_data[,c("id","rt")]

head(match_exp_data)
#read in the library 

#allow for variance in the rt found for the metabolite in the experiment 
#for 2 seconds and find the minimum and max the rt which this metabolite could reasonably have 
#make a new df including columns for the minimum and maximum retention time which could be had 
match_exp_data_w_rt_range <- cbind(match_exp_data, min_rt = match_exp_data[,"rt"] - 2, max_rt = match_exp_data[,"rt"] + 2)


#make a data frame for the matched metabolites 
matched_metabolites <- data.frame(
  exp_ids_that_matched_a_lib_met_rt = I(list()),
  rt = numeric(0),
  matching_metabolite = character(0),
  CAS_Number = character(0),
  Formula = character(0),
  Isomeric_SMILES = character(0),
  InChI = character(0),
  stringsAsFactors = FALSE
)

#loop through the entries in the rt library, compare to each rt range 
for(lib_row in 1:nrow(match_rt_lib)){
  print(c("looking for a match to", "NUMBER:", match_rt_lib[lib_row, "Number"], "METABOLITE:", match_rt_lib[lib_row, "metabolite_identification"]))
  #find the rt of this entry in the library   
  library_rt = (match_rt_lib[lib_row, "retention_time"])
  #make it numeric if the retention time for this metabolite in the library is known 
  if(library_rt != "n.d"){
    library_rt = as.numeric(library_rt)
  }
  
  if(library_rt != "n.d"){
    #loop through the entries in the experimental data
    for(exp_row in 1:nrow(match_exp_data_w_rt_range)){
      
      
      #assign the minimum rt of the metabolite of this iteration to a variable 
      min_exp_rt = as.numeric(match_exp_data_w_rt_range[exp_row, "min_rt"])
      
      
      #assign the maximum rt of the metabolite of this iteration to a variable 
      max_exp_rt = match_exp_data_w_rt_range[exp_row, "max_rt"]
      
      
      
      #if the library rt is within the range of possible rt of this metabolite in the 
      #experimental data then assign it to a column of matches
      #
      
      if (library_rt > min_exp_rt && library_rt < max_exp_rt){
        
        
        # get the experimental id and rt value
        exp_id <- match_exp_data_w_rt_range[exp_row, "id"]
        exp_rt <- match_exp_data_w_rt_range[exp_row, "rt"]
        
        
        
        
        #add the metabolite and its meta data to the dataframe of found matches 
        mat_met_exp_ids_that_matched_a_lib_met_rt_col_index = which(names(matched_metabolites) == "exp_ids_that_matched_a_lib_met_rt")
        
        mat_met_rt_col_index = which(names(matched_metabolites) == "rt")
        
        mat_met_matches_found_col_index = which(names(matched_metabolites) == "matching_metabolite")
        mat_met_db_id_col_index = which(names(matched_metabolites) == "CAS_Number")
        mat_met_Formula_col_index = which(names(matched_metabolites) == "Formula")
        mat_met_Isomeric_SMILES_formula_col_index = which(names(matched_metabolites) == "Isomeric_SMILES")
        mat_met_InChI_formula_col_index = which(names(matched_metabolites) == "InChI")
        
        
        # Construct a new row (using the column names from the library)
        new_row <- data.frame(
          exp_ids_that_matched_a_lib_met_rt = I(list(exp_id)),
          rt = exp_rt,
          matching_metabolite = match_rt_lib[lib_row, "metabolite_identification"],
          CAS_Number = match_rt_lib[lib_row, "CAS_Number"],
          Formula = match_rt_lib[lib_row, "Formula"],
          Isomeric_SMILES = match_rt_lib[lib_row, "Isomeric_SMILES"],
          InChI = match_rt_lib[lib_row, "InChI"],
          stringsAsFactors = FALSE
        )
        
        # Append the new row to the matched_metabolites data frame
        matched_metabolites <- rbind(matched_metabolites, new_row)
        
        
      }#close of if statement when match is found 
      
      
    }#end loop of experimental rows 
  }
}#end of the loop over the lib rows 




#write out the identified metabolites in the matrix data 
write.csv(matched_metabolites, "E:/748/Assessment 3/matched_metabolites.csv")