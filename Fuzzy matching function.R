




fuzzy.match = function(data, id.col, parameter.df, n.results){
  # Author: Daniel Gardiner
  # email: daniel.gardiner@phe.gov.uk
  ###############################################################################
  # Functions needed for analysis
  
  # Perform Fuzzy matching                                                      # 
  
  fuzz.m = function(vec1, vec2 = NULL, FUN, thresh, na.value = 0){
    # function to create a data.frame of raw and adjusted jarowinkler or
    # levenshtein edit distance scores (dependent on FUN argument)
    # NOTE: all missing entries (NAs) will be assigned a score of 0
    #
    # step 1: an an upper triangular matrix is created with entry i,j the 
    #         raw score of the ith entry of vec1 compared to the jth entry of vec2
    # step 2: this matrix is 'melted' to give a data.frame of all raw comparisons 
    # step 3: an adjusted score column is added 
    #
    # Args: vec1:    a vector
    #       vec2:    a vector
    #       FUN:    either jarowinkler or levenshteinSim
    #       thresh: threshold value
    #
    # output: var1: comparitor 1
    #         var2: comparitor 2
    #         raw.score: raw jarowinkler OR levenshtein edit distance score
    #         adj.score: adjusted score (according to thresh argument)
    
    #load packages
    library(RecordLinkage)
    
    library(reshape2)
    
    # if two vectors are supplied determine shorter and longer vector
    if(is.null(vec2)){
      NULL
    } else if (length(vec1) < length(vec2)) {
      short.vec = vec1
      long.vec = vec2
    } else if (length(vec1) > length(vec2)){
      short.vec = vec2
      long.vec = vec1
    } else {
      short.vec = vec1 # note: else vectors are actually of equal length, 
      long.vec = vec2  #       short/long is conventiant naming convention  
    }
    
    # create a square matrix consisting of vec1^2 or vec1 * vec2 NAs   
    if(is.null(vec2)){
      mat = rep(NA, length(vec1)^2)
      mat = matrix(mat, nrow = length(vec1), ncol = length(vec1),
                   dimnames = list(vec1, vec1))  
    } else {  
      mat = rep(NA, length(short.vec) * length(long.vec))
      mat = matrix(rep(NA, length(short.vec) * length(long.vec)), 
                   nrow = length(short.vec), ncol = length(long.vec))
    }
    
    # take the ith entry of vec and generate a jaro winkler score compared to each other
    # entry of vec, store as a row in a matrix, loop over all entries in vec to fill 
    # all rows of the matrix
    if(is.null(vec2)){
      for(i in 1:length(vec1)){
        mat[i, ] = eval(call(FUN, vec1[i], vec1))  
      }
    } else {
      for(i in 1:length(short.vec)){
        mat[i, ] = eval(call(FUN, short.vec[i], long.vec))
      }    
    }
    
    # replace all NA values with 0's as raw score
    mat[is.na(mat)] = 0
    
    # convert matrix to upper triangular form
    mat = ifelse(upper.tri(mat), mat, "REMOVE")
    
    
    # assign row and column names
    if(is.null(vec2)){
      mat = matrix(mat, nrow = length(vec1), ncol = length(vec1),
                   dimnames = list(vec1, vec1))
    } else {
      mat = matrix(mat, nrow = length(short.vec), ncol = length(long.vec),
                   dimnames = list(short.vec, long.vec))
    }
    
    # melt data to give each row as a comparison
    mat.m = melt(mat, id.vars = row.names(mat), measure.vars = colnames(mat))
    
    # keep only those elements on the upper triangular matrix
    mat.m = mat.m[mat.m[, 3] != "REMOVE", ]
    
    # make raw score numeric
    mat.m[, 3] = as.numeric(as.character(mat.m[, 3]))
    
    # add column of adjusted scores
    mat.m[, 4] = ifelse(mat.m[, 3] >= thresh,
                        1 - ((1 - mat.m[, 3])/(1 - thresh)), 0)
    
    # if Var1 or Var2 is NA assign value na.value 
    mat.m[is.na(mat.m$Var1) | is.na(mat.m$Var2), 4] = na.value
    
    # rename columns
    colnames(mat.m) = c("var1", "var2", "raw.score", "adj.score")
    
    # output data.frame from function
    # note: each row represents a single comparison between entries
    mat.m  
  }
  
  ####
  # derive results using the matching function
  
  # format parameters
  
  char.cols = c("match.cols", "func")
  
  parameter.df[, char.cols] = sapply(parameter.df[, char.cols], as.character)
    
  num.cols = c("thresholds", "NA.value", "weightings")
  
  parameter.df[, num.cols] = sapply(parameter.df[, num.cols], as.numeric)
  
  
  # format data
  
  data[, id.col] = as.character(data[, id.col])
  
  data[, parameter.df$match.cols] = sapply(data[, parameter.df$match.cols], as.character)
  
  
  # generate results by applying function to each field
  
  results = NULL
  
  for(i in 1:nrow(parameter.df)){
    
    temp = fuzz.m(vec1 = data[, parameter.df[i, "match.cols"]],
                  FUN =  parameter.df[i, "func"], 
                  thresh = parameter.df[i, "thresholds"],
                  na.value = parameter.df[i, "NA.value"])
    
    temp = parameter.df[i, "weightings"] * temp$adj.score
    
    results = cbind(results, temp)
  }
  
  # sum accross fields to give overall scores for each comparison
  
  overall.score = round(apply(results, 1, sum), 2)
  
  # append overall scores onto each comparison id 
  

  final.results = fuzz.m(vec1 = data[, id.col], 
                         FUN =  "levenshteinSim", 
                         thresh = 1,
                         na.value = 1)[, 1:2]
  
  final.results = data.frame(final.results, overall.score)
  
  # order results by overall score
  
  final.results = final.results[rev(order(final.results$overall.score)), ]
  
  # convert ids to characters
  
  final.results$var1 = as.character(final.results$var1)
  
  final.results$var2 = as.character(final.results$var2)
  
  # generate matches
  
  
  n.results = min(n.results, nrow(final.results))
      
  list.results = vector("list", n.results)
  
  for(i in 1:n.results){
    list.results[[i]] = data.frame(data[data[, id.col] %in% final.results[i, c("var1", "var2")], ],
                                   overall.score = final.results[i, c("overall.score")])
  }
  
  list.results
}



########
# EXAMPLE 

# generate example data

data = data.frame(ID = c(1,  2, 3, 4, 5, 6, 7),
                  Firstname = c("Steve", "Stve",  "Clint", "Julianne",  "Meryl", "Mark", "Maryl"),
                  Surname = c("Carell",	"Corel", NA, "Moore", "Streep", "Ruffalo", "Streep"))


# generate parameter data frame

parameters = data.frame(match.cols = c("Firstname", "Surname"),
                        thresholds = c(0.3, 0.3),
                        NA.value = c(0.1, 0.1),
                        weightings = c(0.6, 0.8),
                        func = "levenshteinSim")

# use function 

fuzzy.match(data, id.col = "ID", parameter.df = parameters, n.results = 5)

