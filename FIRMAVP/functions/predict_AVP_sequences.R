require(e1071)
require(caret)
source("./functions/feature_extractors.R")
#x<-"D:\\PARGT_Windows\\selected_train_test_merged_file.csv"
#y<-"D:\\PARGT_Windows\\input_seq.csv"


predict_results<- function(x,y) {
  library(e1071)
  library(caret)

  #file_n<-x
  
  #file_other<-y
  
  data1 <- x#read.csv(file_n, header = TRUE)
  
  
  cols<-ncol(data1)
  
  data_backup<-data1
  
  nrows_training<-nrow(data1)
  
  data_other <- y#read.csv(file_other, header = TRUE)
  data_other_backup<-data_other
  data_other2<-data_other
  cols2<-ncol(data_other)
  
  DF2<-data_other
  data_other_backup<-data_other
  nrows_testing<-nrow(data_other)
  
  
  data1<-data1[-c(cols)]
  DF1<-data1

  data_norm<-preProcess(data1,method=c("center", "scale"))
  data1<-predict(data_norm, data1)
  
  data1["Output"]<-data_backup[,cols]
  data_other<-predict(data_norm, data_other)
  train<-data1
  ncol(train)
  test<-data_other
  
  kayes<-123
  set.seed(kayes)
  train$Output<-as.factor(train$Output)

  tmodel2<-tune(svm, Output~., data = train, ranges = list(epsilon =0, cost=8),
                probability = TRUE)
  #tmodel2$best.parameters

  mymodel2<-tmodel2$best.model
  results<-predict(mymodel2, test, probability=TRUE)

  res_vec<-c()
  for(i in 1:nrows_testing){
    if(results[i]==1){
      res_vec<-c(res_vec,1)
    }
    else{
      res_vec<-c(res_vec,-1)
    }
  }
  return(list(res_vec, results))
  
}

predict_AVP_sequences <- function(training_file,
                                  testing_file,
                                  fasta_file_path = "./FIRMAVP/data/input_seq.fasta",
                                  raw_sequence = NA,
                                  generate_features = FALSE,
                                  is_html = TRUE){ 
  predictions <- predict_results(training_file,testing_file)
  if(!is.na(raw_sequence)){
    lines = list(raw_sequence)
    sequences_name_input <- " "
  } else{
    lines = readLines(fasta_file_path)
    sequences_name_input= c()
  }
  sequences_input= c()
  for(line in lines){
    if(startsWith(line, ">")){
      sequences_name_input <- append(sequences_name_input, gsub(">", "",line))
    }
    else if(line != ""){
      sequences_input <- append(sequences_input, line)
    }
  }
  # Total number of predicted AVP sequences = 2 
  # AAQ17160.1 vancomycin/teicoplanin A-type resistance protein VanA (plasmid) [Staphylococcus aureus] 
  # CAA94438.1 D-alanine:D-alanine ligase-related protein, partial [Streptococcus equinus]
  
  # count_resistence_sequences <- which(predictions == 1)
  # one_line=paste("Total number of predicted AVP sequences = ",length(count_resistence_sequences)) 
  # sequences_to_write <- c(one_line, sequences_name_input[count_resistence_sequences])
  # if(is_html){
  #   out <- paste(sequences_to_write, collapse = " <br/> ")
  # } else{
  #   out <- paste(sequences_to_write, collapse = "/n")
  # }
  result = data.frame(attributes(predictions[[2]])$probabilities)
  result <- cbind(result, sequences_input)
  result <- cbind(result, sequences_name_input)
  colnames(result) <- c("AVP", "Non-AVP", 'Sequence', 'Peptide')
  return(result)
}



add_new_sequences <- function(input_fasta_filepath = "./FIRMAVP/data/input_seq.fasta",
                              outcome_value = 1){
  new_data <- feature_extraction(input_fasta_filepath)
  new_data$Output <- outcome_value
  return(new_data)
}




