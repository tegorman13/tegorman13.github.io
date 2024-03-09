#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
source(here::here('Functions','fun_cv.R'))
# Define the elements of each column as vectors
programming <- c("R - Rmarkdown - Quarto", "JavaScript","Bash" ,"Python", "MATLAB")
computational <- c("Artificial Neural Networks", "Exemplar Models", "Bayesian Statistics", "Mixed Effect Models", "Dimensionality Reduction")
experimental <- c("Behavioral Tasks","Online Data Collection - JsPsych" ,"Survey Data - Qualtrics" ,"Mechanical Turk" , "MRI", "EEG")

pad_vectors <- function(...){
  args <- list(...)
  max_length <- max(sapply(args, length))
  lapply(args, function(vec) {
    vec[length(vec):max_length] <- ""
    return(vec)
  })
}

padded_vectors <- pad_vectors(programming, computational, experimental)
data <- tibble(Programming = padded_vectors[[1]],
               Computational = padded_vectors[[2]],
               Experimental = padded_vectors[[3]])

pander::pander(data, style = 'rmarkdown', split.table=Inf, justify = 'left')

#fontawesome::fa(name = "r-project", fill = "steelblue")
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
