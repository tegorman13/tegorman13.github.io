pacman::p_load(tidyverse,dplyr,here,purrr,tidyr,glue,yaml,vitae, knitr, kableExtra)

read_yaml_as_tbl <- function(filepath) {
  read_yaml(filepath) %>%
    map_dfr(function(x) {
      sapply(names(x), function(y) {
        # wrap elements in list if length > 1
        if (length(x[[y]]) > 1) {
          x[[y]] <<- list(x[[y]])
        }
      })
      return(as_tibble_row(x))
    }) %>%
    return()
} 
replace_amps <- function(x) {
  return(str_replace_all(x, "&", "\\\\&"))
}
format_md2tex <- function(x) {
  return(x %>% 
           str_replace_all("_(.*)_", "\\\\textit{\\1}") %>%  # italics
           str_replace_all("\\*\\*(.*)\\*\\*", "\\\\textbf{\\1}") %>% # TODO bold
           str_replace_all('\\+', '\\\\plus') %>% 
           replace_amps()
         )
}

split_md_list <- function(x) {
  return(x %>% 
           str_remove('^- ') %>% 
           str_remove("\n$") %>% 
           str_split('\n- ')
         )
}