library(readr)
library(dplyr)
library(tidyr)
library(knitr)
library(stringr)
library(purrr)

## Extract error counts from metaquast vs. valet comparison output
get_comp_summary_vals <- function(comp_summary_file, comp_suspicious_file){
    if(file.size(comp_summary_file) == 0){
        mq_errors = 0; v_valid = 0; valid_overlap = 0
    }else{
        comparison <- read_lines(comp_summary_file)
        # Metaquast errors
        metaquast <- comparison[3] %>% str_split(pattern = "\t")  
        ## mq - errors
        mq_errors <- metaquast[[1]][3] %>% as.numeric()
        
        # Valet
        valet <- comparison[7] %>% str_split(pattern = "\t") 
        ## valet - valid
        v_flagged<- valet[[1]][3] %>% as.numeric()
        
        # Valid overlap
        flagged_overlap <- metaquast[[1]][2] %>% as.numeric()
    }   
    
    if(file.size(comp_suspicious_file) == 0){
        v_flagged = 0; flagged_overlap = 0
    }else{
        comparison <- read_lines(comp_suspicious_file)
        # Metaquast
        metaquast <- comparison[3] %>% str_split(pattern = "\t")  
        # Valet
        valet <- comparison[7] %>% str_split(pattern = "\t") 
        ## valet - valid
        v_valid <- valet[[1]][3] %>% as.numeric()
        # Valid overlap
        valid_overlap <- metaquast[[1]][2] %>% as.numeric()
    }
    data_frame(mq_errors, v_flagged, flagged_overlap, v_valid, valid_overlap)
}


## extract table summary values from metaquast output report.tsv file
get_asm_val <- function(f_line){
    f_line %>% str_split("  ") %>% flatten() %>% .[2] %>% as.numeric()
}
get_assembly_summary_vals <- function(assembly_file){
    assembly_metrics <- read_tsv(assembly_file, col_names = FALSE)
    assembly_size <- assembly_metrics %>% 
        filter(X1 == "Total length (>= 5000 bp)") %>% .$X2 %>% 
            as.numeric() %>% prettyNum(big.mark=",")
    n_ctgs <-  assembly_metrics %>% 
        filter(X1 == "# contigs (>= 5000 bp)") %>% .$X2 %>% 
            as.numeric() %>% prettyNum(big.mark=",")
    n50 <- assembly_metrics %>% 
        filter(X1 == "N50") %>% .$X2 %>% as.numeric() %>% prettyNum(big.mark=",")
    na50 <- assembly_metrics %>% 
        filter(X1 == "NA50") %>% .$X2 %>% as.numeric() %>% prettyNum(big.mark=",")
    mpk <- assembly_metrics %>% 
        filter(X1 == "# mismatches per 100 kbp") %>% .$X2 %>% as.numeric()
    data_frame(assembly_size, n_ctgs, n50, na50, mpk)
}

## extract assembly metrics and error comparison values from output files
get_comp_assembly_tbl_values <- function(assembler, results_dir){
    assembly_file <- paste0(results_dir,"/", assembler, "/report.tsv")
    comp_file <- paste0(results_dir,"/COMP_QUAST-w-reads-5k-", 
                        assembler, 
                        "/comparison_reference.summary.results")
    
    comp_summary_file <- paste0(results_dir,"/COMP_QUAST-w-reads-5k-", 
                                assembler, 
                                "/comparison_reference.summary.results")
    comp_suspicious_file <- paste0(results_dir,"/COMP_QUAST-w-reads-5k-", 
                                   assembler, 
                                   "/comparison_reference.suspicious.results")
    data_frame(assembler) %>% 
        bind_cols(get_assembly_summary_vals(assembly_file)) %>% 
        bind_cols(get_comp_summary_vals(comp_summary_file, 
                                        comp_suspicious_file))
}

## reformat table for supplemental material
format_tbl <- function(comp_tbl){
    comp_tbl %>%
        select(assembler, assembly_size, n_ctgs, n50, na50, mpk, mq_errors, 
               v_flagged, flagged_overlap, v_valid, valid_overlap) %>% 
        rename(Assembler = assembler, `Len` = assembly_size, Ctgs = n_ctgs,
               `N50` = n50, `NA50` = na50,
               `MPK` = mpk,
               `MQ Errs` = mq_errors,
               `Flg` = v_flagged, `Flg Overlap` = flagged_overlap,
               `Vld` = v_valid, `Vld Overlap` = valid_overlap) %>% 
        gather(Metric, value,-Assembler) %>% spread(Assembler,value)
}

# function to batch process files and return single dataframe
make_comparison_tbl <- function(results_dir, assemblers){
    assemblers %>% map(get_comp_assembly_tbl_values, results_dir) %>% 
        bind_rows() %>% format_tbl()
}