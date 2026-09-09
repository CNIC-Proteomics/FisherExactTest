# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Author: Samuel Lozano Juárez 
# Date: 03/08/2026
# Institution: Spanish National Center for Cardiovascular Research (CNIC)
# Contact: samuel.lozano@cnic.es
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# IMPORTATIONS AND SETUP ----
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

library(future.apply)
library(dplyr)
library(tidyr)
library(data.table)
library(stringr)
library(DiscreteFDR)
library(optparse)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# PARAMETERS READING ----
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Define arguments
option_list <- list(
  make_option(
    c("-i", "--input"),
    type = "character",
    help = "Tabla TSV con la salida de iSanXoT o DIANN",
    metavar = "FILE"
  ),
  make_option(
    c("-g", "--groups"),
    type = "character",
    help = "Tabla TSV con los grupos y las columnas de cada grupo",
    metavar = "FILE"
  ),
  make_option(
    c("-c", "--col_name"),
    type = "character",
    help = "Nombre de la columna a utilizar",
    metavar = "STRING"
  )
)

# Parsear argumentos
opt <- parse_args(OptionParser(option_list = option_list))

# Comprobar que se han proporcionado ambos
if (is.null(opt$input)) {
  stop("Error: debes proporcionar un archivo con -i")
}

if (is.null(opt$groups)) {
  stop("Error: debes proporcionar un archivo con -g")
}

# Leer tablas
input_table <- fread(
  opt$input,
  sep = "\t",
  check.names = F
)

groups_table <- fread(
  opt$groups,
  sep = "\t",
  check.names = F
)

#define the groups
grupos <- lapply(groups_table, function(x) {
  x[!is.na(x) & x != ""]
})

names(grupos) <- colnames(groups_table)

if (length(grupos)<2) {
  stop("Error: el número de grupos es inferior a 2. Por favor chequea la groups_table")
}

#transform the dataframe to modify it by reference (so it won't be loaded in memory)
setDT(input_table)

#define the num of columns
sample_cols <- unlist(grupos, use.names = FALSE)
N <- length(sample_cols)

#create the presence matrix
pres_mat <- !is.na(as.matrix(input_table[, ..sample_cols]))
storage.mode(pres_mat) <- "integer" #binarize (0-1) the matrix

#get the total number of presences by peptide
K <- rowSums(pres_mat)
ids <- input_table[[opt$col_name]]

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# USEFUL FUNCTIONS ----
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# índices de columnas de cada grupo dentro de pres_mat
group_idx <- lapply(grupos, function(x)
  match(x, sample_cols))

# function to compute probability table
table_prob <- function(tab) {
  rs <- rowSums(tab)
  cs <- colSums(tab)
  N  <- sum(tab)
  exp(
    sum(lfactorial(rs)) + sum(lfactorial(cs)) -
      lfactorial(N) - sum(lfactorial(tab))
  )
}

#create the discrete supports (for Discrete BH adjust)
build_supports <- function(n_g, N) {
  Ks <- 0:N #crate all the possible K values (remember, K is the number of presences in the population)
  greater <- vector("list", length(Ks))
  less <- vector("list", length(Ks))
  for (k in Ks) { #iter along all the possible values of K (presence in the whole population)
    m <- k; n <- N - k
    xmin <- max(0, n_g - n) #the min number of presences in the group is 0 or the number of the elements in the group - the total absences in population
    xmax <- min(n_g, m) #the max number of presences is the minimum between the total of elements in the group and the total of presences in population
    xs <- xmin:xmax
    greater[[k + 1]] <- sort(unique(phyper(xs - 1, m, n, n_g, lower.tail = FALSE))) #possible p-values given the xmin-xmax range for enrichment 
    less[[k + 1]] <- sort(unique(phyper(xs,     m, n, n_g, lower.tail = TRUE))) #possible p-values given the xmin-xmax range for depletion 
  }
  list(greater = greater, less = less)
}

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# ONLY TWO GROUPS ----
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# if there're only two groups, onlye the hypergeom test is needed (not the multicategory)
if (length(grupos)==2){
  #function to analyze a whole group
  analizar_grupo <- function(g_name) {
    cols_g <- grupos[[g_name]]
    n_g <- length(cols_g)
    
    x <- rowSums(pres_mat[, cols_g, drop = FALSE])    # presences in the group
    
    # raw p-valores, vectorizados (equivalent to fisher.test 2x2 one-sided)
    p_greater <- phyper(x - 1, K, N - K, n_g, lower.tail = FALSE)
    
    #build the support vectors
    sup <- build_supports(n_g, N) #recover the possible theoretical pvalues
    pCDF_greater <- sup$greater[K + 1]
    
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    # FDR 0.05
    #adjust the FDR (apply an alpha of 0.05)
    dbh_greater_005 <- discrete.BH(p_greater, pCDF_greater, direction = "su", alpha = 0.05)
    #and create the pass/not pass vector
    pass_DBH_005_enrichment <- rep(FALSE, length(ids))
    pass_DBH_005_enrichment[dbh_greater_005$Indices] <- TRUE
    
    data.table(
      ID = ids,
      grupo = g_name,
      P.Val.Enrich = p_greater,
      Signif_Enrich = pass_DBH_005_enrichment
    )
  }
  
  #merge all the groups analysis
  plan(multisession, workers = min(2, future::availableCores()))
  resultados_lista <- future_lapply(names(grupos), analizar_grupo, future.seed = TRUE)
  plan(sequential)
  
  resultados_1vsall <- rbindlist(resultados_lista)
  one_vs_all_wide <- dcast(
    resultados_1vsall,
    ID ~ grupo,
    value.var = setdiff(names(resultados_1vsall), c("ID", "grupo"))
  )
  combined <- input_table %>% left_join(one_vs_all_wide, by = setNames("ID", opt$col_name))
  path2save <- gsub("\\.tsv","_FET.tsv",opt$input)
  fwrite(combined, path2save, sep = "\t")
}



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# MORE THAN TWO GROUPS ----
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

if (length(grupos)>2){
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  # FISHER EXACT TEST ONE GROUP VS THE REST
  
  #analyze the whole group
  analizar_grupo <- function(g_name) {
    cols_g <- grupos[[g_name]]
    n_g <- length(cols_g)
    
    x <- rowSums(pres_mat[, cols_g, drop = FALSE])    # presences in the group
    
    # raw p-valores, vectorizados (equivalent to fisher.test 2x2 one-sided)
    p_greater <- phyper(x - 1, K, N - K, n_g, lower.tail = FALSE)
    p_less <- phyper(x, K, N - K, n_g, lower.tail = TRUE)
    
    #build the support vectors
    sup <- build_supports(n_g, N) #recover the possible theoretical pvalues
    pCDF_greater <- sup$greater[K + 1]
    pCDF_less    <- sup$less[K + 1]
    
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    # FDR 0.05
    #adjust the FDR (apply an alpha of 0.05)
    dbh_greater_005 <- discrete.BH(p_greater, pCDF_greater, direction = "su", alpha = 0.05)
    dbh_less_005 <- discrete.BH(p_less, pCDF_less, direction = "su", alpha = 0.05)
    #and create the pass/not pass vector
    pass_DBH_005_enrichment <- rep(FALSE, length(ids))
    pass_DBH_005_enrichment[dbh_greater_005$Indices] <- TRUE
    
    pass_DBH_005_depletion <- rep(FALSE, length(ids))
    pass_DBH_005_depletion[dbh_less_005$Indices] <- TRUE
    
    data.table(
      ID = ids,
      grupo = g_name,
      P.Val.Enrich = p_greater,
      P.Val.Deple = p_less,
      Signif_Enrich = pass_DBH_005_enrichment,
      Signif_Deple = pass_DBH_005_depletion
    )
  }
  
  #merge all the groups analysis
  plan(multisession, workers = min(4, future::availableCores()))
  resultados_lista <- future_lapply(names(grupos), analizar_grupo, future.seed = TRUE)
  plan(sequential)
  resultados_1vsall <- rbindlist(resultados_lista)
  
  one_vs_all_wide <- dcast(
    resultados_1vsall,
    ID ~ grupo,
    value.var = setdiff(names(resultados_1vsall), c("ID", "grupo"))
  )
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  # FISHER EXACT TEST MULTICATEGORY (DISTRIBUTION VS NULL HYPOTHESIS)
  
  group_sizes <- sapply(grupos, length)
  #function for fisher exact test in one row
  analizar_pgm <- function(i){
    pres <- pres_mat[i, ]
    presentes <- sapply(group_idx, function(idx) sum(pres[idx]))
    ausentes  <- group_sizes - presentes
    tabla <- rbind(
      Present = presentes,
      Absent  = ausentes
    )
    ## Fisher exacto 2x4
    p <- fisher.test(tabla, simulate.p.value = F)$p.value
    ## probabilidad de la tabla observada bajo margenes fijos
    p_obs <- table_prob(tabla)
    ## mid-p value
    p_mid <- p - p_obs/2
    list(
      ID = ids[i],
      p = p,
      pmid = p_mid
    )
  }
  
  plan(multisession,workers = min(future::availableCores(), 6))
  raw_pvalues <- future_lapply(seq_len(nrow(pres_mat)),analizar_pgm,future.seed = TRUE)
  plan(sequential)
  
  pvalues <- vapply(raw_pvalues, `[[`, numeric(1), "p")
  pmid <- vapply(raw_pvalues, `[[`, numeric(1), "pmid")
  adj_pmid <- p.adjust(pmid, method = "BH")
  
  resultados_multicat <- data.table(
    ID = ids,
    P.Val = pvalues,
    P.Mid = pmid,
    Adj.P.Val = adj_pmid
  )
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  # COMBINE THE TWO TESTS
  combined <- input_table %>% left_join(one_vs_all_wide, by = setNames("ID", opt$col_name)) %>% left_join(resultados_multicat, by = setNames("ID", opt$col_name))
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  # SUBSAMPLE TO KEEP ONLY THE SIGNIFICATIVE ROWS
  significatives <- combined %>% filter(if_any(contains("Signif"), ~ .x == TRUE) |Adj.P.Val < 0.05)
  
  #transform the dataframe to modify it by reference (so it won't be loaded in memory)
  setDT(significatives)
  
  #create the presence matrix
  pres_mat <- !is.na(as.matrix(significatives[, ..sample_cols]))
  storage.mode(pres_mat) <- "integer" #binarize (0-1) the matrix
  
  #get the total number of presences by peptide
  K <- rowSums(pres_mat)
  sig_ids <- significatives[[opt$col_name]]
  
  #create the pairwise comparisons
  comparisons <- combn(names(grupos), 2, simplify = FALSE)
  
  #create the function to analyse one comparison
  analizar_comparacion <- function(comp){
    
    g1 <- comp[1] #select the first group of the comparison
    g2 <- comp[2] #select the second group of the comparison
    
    cols_g1 <- grupos[[g1]]
    cols_g2 <- grupos[[g2]]
    
    n1 <- length(cols_g1)
    n2 <- length(cols_g2)
    
    # Total population = only the two compared groups
    N <- n1 + n2
    
    # Presences in group 1
    x <- rowSums(pres_mat[, cols_g1, drop = FALSE])
    
    # Presences in the whole population (group1 + group2)
    K <- x + rowSums(pres_mat[, cols_g2, drop = FALSE])
    
    #----------------------------------------------------------
    # Raw p-values
    #----------------------------------------------------------
    
    p_greater <- phyper(x - 1, K, N - K, n1, lower.tail = FALSE)
    p_less <- phyper(x, K, N - K, n1, lower.tail = TRUE)
    
    #----------------------------------------------------------
    # Theoretical p-value supports
    #----------------------------------------------------------
    
    sup <- build_supports(n1, N)
    
    pCDF_greater <- sup$greater[K + 1]
    pCDF_less    <- sup$less[K + 1]
    
    #----------------------------------------------------------
    # FDR 0.05
    #----------------------------------------------------------
    
    dbh_greater_005 <- discrete.BH(
      p_greater,
      pCDF_greater,
      direction = "su",
      alpha = 0.05
    )
    
    dbh_less_005 <- discrete.BH(
      p_less,
      pCDF_less,
      direction = "su",
      alpha = 0.05
    )
    
    pass_DBH_005_enrichment <- rep(FALSE, length(sig_ids))
    pass_DBH_005_enrichment[dbh_greater_005$Indices] <- TRUE
    
    pass_DBH_005_depletion <- rep(FALSE, length(sig_ids))
    pass_DBH_005_depletion[dbh_less_005$Indices] <- TRUE
    
    #----------------------------------------------------------
    # Output
    #----------------------------------------------------------
    
    data.table(
      ID = sig_ids,
      comparison = paste0(g1, "_", g2),
      
      P.Val.Enrich = p_greater,
      P.Val.Deple = p_less,
      
      Signif_Enrich = pass_DBH_005_enrichment,
      Signif_Deple = pass_DBH_005_depletion
    )
  }
  
  plan(multisession, workers = min(length(comparisons),future::availableCores()))
  resultados_lista <- future_lapply(
    comparisons,
    analizar_comparacion,
    future.seed = TRUE
  )
  
  plan(sequential)
  resultados_1vs1 <- rbindlist(resultados_lista)
  
  one_vs_one_wide <- dcast(
    resultados_1vs1,
    ID ~ comparison,
    value.var = setdiff(names(resultados_1vs1), c("ID", "comparison"))
  )
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  # MERGE ALL THE RESULTS
  
  combined_all <- combined %>% left_join(one_vs_one_wide, by = setNames("ID", opt$col_name))
  
  path2save <- gsub("\\.tsv","_FET.tsv",opt$input)
  fwrite(combined_all, path2save, sep = "\t")
}
