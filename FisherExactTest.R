# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Author: Samuel Lozano Juárez 
# Date: 24/09/2026
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
library(DiscreteTests)
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
  ),
  make_option(
    c("-f", "--fdr"),
    type = "numeric",
    default = 0.05,
    help = "FDR threshold (default: 0.05)",
    metavar = "NUMBER"
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

fdr_threshold <- opt$fdr

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

#create the discrete supports for more than 2 groups (for Discrete BH adjust)
build_supports_fisher_2xK <- function(n_g, tolerance = 1e-12, verbose = TRUE) {
  
  # ------------------------------------------------------------
  # n_g: vector with the sample size of each group
  #
  # Returns a list where: supports[[M + 1]]
  #
  # contains all possible two-sided Fisher p-values for a 2 x K table with:
  #   - group sizes = n_g
  #   - total number of presences = M
  #
  # The two-sided p-value is defined as: sum P(T) for all tables T whose probability is <= P(observed table)
  #
  # which is the definition used by fisher.test() for two-sided exact tests.
  # ------------------------------------------------------------
  
  
  # Basic checks
  if (length(n_g) < 2) {
    stop("n_g must contain at least two groups.")
  }
  
  if (any(n_g < 0) || any(n_g != floor(n_g))) {
    stop("n_g must contain non-negative integers.")
  }
  
  N <- sum(n_g)
  K <- length(n_g)
  
  if (N == 0) {
    stop("The total sample size must be > 0.")
  }
  
  
  # ------------------------------------------------------------
  # Generate all possible vectors x = (x1,...,xK) such that:
  #
  #   sum(x) = M
  #   0 <= xj <= nj
  #
  # We do this recursively rather than constructing a huge Cartesian product.
  # ------------------------------------------------------------
  
  generate_tables <- function(n_g, M) {
    
    # Number of groups
    K <- length(n_g)
    
    # Vector that will temporarily store the number of presences assigned to each group.
    #
    # Initially: x = (0, 0, ..., 0)
    x <- integer(K)
    
    # List where we will store every valid combination of presences.
    tables <- list()
    
    # Counter used to indicate the position of the next valid combination in 'tables'.
    counter <- 0L
    
    recurse <- function(j, remaining) {
      # this function decides how many presences we assign to the group j, given the "remaining" presences that are still remaining to be assigned
      if (j == K) { # Last group is completely determined
        value <- remaining
        
        if (value >= 0 && value <= n_g[j]) { #we check that the number of presences is not 0, and is not grater than the number of samples in the group
          
          x[j] <<- value
          counter <<- counter + 1L
          tables[[counter]] <<- x
        }
        
        return(invisible(NULL))
      }
      
      # Minimum and maximum value that x[j] can take IF WE ARE NOT IN THE LAST GROUP

      min_x <- max( #the number of presences in the group j cannot be less than 0, and must be enough so the remaining ones fit in the other groups (given the number of samples of the rest of the groups)
        0L,
        remaining - sum(n_g[(j + 1):K])
      )
      
      max_x <- min( #the maximum number of presences is restricted by the remaining presences to be assigned, and the number of samples in this group
        n_g[j],
        remaining
      )
      
      if (min_x > max_x) { #if min is greater than max, then no solution is feasible for this branch 
        return(invisible(NULL))
      }
      
      for (value in min_x:max_x) { #in this case we already know that min_x <= max_x, we try then all the possible values from min to max
        
        x[j] <<- value
        
        recurse( #for each possible value, we distribute the remaining values among the next group, using recursivity
          j + 1L,
          remaining - value
        )
      }
      
      invisible(NULL)
    }
    
    recurse(1L, M) #this is the initial call
    
    do.call(rbind, tables)
  }
  
  
  # ------------------------------------------------------------
  # Precompute log-combinations:
  #
  # log(C(n, x))
  #
  # This avoids repeatedly calculating choose(), and is numerically more stable for moderately large n.
  # ------------------------------------------------------------
  
  log_choose <- lapply(
    n_g,
    function(n) {
      lchoose(n, 0:n)
    }
  )
  
  
  # ------------------------------------------------------------
  # Calculate the support for one value of M
  # ------------------------------------------------------------
  
  build_one_support <- function(M) {
    
    # All possible tables for this M (M=total number of presences, n_g=size of each group)
    x_mat <- generate_tables(n_g, M)
    
    # log probability of each table:
    #
    # log P(T) = sum_j log(C(n_j, x_j)) - log(C(N, M))
    
    log_prob <- numeric(nrow(x_mat))
    
    #we get the probability of each table
    for (j in seq_len(K)) { #iter over each group
      log_prob <- log_prob + log_choose[[j]][x_mat[, j] + 1L]
    }
    
    log_prob <- log_prob - lchoose(N, M)
    
    # Convert to probabilities
    prob <- exp(log_prob)
    
    # ----------------------------------------------------------
    # Every possible table has a Fisher two-sided p-value:
    #
    # P-value(T) = sum_{T': P(T') <= P(T)} P(T')
    #
    # ----------------------------------------------------------
    
    # Sort probabilities
    ord <- order(prob)
    prob_sorted <- prob[ord]
    
    # Cumulative probability
    cum_prob <- cumsum(prob_sorted)
    
    # Because several tables can have exactly the same probability, all tables with equal probability must receive the same p-value.
    #
    # Numerical tolerance is used because probabilities are floating-point numbers.
    
    #create the vector where we'll store the p-values
    p_sorted <- numeric(length(prob_sorted))
    start <- 1L
    
    
    while (start <= length(prob_sorted)) { #iter over all the possible probabilities
      
      current <- prob_sorted[start] #get the actual probability
      end <- start
      
      #check if there is any probability that has the same value as the actual probability
      while (end < length(prob_sorted) && abs(prob_sorted[end + 1L] - current) <= tolerance * max(1, abs(current))) {
        end <- end + 1L
      }
      
      # The p-value includes all tables with probability <= current probability.
      p_value <- cum_prob[end]
      p_value <- min(p_value, 1) #adjust to 1 as the maximum pvalue, to avoid possible cummulative errors
      
      p_sorted[start:end] <- p_value
      
      start <- end + 1L
    }
    
    # Return the unique possible p-values
    sort(unique(p_sorted))
  }
  
  
  # ------------------------------------------------------------
  # Build supports for every possible number of presences
  # ------------------------------------------------------------
  
  supports <- vector("list", N + 1L)
  
  names(supports) <- paste0("M_", 0:N)
  
  if (verbose) {
    message(
      "Building Fisher 2xK supports for ",
      K, " groups and N = ", N, " samples."
    )
  }
  
  for (M in 0:N) {
    
    if (verbose && (M %% max(1, floor(N / 10)) == 0)) {
      message("  M = ", M, "/", N)
    }
    
    supports[[M + 1L]] <- build_one_support(M)
  }
  
  return(supports)
}

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# ONLY TWO GROUPS ----
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# if there're only two groups, only the hypergeom test is needed
if (length(grupos)==2){
  g1 <- grupos[[1]]
  g2 <- grupos[[2]]
  
  n1 <- length(g1)
  n2 <- length(g2)
  
  # Number of presences in each group
  x1 <- rowSums(pres_mat[, g1, drop = FALSE])
  x2 <- rowSums(pres_mat[, g2, drop = FALSE])
  
  # Build one 2x2 contingency table per peptide:
  #
  #              Group 1   Group 2
  # Present         x1        x2
  # Absent        n1-x1     n2-x2
  message("Building contingency tables and performing Fisher's Exact Test")
  fisher_tables <- cbind(x1,n1 - x1,x2,n2 - x2)
  
  # Fisher exact test + discrete p-value supports
  fisher_results <- fisher_test_pv(fisher_tables,alternative = "two.sided",exact = TRUE)
  
  # Extract p-values and supports
  p_values <- fisher_results$get_pvalues()
  p_supports <- fisher_results$get_pvalue_supports()
  
  message("Applying discrete FDR")
  # Discrete BH step-down
  dbh <- discrete.BH(p_values,p_supports,direction = "sd")
  
  # Adjusted pvalues
  adj_pval <- dbh$Adjusted
  
  # Percentage of presence in each group
  percent_1 <- x1/n1
  percent_2 <- x2/n2
  
  # Get the LPS
  lps <- log2(adj_pval) * sign(percent_2-percent_1)
  
  resultados <- data.table(
    ID = ids,
    percent_1,
    percent_2,
    P.Val = p_values,
    Adj.P.Val = adj_pval,
    LPS = lps
  )
  
  setnames(
    resultados,
    c("percent_1", "percent_2"),
    paste0("Completeness_", names(grupos))
  )
  
  # Merge with original table
  combined <- input_table %>%
    left_join(
      resultados,
      by = setNames("ID", opt$col_name)
    )
  
  path2save <- gsub("\\.tsv", "_FET.tsv", opt$input)
  
  message("Saving results at ", path2save)
  fwrite(
    combined,
    path2save,
    sep = "\t"
  )
}



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# MORE THAN TWO GROUPS ----
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

if (length(grupos)>2){
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  ## FISHER EXACT TEST MULTICATEGORY (DISTRIBUTION VS NULL HYPOTHESIS) ----
  message("Building contingency tables and performing Global Fisher's Exact Test")
  #get the size of each group
  group_sizes <- sapply(grupos, length)
  #build the  vector support for all the possible number of presences (0,1,2...sum(group_sizes))
  supports <- build_supports_fisher_2xK(group_sizes)
  
  # Function for Fisher exact test in one row
  analizar_pgm <- function(i) {
    
    pres <- pres_mat[i, ]
    # Number of present observations in each group
    presentes <- sapply(group_idx,function(idx) sum(pres[idx]))
    
    # Number of absent observations in each group
    ausentes <- group_sizes - presentes
    
    # 2 x K contingency table
    tabla <- rbind(
      Present = presentes,
      Absent  = ausentes
    )
    
    # Total number of presences
    M <- sum(presentes)
    
    # Fisher exact p-value
    p <- fisher.test(tabla,simulate.p.value = F)$p.value
    
    # Support for this set of margins
    support <- supports[[M + 1L]]
    
    list(
      ID = ids[i],
      p = p,
      support = support
    )
  }
  
  plan(multisession,workers = min(future::availableCores(), length(grupos)))
  #execute the fisher test for each row
  raw_pvalues <- future_lapply(
    seq_len(nrow(pres_mat)),
    analizar_pgm,
    future.seed = TRUE
  )
  
  plan(sequential)
  
  #extract the pvalues and the support vectors
  p_values <- sapply(raw_pvalues,`[[`,"p")
  p_supports  <- lapply(raw_pvalues,`[[`,"support")
  
  message("Applying discrete FDR for Global FET")
  #apply the discrete BH adjust (step-down)
  dbh <- discrete.BH(p_values,p_supports,direction = "sd")
  
  resultados_multicat <- data.table(
    ID = ids,
    P.Val = p_values,
    Adj.P.Val = dbh$Adjusted
  )
  
  # Select proteins significant in the global multicategory test
  sig_idx <- which(dbh$Adjusted < fdr_threshold)
  
  # Presence matrix restricted to globally significant proteins
  pres_mat_sig <- pres_mat[sig_idx, , drop = FALSE]
  
  # IDs of globally significant proteins
  ids_sig <- ids[sig_idx]
  
  if (length(ids_sig) > 0){
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    ## FISHER EXACT TEST ONE GROUP VS THE REST ----
    
    message("Building contingency tables and performing One vs All Fisher's Exact Test")
    message("Applying discrete FDR for One vs All FET")
    #analyze the whole group
    analizar_grupo <- function(g_name) {
      g <- grupos[[g_name]]
      
      # Number of samples in the group and in the whole dataset
      n_g <- length(g)
      N <- ncol(pres_mat_sig)
      n_rest <- N - n_g
      
      # Number of presences in the current group
      x <- rowSums(pres_mat_sig[, g, drop = FALSE])
      # Total number of presences across all samples
      K <- rowSums(pres_mat_sig)
      # Number of presences in the rest of the samples
      x_rest <- K - x
      
      # Build one 2x2 contingency table per PGM
      fisher_tables <- cbind(x,n_g - x,x_rest,n_rest - x_rest)
      
      # Fisher exact test + discrete p-value supports
      fisher_results <- fisher_test_pv(fisher_tables, alternative = "two.sided",exact = TRUE)
      
      # Extract raw p-values and their corresponding supports
      p_values <- fisher_results$get_pvalues()
      p_supports <- fisher_results$get_pvalue_supports()
      
      # Discrete BH step-down
      dbh <- discrete.BH(p_values,p_supports,direction = "sd")
      # Adjusted p-values
      adj_pval <- dbh$Adjusted
      
      # Percentage of samples in the current group where the PGM is present
      completeness_g <- x / n_g
      # Percentage of samples in the rest where the PGM is present
      completeness_rest <- x_rest / n_rest
      
      lps <- log2(adj_pval) * sign(completeness_rest - completeness_g)
      
      
      # -------------------------------------------------------------------------------
      # Results
      # -------------------------------------------------------------------------------
      
      data.table(
        ID = ids_sig,
        grupo = g_name,
        P.Val = p_values,
        Adj.P.Val = adj_pval,
        LPS = lps
      )
    }
    
    
    #merge all the groups analysis
    plan(multisession, workers = min(length(grupos), future::availableCores()))
    resultados_lista <- future_lapply(names(grupos), analizar_grupo, future.seed = TRUE)
    plan(sequential)
    
    resultados_1vsall <- rbindlist(resultados_lista)
    
    one_vs_all_wide <- dcast(
      resultados_1vsall,
      ID ~ grupo,
      value.var = c("P.Val", "Adj.P.Val", "LPS")
    )
    
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    # ONE VS ONE ----
    #create the pairwise comparisons
    comparisons <- combn(names(grupos), 2, simplify = FALSE)
    
    message("Building contingency tables and performing One vs One Fisher's Exact Test")
    message("Applying discrete FDR for One vs One FET")
    #create the function to analyse one comparison
    analizar_comparacion <- function(comp){
      
      g1 <- comp[1] #select the first group of the comparison
      g2 <- comp[2] #select the second group of the comparison
      
      cols_g1 <- grupos[[g1]]
      cols_g2 <- grupos[[g2]]
      
      n1 <- length(cols_g1)
      n2 <- length(cols_g2)
      
      # Number of presences in each group
      x1 <- rowSums(pres_mat_sig[, cols_g1, drop = FALSE])
      x2 <- rowSums(pres_mat_sig[, cols_g2, drop = FALSE])
      
      # Fisher exact test
      fisher_tables <- cbind(x1,n1 - x1,x2,n2 - x2)
      fisher_results <- fisher_test_pv(fisher_tables,alternative = "two.sided",exact = TRUE)
      
      # Raw p-values
      p_values <- fisher_results$get_pvalues()
      # Theoretical p-value supports
      p_supports <- fisher_results$get_pvalue_supports()
      
      # Discrete BH step-down
      dbh <- discrete.BH(p_values,p_supports,direction = "sd")
      adj_pval <- dbh$Adjusted
      
      # Coverage
      completeness_g1 <- x1 / n1
      completeness_g2 <- x2 / n2
      
      # LPS
      lps <- log2(adj_pval) * sign(completeness_g2 - completeness_g1)
      
      # Output
      data.table(
        ID = ids_sig,
        comparison = paste0(g1, "_", g2),
        P.Val = p_values,
        Adj.P.Val = adj_pval,
        LPS = lps
      )
    }
    
    plan(multisession,workers = min(length(comparisons),future::availableCores()))
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
      value.var = c("P.Val", "Adj.P.Val", "LPS"),
      sep = "_"
    )
  } else {
    one_vs_all_wide <- data.frame(ID = character())
    one_vs_one_wide <- data.frame(ID = character())
  }
  
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  # MERGE ALL THE RESULTS
  coverage <- data.table(
    ID = ids,
    sapply(
      names(grupos),
      function(g_name) {
        g <- grupos[[g_name]]
        rowSums(pres_mat[, g, drop = FALSE]) / length(g)
      }
    )
  )
  
  setnames(
    coverage,
    old = names(grupos),
    new = paste0("Completeness_", names(grupos))
  )
  
  all_together <- input_table |>
    left_join(coverage, by = setNames("ID", opt$col_name)) |>
    left_join(resultados_multicat, by = setNames("ID", opt$col_name)) |>
    left_join(one_vs_all_wide, by = setNames("ID", opt$col_name)) |>
    left_join(one_vs_one_wide, by = setNames("ID", opt$col_name))
  
  path2save <- gsub("\\.tsv","_FET.tsv",opt$input)
  message("Saving results at ", path2save)
  fwrite(all_together, path2save, sep = "\t")
}
