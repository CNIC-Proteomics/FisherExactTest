# FisherExactTest
R script implementing Fisher's Exact Test (hypergeometric test) to assess statistical significance in peptide/protein presence/absence patterns. Supports comparisons across two or more groups.

## Dependencies

 - `future.apply` v1.20.2
 - `dplyr` v1.2.1
 - `tidyr` v1.3.2
 - `data.table` 1.18.6.1
 - `stringr` 1.6.0
 - `DiscreteFDR` v2.1.1
 - `optparse` 1.8.2

## Usage

```bash
Rscript FisherExactTest.R -i <input_table> -g <groups_table> -c <col_name>
```
   > ⚠️ _R bin folder should be in your environment PATH in order to execute R scripts from command line_

## Parameters

* **input_table**: tabular file (TSV) containing peptide/protein-level information, with samples in the columns and peptides/proteins in the rows. This is typically the output of DIA-NN or iSanXoT. Cells corresponding to NA values must be left empty (do not fill them with 0, blank text, or anything similar).
* **groups_table**: TSV file with the sample-to-group correspondence. It must have as many columns as there are groups in input_table, with each sample placed under the column corresponding to its group.
* **col_name**: name of the column in input_table to be used as the row identifier (typically `Precursor.Id` or `Protein.Group` in DIA-NN, or `pgm`, `p`, or `q` in iSanXoT).

## Output

After running the script, a TSV file will be generated in the same location as input_table with the `_FET` suffix. This file will contain the same columns as input_table plus the following columns:

* **P.Val.Enrich_\<group\>**: contains the unadjusted p-value determining whether the corresponding group is enriched (i.e., has greater presence in that group).
* **Signif_Enrich_\<group\>**: binary variable that takes the value TRUE when the enrichment p-value for that group passes an FDR of 0.05.

If there are more than two groups, the following additional columns will also be provided:

* **P.Val.Deple_\<group\>**: contains the unadjusted p-value determining whether the corresponding group is depleted (i.e., has greater absence in that group).
* **Signif_Deple_\<group\>**: binary variable that takes the value TRUE when the depletion p-value for that group passes an FDR of 0.05.
* **P.Val**: p-value for the "multi-category" Fisher Exact Test (i.e., it compares the presence/absence distribution against the null distribution and determines whether there is a different pattern in our distribution).
* **P.Mid**: the mid-p-value, a continuous approximation for p-values from discrete tests (such as the Fisher Exact Test), which allows p-value correction methods to be applied.
* **Adj.P.Val**: mid-p-value adjusted using the Benjamini-Hochberg method.
* **P.Val.Enrich/Deple_\<groupA\>_\<groupB\>**: for those rows (peptides/proteins) with an Adj.P.Val below 0.05, or with any Signif_Enrich/Deple_\<group\> column equal to TRUE, a Fisher Exact Test is performed across all possible pairwise group combinations. For each comparison, a p-value is provided for the enrichment of the first group (A) over the second (B), as well as for the depletion of the first group (A) relative to the second (B).
* **Signif_Enrich/Deple_\<groupA\>_\<groupB\>**: binary variable that takes the value TRUE when the enrichment/depletion p-value for that groupA-groupB comparison passes an FDR of 0.05.

