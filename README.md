# FisherExactTest
R script implementing Fisher's Exact Test (hypergeometric test) to assess statistical significance in peptide/protein presence/absence patterns. Supports comparisons across two or more groups.

## Dependencies

 - `future.apply` v1.20.2
 - `dplyr` v1.2.1
 - `tidyr` v1.3.2
 - `data.table` 1.18.6.1
 - `stringr` 1.6.0
 - `DiscreteFDR` v2.1.1
 - `DiscreteTests` v0.5.1
 - `optparse` 1.8.2

## Usage

```bash
Rscript FisherExactTest.R -i <input_table> -g <groups_table> -c <col_name> -f <FDR threshold>
```
   > ⚠️ _R bin folder should be in your environment PATH in order to execute R scripts from command line_

## Parameters

* **input_table**: tabular file (TSV) containing peptide/protein-level information, with samples in the columns and peptides/proteins in the rows. This is typically the output of DIA-NN or iSanXoT. Cells corresponding to NA values must be left empty (do not fill them with 0, blank text, or anything similar).
* **groups_table**: TSV file with the sample-to-group correspondence. It must have as many columns as there are groups in input_table, with each sample placed under the column corresponding to its group.
* **col_name**: name of the column in input_table to be used as the row identifier (typically `Precursor.Id` or `Protein.Group` in DIA-NN, or `pgm`, `p`, or `q` in iSanXoT).
* **fdr**: FDR threshold that will be applied to Global Fisher's Exact Test to select which rows are selected to perform One vs All and One vs One.

## Output

After running the script, a TSV file will be generated in the same location as input_table with the `_FET` suffix.

If the number of **groups is 2**, the following columns will be added at the end of the table:
* P.Val: ss
* Adj.P.Val:
* Completeness_<groupA>:
* Completeness_<groupB>:
* LPS: 
