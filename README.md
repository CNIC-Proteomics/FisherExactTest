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

* `input_table`: tabular file (TSV) containing peptide/protein-level information, with samples in the columns and peptides/proteins in the rows. This is typically the output of DIA-NN or iSanXoT. Cells corresponding to NA values must be left empty (do not fill them with 0, blank text, or anything similar).
* `groups_table`: TSV file with the sample-to-group correspondence. It must have as many columns as there are groups in input_table, with each sample placed under the column corresponding to its group.
* `col_name`: name of the column in input_table to be used as the row identifier (typically `Precursor.Id` or `Protein.Group` in DIA-NN, or `pgm`, `p`, or `q` in iSanXoT).
* `fdr`: FDR threshold that will be applied to Global Fisher's Exact Test to select which rows are selected to perform One vs All and One vs One.

## Output

After running the script, a TSV file will be generated in the same location as input_table with the `_FET` suffix.

If the number of **groups is 2**, the following columns will be added at the end of the table:
* `Completeness_<groupA>`: percentage of samples in group A containing information for the given feature.
* `Completeness_<groupB>`: percentage of samples in group B containing information for the given feature.
* `P.Val`: two-sided Fisher's exact test p-value assessing whether there is an association between the presence/absence pattern and the defined groups.
* `Adj.P.Val`: p-value adjusted for multiple testing using the discrete FDR step-down procedure.
* `LPS`: signed log2-transformed adjusted p-value. Positive values indicate that the feature is more represented in group A, whereas negative values indicate that the feature is more represented in group B. A value of 0 indicates no difference between groups.

If there are **more than 2 groups**, the following columns will be added at the end of the table:
* `Completeness_<group>`: percentage of samples in the corresponding group containing information for the given feature.
* `P.Val`: two-sided Fisher's exact test p-value assessing whether there is an association between the presence/absence pattern and the defined groups.
* `Adj.P.Val`: p-value adjusted for multiple testing using the discrete FDR step-down procedure.
  
For features with a significant Adj.P.Val (according to the user-defined FDR threshold), the following additional columns are included:

* `P.Val_<group>`: two-sided Fisher's exact test p-value comparing the presence/absence pattern in the specified group against all other groups combined.
* `Adj.P.Val_<group>`: p-value from P.Val_<group> adjusted for multiple testing using the discrete FDR step-down procedure.
* `LPS_<group>`: signed log2-transformed adjusted p-value. Positive values indicate that the feature is more represented in the specified group than in the other groups, whereas negative values indicate that the feature is underrepresented in the specified group. A value of 0 indicates no difference.
* `P.Val_<groupA>_<groupB>`: two-sided Fisher's exact test p-value assessing whether there is a significant difference in the presence/absence pattern between group A and group B.
* `Adj.P.Val_<groupA>_<groupB>`: p-value from P.Val_<groupA>_<groupB> adjusted for multiple testing using the discrete FDR step-down procedure.
* `LPS_<groupA>_<groupB>`: signed log2-transformed adjusted p-value. Positive values indicate that the feature is more represented in group A than in group B, whereas negative values indicate that the feature is more represented in group B than in group A. A value of 0 indicates no difference.
