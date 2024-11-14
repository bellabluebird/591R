# Assignment Counts Data

## Notes for this assignment

For this and all future assignments, you will now have to code your markdown 
file in addition to the main.R functions. The report.Rmd file still has the same 
background and assignment information included and the R code blocks are created, 
just like in previous assignments, but from now on they will devoid of code.

Don't worry, there are instructions in report.Rmd to guide you through the code 
blocks!

## Problem Statement

When dealing with gene counts in an mRNA-seq dataset, it is important to
normalize the data before performing any analyses so you can make accurate
comparisons of gene expression between samples. There are a number of different
methods you can use to accomplish and the method you use will depend on the
kinds of samples you have and the analyses you want to perform. This first part
of this assignment will guide you through two commonly used methods: Counts Per
Million and DESeq2 Normalization.

Once you normalize your data you can do exploratory analysis and construct the
appropriate graphs. As bioinformaticians, you will need to be able present your
data in tidy reports with generated plots- sometimes between blocks of text. R
Markdown is a good tool to accomplish this with, allowing you to display tables
and plots alongside the bulk of your writing. You've had some exposure to R
Markdown in previous assignments, now you will start writing your own code to
generate the visualizations for your report. It will be important to remember
that you need to keep your code separate depending on its function: the function
declarations and implementation- your code that will work behind-the-scenes to
perform data manipulation and process inputs without user interaction on the
"back end"- should stay in your `main.R` file while your calls to construct and
display the tibbles and visuals created in the back end should be put in your
`report.Rmd` file- despite the user not directly manipulating or editing the
outputs you generate, they are still interacting with the display by viewing it
and therefore the lines of code that call the functions to create these outputs
belong in your "front end."

## Required Readings

Please read sections 9.7.5 to 9.7.9. This assignment will also draw on concepts from chapters 6, 7, 8.

## Learning Objectives
- Data normalization methods (CPM, DESeq2)
- Plotting data in Tidyverse
- Displaying results and compiling reports in R Markdown
- Application of Tidyverse functions

## Skill List
- DESeq normalization, referencing the linked Bioconductor vignette if needed
- Tibble manipulation
- Creating different types of graphs in ggplot2
- Running PCA

## Instructions
You will be given the following files:

```
main.R
test_main.R
report.Rmd
verse_counts.tsv
example_report.html

```
Like previous assignments, the bulk of the work will be done in the skeleton
file `main.R`, but unlike previous assignments the provided `report.Rmd` is also
a skeleton file. Both documents will provide the information needed to complete
the assignment. Please ensure that your inputs and outputs match the
specifications listed- or at the very least can properly handle inputs as
instructed and return the expected output since that is how your functions will
be tested.

To help you complete the assignment, you are being provided a sample report and
a set of test functions that are similar to the ones we will use to grade. The
sample report is provided in the appropriately named `sample_report.html` and
the tests for you to use can be found in `test_main.R`. Please note that while
the sample report does not display any code, it is okay and- preferred, really-
for your report to _display the code blocks_ like they have done in previous
assignments.

### Tasks

1. Implement the functions as described in `main.R`.  
2. [optional but highly recommended] Test your functions using `test_main.R
3. Fill in `report.Rmd` as instructed. You will need the functions you implemented in `main.R`
4. Knit `report.Rmd` and create your `report.html`

## Deliverables

1. `main.R` with all of the functions completed
2. `report.Rmd` filled in as instructed
2. `report.html` knitted from your `report.Rmd`


## Function Details

### 1. `read_data()`

Load a tsv located at a specified location

Input (1): string path to file

Output: (g x m) tibble

Details: This function will be tested for handling various input strings and
returning the proper output: dimensions, column names, if column 'gene' is
located in the first column, and output type tibble. Your test_main will include
additional testing of output column types and lack of row names for your
reference, to help catch errors early on.

### 2. `filter_zero_var_genes()`

Filter out genes with zero variance

Input (1): (g x m) tibble

Output: (n x m) tibble

Details: This function will be tested for handling an input (g x m) tibble and
returning the proper output: column names, if column 'gene' is located in the
first column, names of 'genes' returned, and output type tibble. Your test_main
will include additional testing for row consistency- ensuring that the sample
data still correlates with the gene names- for your reference and to help catch
errors early on

### 3. `timepoint_from_sample()`

Extract time point information from sample name

Input (1): string (length 5) of sample name in format `v[A-Z][a-z,1-9]_[1-9]`
(In other words: `v[α][β]_[γ]`, where α is any capital modern English letter, β
is any lower case modern English letter OR any Arabic number 0-9, and γ is any
Arabic number from 0-9)

Output: string (length 2) of substring `[A-Z][a-z,1-9]` from sample name:

Details: This function will be tested for handling various strings of length 5
in the form `v[A-Z][a-z,1-9]_[1-9]` and outputting the proper string, preserving
letter case where letter case is provided.

### 4. `sample_replicate()`

Grab sample replicate number from sample name

Input (1): string (length 5) of sample name in format
`v[A-Z][a-z,1-9]_[1-9]` (In other words: `v[α][β]_[γ]`, where α is any capital
modern English letter, β is any lower case modern English letter OR any Arabic
number 0-9, and γ is any Arabic number from 0-9)

Output: string (length 1 of substring `[1-9]` from sample name:
`v[A-Z][a-z,1-9]_[1-9]` (In other words: [γ] from `v[α][β]_[γ]`)

Details: This function will be tested for handling various strings of length 5
in the form `v[A-Z][a-z,1-9]_[1-9]` and outputting the proper character string.

### 5. `meta_info_from_labels()`

Generate sample-level metadata from sample names and stores the data into a
tibble. Will include columns named "sample", "timepoint", and "replicate" that
store sample names, sample time points, and sample replicate, respectively.


Input (1): Character vector of length `_S_ ` of sample names with column names "sample", "timepoint", and "replicate"

Output: a `(_S_ x 3)` tibble

Details: This function will be tested for handling of a character vector of
length `_S_`, where each element in the vector is a string with a length of 5 in
the form of `v[A-Z][a-z,1-9]_[1-9]`, and properly outputting a `(_S_ x 3)`
tibble with columns named "sample", "timepoint", and "replicate" and rows that
correspond with the input samples. Your test_main will include additional
testing for the order of elements in column 'sample''s correspondence with the
order of elements in the input vector for your reference. The column types of
your output tibble will not be tested.

### 6. `get_library_size()`

Calculate total read counts for each sample in a counts dataset.

Input (1): a (n x m) tibble of raw read counts

Output: tibble of read totals from each sample. a tibble can be `(1 x _S_)` 
with sample names as columns names OR `(_S_ x 2)` where sample name is in the 
first column and library size is the second column

Details: This function will be tested for the return of a tibble that have sample 
names which correspond with the appropriate library size.

### 7. `normalize_by_cpm()`

Normalize raw counts data to counts per million using (counts) /
(sample_library_size) * 10^6

Input (1): a (n x m) tibble of raw read counts

Output: a (n x m) tibble with read count normalized to counts per million

Details: This function will be tested to handle a (n x m) tibble. Its output
will be tested for dimensions, column names, location of string and numeric
column(s), performance of the desired equation on numeric columns, and that gene
names still correspond to their rows.

### 8. `deseq_normalize()`

Normalize raw counts data using DESeq2

Input (1): a (n x m) tibble of raw read counts

Output: a (n x m) tibble of DESeq2 normalized counts data

Details:This function will be tested to handle a (n x m) tibble. Its output will
be tested for dimensions, column names, location of string and numeric
column(s), performance of the desired equation on numeric columns, and that gene
names still correspond to their rows.

### 9. `plot_pca()`

Input (3): a `(n x _S_)` tibble of data, a `(_S_ x 3)` tibble of sample-level
meta information, and a string

Output: a ggplot scatter plot showing each sample, with PC1 on x-axis and PC2 on y-axis

Details: The output of this will be tested for the appropriate test run and PCs
used to plot. It may be visually inspected as part of your grade

### 10. `plot_sample_distributions()`

Input (3): a `(n x _S_)` tibble of data, a boolean to determine whether to scale
the 'y' axis to log10 values, and a string

Output: a ggplot boxplot that shows gene count distributions

Details: This function will be tested on a `(n x _S_)` tibble. It will be tested
for functionality of its inputs, handling of data, expected graph elements, and
graph type. It may also be visually inspected as part of your grade

### 11. `plot_variance_vs_mean()`

Input (3): a `(n x _S_)` tibble of data, a boolean to determine whether to scale
the 'y' axis to log10 values, and a string

Output: a ggplot scatter plot where the x-axis is the rank of gene ordered by
mean count over all samples, and the y-axis is the observed variance of the
given gene. Each dot should have their transparency increased. The scatter plot
should also be accompanied by a line representing the average mean and variance
values


Details: This function will be tested on a `(n x _S_)` tibble. It will be tested
for functionality of its inputs, handling of data, expected graph elements, and
graph type. It may also be visually inspected as part of your grade

## Hints

* Make sure you don't have any rownames on any of the tibbles your function(s)
 return. You should get a warning in your R console if you do.

* Sometimes Tydyverse is the best tool for the job, sometimes it isn't

* Some tasks are easier to accomplish with dataframes than tibbles. Not everyone
  will use these methods, but if you do you are welcome to use dataframes within
  a function as long as the inputs and outputs are the ones specified (returning
  a dataframe instead of a tibble may points).

* The Bioconductor vignette for DESeq2 is linked in part 3 of your `report.Rmd`.
  [It can also be found
  here](http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#data-transformations-and-visualization))

* In a previous version of this assignment, the output of
  `filter_zero_var_genes()` was transformed from `(n x m)` to `(n x _s_)` and
  called the `counts_matrix`. One of the TAs thought use of 'matrix' would be
  confusing (since we want you to use datasets of type `tibble`) so it was
  re-worded, but this piece of trivia might be helpful for interpreting the
  Bioconductor vignette- or if you see any sort of 'matrix' referenced in
  `report.Rmd` because it got overlooked during editing.

* There may be occasions where you will need to reshape your data from a wide
 format to a long format

* Symbols used in `main.R`~ `g`: initial number of Genes, `m`: initial number of
  columns expected when you import `verse_counts.tsv`, `n`: number of genes
  expected after you filter in part 1b, `_S_`: number of Samples
