# Common loci underlie natural variation in diverse toxin responses

This repository contains all the processed data and scripts required to recapitulate all results from the paper.

## Supplementary Data Files

### FileS1 - dosedata
Contains results after [*easysorter*](http://github.com/andersenlab/easysorter) processing of the dose response HTA assays for all 16 toxins. Consists of a R dataframe with condition, control, strain, trait, and pruned phenotype.

### FileS2 - allRIAILregressed
Contains residual phenotypic values for each RIAIL for each trait after [*easysorter*](http://github.com/andersenlab/easysorter) processing of linkage mapping HTA assays for all 16 toxins. Consists of an R dataframe with condition, control, strain, trait, and regressed phenotype.

### FileS3 - allAnnotatedLods
Contains the annotated QTL and confidence intervals identified through linkage mapping using the [*linkagemapping*](http://github.com/andersenlab/linkagemapping) package. Consists of an R dataframe with the SNV marker name (marker), chromosome (chr), position (pos), trait, LOD score (lod), genome-wide error rate threshold of significance (threshold), iteration of the forward search (iteration), variance explained of the QTL (var_exp), effect size of the QTL (eff_size), left and right marker and marker position of the 95% confidence interval defined by a 1.5-LOD drop (ci_l_marker, ci_l_pos, ci_r_marker, ci_r_pos).

### FileS4 - scantwosummary
Contains the summary of a two-factor genome scan performed using *Rqtl* for all traits with a significant QTL identified with linkage mapping.

### FileS5 - varianceComponents
Contains the broad-sense heritability estimates as well as additive and interactive components of heritability for each trait.

### FileS6 - correlated_traits
Includes QTL for all mapped traits, the representative trait for overlapping QTL within each toxin, and the correlation coefficient between the representative trait and each correlated trait calculated from the dose response data.

### FileS7 - uniqueQTL
A subset of FileS3 (allAnnotatedLods) for only the 114 distinct QTL identified.

### FileS8 - WGS.vcf
A VCF file for all NILs and CSSs mentioned in this manuscript.

### FileS9 - allNILCSSregressed
Contains residual phenotypic values for each NIL and CSS strain (including parents) for each trait after [*easysorter*](http://github.com/andersenlab/easysorter) processing of NIL/CSS HTA assays for all toxins tested. Consists of an R dataframe with condition, control, strain, trait, and regressed phenotype.

### FileS10 - css_nil_stats
Contains the statistical significance for all pairwise combinations of strains tested for each trait.

### FileS11 - assays_category
Contains the assay categorization for all traits tested with the NIL and CSS strains.
