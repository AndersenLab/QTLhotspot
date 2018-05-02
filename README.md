# Common loci underlie natural variation in diverse toxin responses

This repository contains all the processed data and scripts required to recapitulate all results from the paper.

## Figures
- **20180419_main_figures.pdf** All main figures from the text.
- **supplementing_figures.pdf** All supplemental figures from the text.

## Scripts

- **20180419_main_figures.Rmd** uses the supplemental data provided to reproduce all main figures from the text.
- **supplementing_figures.Rmd** uses the supplemental data provided to reproduce all supplemental figures from the text.


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
Contains the broad-sense heritability estimates (broadsense_H2) as well as additive (narrowsense_h2) and interactive (interaction_VE) components of heritability for each trait. The standard error for additive (narrowsense_h2_SE) and interactive (interaction_SE) components are also available.

### FileS6 - correlated_traits
Includes QTL for all mapped traits (chromosome, condition, and trait), the representative trait for overlapping QTL within each toxin (representative_trait), and the correlation coefficient for the trait and the representative trait calculated from the dose response data (trait_correlation).

### FileS7 - uniqueQTL
A subset of FileS3 (allAnnotatedLods) for only the 114 distinct QTL identified.

### FileS8 - WGS.vcf
A VCF file for all NILs and CSSs mentioned in this manuscript.

### FileS9 - allNILCSSregressed
Contains residual phenotypic values for each NIL and CSS strain (including parents) for each trait after [*easysorter*](http://github.com/andersenlab/easysorter) processing of NIL/CSS HTA assays for all toxins tested. Consists of an R dataframe with condition, control, strain, trait, and regressed phenotype.

### FileS10 - css_nil_stats
Contains the statistical significance for all pairwise combinations of strains tested for each trait calculated using Tukey's Honest Significant Difference Test. Consists of an R dataframe with columns as defined below:
- **condition** - toxin
- **trait** - trait
- **par_sig** - statistical significance between the parents, N2 and CB4856
- **N2_res** - *TRUE* if N2 is more resistant than CB4856 (calculated by comparing strain medians)
- **N2nil_sig_CB** - statistical significance between N2 > CB4856 NIL and CB4856 parent
- **N2nil_recap** - *TRUE* if N2 > CB4856 NIL recapitulates the expected phenotype of N2 (*e.g.* is more resistant than CB4856 if N2 is the resistant parent strain)
- **N2nil_sig_N2** - statistical significance between N2 > CB4856 NIL and N2 parent
- **CBnil_sig_N2** - statistical significance between CB4856 > N2 NIL and N2 parent
- **CBnil_recap** - *TRUE* if CB4856 > N2 NIL recapitulates the expected phenotype of CB4856 (*e.g.* is more resistant than N2 if CB4856 is the resistant parent strain)
- **CBnil_sig_CB** - statistical significance between CB4856 > N2 NIL and CB4856 parent
- **nils_sig** - statistical significance between the NIL strains
- **N2nil_res** - *TRUE* if N2 > CB4856 NIL is more resistant than CB4856 > N2 NIL
- **exp** - name of the assay: IVL - hotspot on the center of chromosome IV; IVR - hotspot on right of chromosome IV; V - hotspot on center of chromosome V (NIL); CSSV - hotspot on center of chromosome V (CSS)
- **chr** - chromosome of assay (IV or V)
- **N2nil** - name of the N2 > CB4856 NIL for this assay
- **CBnil** - name of the CB4856 > N2 NIL for this assay
- **N2** - median phenotype for N2 parent strain
- **CB** - median phenotype for CB4856 parent strain
- **N2nil_value** - median phenotype for N2 > CB4856 NIL strain
- **CBnil_value** - median phenotype for CB4856 > N2 NIL strain

### FileS11 - assays_category
Contains the assay categorization for all traits tested with the NIL and CSS strains. Consists of an R dataframe with condition, trait, hotspot/experiment (exp, IVR = hotspot on the right of IV, IVL = hotspot on the center of IV, V = hotspot on the center of V (NIL), CSSV = hotspot on center of V (CSS)), category defined by just the NIL or the CSS assay (primary_category) and category defined by combining the NIL and CSS assay, if applicable (only traits that mapped to chromosome V) (secondary_category).
