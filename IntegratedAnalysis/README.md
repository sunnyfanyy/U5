## Time-series expression analysis 

**For the analysis of expression levels and relative changes under dark-to-light transition at different time points, it is noted that the Dark-to-Light 6-hour sample (DTL-6h) was prepared using an rRNA-depleted library construction method, while the remaining dark-to-light time points (12h, 24h, 48h) utilized mRNA-seq library preparation. To ensure methodological consistency, the corresponding Dark control samples for each library type were used as baselines. Fold changes in gene expression at each dark-to-light transition time point were calculated relative to their respective Dark controls based on the matched library preparation protocols.**

**Key Clarifications:**

1. Methodological Alignment:
   - The DTL-6h sample (rRNA-depleted) and its Dark control were processed using the same rRNA-depletion protocol.
   - For mRNA-seq time points (12h, 24h, 48h), Dark controls were generated via mRNA-seq to maintain comparability.
   - This approach minimizes technical biases arising from library preparation differences.
2. Rationale for Control Selection:
   - Cross-method comparisons (e.g., rRNA-depleted vs. mRNA-seq controls) could introduce confounding factors due to differences in transcript capture efficiency and background noise.
   - Matched controls ensure that observed expression changes reflect biological responses rather than methodological artifacts.
3. Statistical Implications:
   - Fold change calculations were normalized within each library type to account for potential batch effects.

This strategy adheres to best practices in transcriptomic studies where heterogeneous protocols are employed across experimental conditions.