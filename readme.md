# Demo for Self-RRM

This repository contains work-in-progress demo code for the Self-reference region method (Self-RRM) for analysing DCE-MRI.
The Self-RRM relies on the reference region model, such as [CERRM](https://github.com/MPUmri/ERRM), where the reference tissue parameters are estimated with [RRIFT](https://github.com/MPUmri/RRIFT).
What sets Self-RRM apart from previous work is that the reference tissue is automatically identified from the tissue of interest (i.e. tumour), whereas the previous approach requires a manually identified reference region (e.g. muscle).

This demo applies the Self-RRM to DCE-MRI data acquired from a patient with glioblastoma multiforme.

- b01_process.m
    + Fits the tumour data using:
        * Self-RRM
        * Reference region model using muscle as the reference region (conventional approach)
        * Extended Tofts model using the measured arterial input function
- b02_showResults.m
    + Shows parameter maps for the three fits, along with the concordance correlation coefficients


## Acknowledgements

- Glibolastoma multiforme data was obtained from the Cancer Genome Atlas Glioblastoma Multiforme [(TCGA-GBM)](https://wiki.cancerimagingarchive.net/display/Public/TCGA-GBM) collection on the Cancer Imaging Archive (TCIA)
- Code for NMF-based hierarchical clustering was written by [Nicolas Gillis](https://sites.google.com/site/nicolasgillis/code)
