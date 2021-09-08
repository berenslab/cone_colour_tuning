The used SBI method is a version of SNPE-C developped by the Mackelab.
The forked and modified version for adaptive observations cand be found here:
https://github.com/berenslab/delfi/tree/adaptive_obs

The following notebooks depend on this repo. 
Additionally the analysis depends on the draws from the posteriors, which are not included in this repo (>1GB), but they can be send upon request.

## for main Fig. 4 and S4
- data: "data_reextracted.hdf5"
- for SNPE fitting: model_fit_SNPE_v0-mode_x.ipynb
- for analyzing posteriors:
    model_SNPE_analyse_v2.ipynb (posterior checking and deciding which to choose)
- for sampling from posterior, and generating simulations:
    model_SNPE_evaluation_v2.ipynb
- for plotting and final comparison: best_model_comparison_v2.ipynb



## for regionwise models (cf. STable S1):
- model_fit_SNPE_v0-mode_x-regionwise.ipynb (for fitting)
- model_fit_SNPE_analyse_v2-regionwise.ipynb (for finding best posterior)
- model_fit_SNPE_evaluation_v2-regionwise.ipynb (for sampling and evaluating model from posterior)
- best_model_comparison_v2-regionwise.ipynb (for plotting and final comparison)
