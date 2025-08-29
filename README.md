
# IPMsquared: An integrated integral projection model to disentange size-structured harvest and natural mortality

<img src="greencrab_drawing.png" align="right" height="200" style="float:right; height:200px;" dpi="700"/>

This repository contains code and data for running an IPM^2, a
combination of an integrated population model and integral projection
model, to understand size-structured population dynamics of the invasive
European green crab. This repository also contains the manuscript and
supplemental information associated with the manuscript.

## File structure and contents

    IPMsquared
    │   README.md  
    │
    └───code
    │   │
    │   │ IPM_code.R # main code file for running the IPM^2
    │   │ IPM_code_seadrift # supplemental code file for running the IPM^2 with the alternative mark-recapture dataset (D2)    
    │   │ IPM_simulations.R # code for running stochastic simulations with posterior samples to approximate the stable size structure in response to varying removal efforts
    │   │ posterior_predictive_check.R # code for generating posterior predictive samples and calculating Bayesian p-values with the posterior generated in the primary analysis
    │   │ posterior_predictive_check_seadrift.R # code for generating posterior predictive samples and calculating Bayesian p-values with the posterior generated in the supplemental analysis with the alternative mark-recapture dataset (D2)
    │   │
    │   └─────model_selection
    │         │   IPM_code_modelX.R # code for conducting the WAIC calculation with model X in the model selection procedure (Appendix 6)
    │         │   savedsamples_modelX.rds # saved posterior samples associated with model X in the model selection procedure (Appendix 6)
    │   
    └───data
    │   │
    │   └─────model_data # all data and constants used in the IPM^2 model; description of each file found in folder README
    │   │
    │   └─────D1_catch_data # catch data at D1 used to generate Figure 1B
    │   │
    │   └─────D2_growth_data # size-at-age data (D2)
    │   │
    │   └─────D3_mark_recapture # mark-recapture data (D3) and code used to generate .rds files found in /model_data
    │   │
    │   └─────drayton_harbor_map_data # data used to generate map of Drayton Harbor in Appendix 1
    │   │
    │   └─────posterior_predictive_check # data associated with the posterior predictive check with the primary analysis and /seadrift folder contains data associated with the posterior predictive check with the alternative mark-recapture dataset (D3)
    │   │
    │   └─────posterior_samples # posterior samples from the primary IPM^2 model fit and from the supplemental model fit with the alternative mark-recapture dataset (D3)
    │   │
    │   └─────simulations # data associated with the stochastic simulations used to forecast equilibrium population size
    │
    └───figures # code used to generate primary and supplemental figures, as well as resulting .png figure files
    │
    └───manuscript
        │
        └─────IPMsquared_main.tex # latex file of main manuscript text
        │
        └─────Main_text.pdf # pdf file of main manuscript text
        │
        └─────title_page.tex # latex file of manuscript title page
        │
        └─────Title_page.pdf # pdf file of manuscript title page
        │
        └─────maintext.bib # bibliography file for main text
        │
        └─────AppendixX.Rmd # Markdown file to generate Appendix X
        │
        └─────AppendixX.pdf # pdf file of Appendix X
        │
        └─────appendixX.bib # bibliography file associated with Appendix X

## Bayesian implementation

The model is specified using NIMBLE (download instructions given here:
<https://r-nimble.org/download>), which adopts and extends BUGS as a
Bayesian modeling language. We also use the Markov chain Monte Carlo
algorithms written in the NIMBLE software to generate posterior samples,
and use the NIMBLE model objective to generate posterior predictive
samples as part of our model checking procedure.

We also make use of NIMBLE functions (see the [main model code
file](https://github.com/abigailkeller/IPMsquared/blob/main/code/IPM_code.R))
to make our code more modular and readable and to create custom
distributions (i.e., the Dirichlet-Multinomial).

## Citation

Releases of this repository are archived and issued DOIs via Zenodo:
[https://zenodo.org/records/10462269](https://zenodo.org/records/10572340)

The citation of the latest version is:

## Contact

Abby Keller: <agkeller@berkeley.edu>
