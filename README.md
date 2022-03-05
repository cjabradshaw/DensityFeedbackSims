# Density-Feedback Simulations

<img align="right" src="www/Thylacoleo carnifex.png" alt="thylacine" width="300" style="margin-top: 20px">

Stochastic simulations of population abundance with known component density feedback on survival to test for ability to return ensemble feedback signal

R code accompanying paper (pre-print; currently in review elsewhere):

<a href="http://scholar.google.com.au/citations?sortby=pubdate&hl=en&user=1sO0O3wAAAAJ&view_op=list_works">Bradshaw, CJA</a>, <a href="https://scholar.google.com.au/citations?hl=en&user=-BSGg1MAAAAJ">S Herrando-Pérez</a>. 2021. Density-independent processes decouple component and ensemble density feedbacks. <em><strong>bioRχiv</strong></em> doi:<a href="http://doi.org/10.1101/2021.09.19.460939">10.1101/2021.09.19.460939</a>

## Abstract
Analysis of long-term trends in abundance provide insights into population dynamics. Population growth rates are the emergent interplay of fertility, survival, and dispersal, but the density feedbacks on some vital rates (component) can be decoupled from density feedback on population growth rates (ensemble). However, the mechanisms responsible for this decoupling are poorly understood. We simulated component density feedbacks on survival in age-structured populations of long-living vertebrates and quantified how imposed nonstationarity (density-independent mortality and variation in carrying capacity) modified the ensemble feedback signal estimated from logistic-growth models to the simulated abundance time series. The statistical detection of ensemble density feedback was largely unaffected by density-independent processes, but catastrophic and proportional mortality eroded the effect of density-dependent survival on ensemble-feedback strength more strongly than variation in carrying capacity. Thus, phenomenological models offer a robust approach to capture density feedbacks from nonstationary census data when density-independent mortality is low.

<br>
Prof <a href="http://scholar.google.com.au/citations?sortby=pubdate&hl=en&user=1sO0O3wAAAAJ&view_op=list_works">Corey J. A. Bradshaw</a> <br>
<a href="http://globalecologyflinders.com" target="_blank">Global Ecology</a>, <a href="http://flinders.edu.au" target="_blank">Flinders University</a>, Adelaide, Australia <br>
September 2021 (updated February 2022) <br>
<a href=mailto:corey.bradshaw@flinders.edu.au>e-mail</a> <br>

## INSTRUCTIONS

1. Run initial base population models (<code>baseModels.R</code>; in the <a href="https://github.com/cjabradshaw/DensityFeedbackSims/tree/main/basemodels"><em>basemodels</em></a> folder)

<img align="right" src="www/Palorchestes.png" alt="thylacine" width="250" style="margin-top: 20px">
Models are based on the demography of the following genera/species:
- <strong>VOMBATIFORM HERBIVORES</strong>: <a href="https://australian.museum/learn/australia-over-time/extinct-animals/diprotodon-optatum/"><i>Diprotodon</i></a><sup>†</sup>, <a href="https://australian.museum/learn/animals/mammals/palorchestes-azeal/"><i>Palorchestes</i></a><sup>†</sup>, <a href="http://www.megafauna.com.au/view/megafauna/zygomaturus-trilobus"><i>Zygomaturus</i></a><sup>†</sup>, <a href="http://www.seamonsters.qm.qld.gov.au/sitecore/content/QM%20Micro/Project%20DIG/Home/research/tropical-megafauna/species/phascolonus"><i>Phascolonus</i></a><sup>†</sup>, <a href="https://australian.museum/learn/animals/mammals/common-wombat/"><i>Vombatus ursinus</i></a>
<img align="left" src="www/Simosthenurus occidentalis.png" alt="thylacine" width="250" style="margin-top: 20px">
- <strong>MACROPODIFORM HERBIVORES</strong>: <a href="http://www.seamonsters.qm.qld.gov.au/sitecore/content/QM%20Micro/Project%20DIG/Home/research/tropical-megafauna/species/protemnodon"><i>Protemnodon</i></a><sup>†</sup>, <a href="https://australian.museum/learn/animals/mammals/red-kangaroo/"><i>Osphranter rufus</i></a>, <a href="https://en.wikipedia.org/wiki/Sthenurus"><i>Sthenurus</i></a><sup>†</sup>, <a href="http://www.megafauna.com.au/view/megafauna/simosthenurus-occidentalis"><i>Simosthenurus</i></a><sup>†</sup>, <a href="https://australian.museum/learn/australia-over-time/extinct-animals/procoptodon-goliah/"><i>Procoptodon</i></a><sup>†</sup>, <a href="https://en.wikipedia.org/wiki/Sthenurinae"><i>Metasthenurus</i></a><sup>†</sup>, <a href="https://bie.ala.org.au/species/urn:lsid:biodiversity.org.au:afd.taxon:4bd05bcb-614d-40b0-b81f-75ac14ea4afd"><i>Notamacropus</i></a>
- <strong>LARGE BIRDS</strong>: <a href="https://australian.museum/learn/australia-over-time/extinct-animals/genyornis-newtoni/"><i>Genyornis</i></a><sup>†</sup>, <a href="https://www.birdlife.org.au/bird-profile/emu"><i>Dromaius novaehollandiae</i></a>, <a href="https://www.birdlife.org.au/bird-profile/australian-brush-turkey"><i>Alectura lathami</i></a>
- <strong>CARNIVORES</strong>: <a href="https://australian.museum/learn/animals/mammals/tasmanian-devil/"><i>Sarcophilus</i></a>, <a href="https://australian.museum/learn/australia-over-time/extinct-animals/the-thylacine/"><i>Thylacinus</i></a><sup>†</sup>, <a href="https://australian.museum/learn/animals/mammals/thylacoleo-carnifex/"><i>Thylacoleo</i></a><sup>†</sup>, <a href="https://australian.museum/learn/animals/mammals/spotted-tailed-quoll/"><i>Dasyurus</i></a>
- <strong>MONOTREMES</strong>: <a href="https://www.artistwd.com/joyzine/australia/articles/megafauna/megalibgwilia_ramsayi.php"><i>Megalibgwilia</i></a><sup>†</sup>, <a href="https://www.bushheritage.org.au/species/echidna"><i>Tachyglossus</i></a>

<sup><small>† extinct</small></sup>

More details and justification of the model components can be found in:
- Bradshaw, CJA, CN Johnson, J Llewelyn, V Weisbecker, G Strona, F Saltré. 2021. <a href="http://doi.org/10.7554/eLife.63870">Relative demographic susceptibility does not explain the extinction chronology of Sahul’s megafauna</a>. <em>eLife</em> 10: e63870. doi:10.7554/eLife.63870
- and its associated Github repository <a href="https://github.com/cjabradshaw/MegafaunaSusceptibility">here</a>
    
2. Run projection scenarios using outputs from Step 1 (the following scripts are in the <a href="https://github.com/cjabradshaw/DensityFeedbackSims/tree/main/projectionscenarios"><em>projectionscenarios</em></a> folder):
    - <code>stableCatastrophe.R</code>: stable, fixed carrying capacity (<em>K</em>) with generationally scaled catastrophe
    - <code>pulseCatastrophe.R</code>: stable, fixed <em>K</em> with 90% mortality pulse disturbance at 20<em>G</em> (generations); generationally scaled catastrophe
    - <code>r001-01Catastrophe.R</code>: stable, fixed <em>K</em> with density-independent mortality to cause mean population growth rate (<em>r</em>) to be -0.001 or -0.01 (user choice); generationally scaled catastrophe
    - <code>KstochCatastrophe.R</code>: stable, mean stochastic <em>K</em> (constant variance); generationally scaled catastrophe
    - <code>KstochVarIncCatastrophe.R</code>: stable, mean stochastic <em>K</em> (variance doubles over projection interval); generationally scaled catastrophe
    - <code>KstochDeclCatastrophe.R</code>: declining (-0.001), mean stochastic <em>K</em> (constant variance); generationally scaled catastrophe
    - <code>NoDFstableCatastrophe.R</code>: stable, mean trajectory with no density feedback (no <em>K</em>)
3. Import resulting .RData files from Step 2 to estimate evidence and strength of ensemble density feedback (<code>logisticGrowthFit.R</code>; in the <a href="https://github.com/cjabradshaw/DensityFeedbackSims/tree/main/ensembleDF"><em>ensembleDF</em></a> folder)
4. Import resulting .RData files from Step 2 to estimate stationarity metrics (the following scripts are in the <a href="https://github.com/cjabradshaw/DensityFeedbackSims/tree/main/stationarity"><em>stationarity</em></a> folder):
    - <code>StatStableCatastrophe.R</code>: stable, fixed <em>K</em> with generationally scaled catastrophe
    - <code>StatPulseCatastrophe.R</code>: stable, fixed <em>K</em> with 90% mortality pulse disturbance at 20<em>G</em>; generationally scaled catastrophe
    - <code>Statr001Catastrophe.R</code>: stable, fixed <em>K</em> with density-independent mortality to cause mean <em>r</em> = -0.001; generationally scaled catastrophe
    - <code>Statr01Catastrophe.R</code>: stable, fixed <em>K</em> with density-independent mortality to cause mean <em>r</em> = -0.01; generationally scaled catastrophe
    - <code>StatKstochCatastrophe.R</code>: stable, mean stochastic <em>K</em> (constant variance); generationally scaled catastrophe
    - <code>StatKstochVarIncCatastrophe.R</code>: stable, mean stochastic <em>K</em> (variance doubles over projection interval); generationally scaled catastrophe
    - <code>StatKstochDeclCatastrophe.R</code>: declining, (-0.001) mean stochastic <em>K</em> (constant variance); generationally scaled catastrophe
    - <code>StatStableNoDFCatastrophe.R</code>: stable, mean trajectory with no density feedback (no <em>K</em>)
<img align="right" src="www/Diprotodon 2.png" alt="thylacine" width="300" style="margin-top: 20px">
5. Check a few assumptions (in the <a href="https://github.com/cjabradshaw/DensityFeedbackSims/tree/main/assumptionscheck"><em>assumptionscheck</em></a> folder)
    - <code>effectIncVarjuvSurv.R</code>: this script increases the variation in juvenile survival relative to adults (3x, declining linearly to equal adult variation by age at first reproduction) (<em>Diprotodon</em> only)
    -  <code>DFfertDFsurv.R</code>: this script adds a feedback mechanism to fertility in addition to survival to examine how dispersing the feedback signal among > 1 vital rates affects the phenomenological model results (<em>Diprotodon</em> only)


## Requires the following libraries
- <code>dplyr</code>
- <code>plotly</code>
- <code>expss</code>
- <code>car</code>
- <code>Hmisc</code>
- <code>cluster</code>
- <code>bootstrap</code>
- <code>data.table</code>
- <code>ggplot2</code>
- <code>ggridges</code>
- <code>ggpubr</code>
- <code>biostat</code>
- <code>reshape2</code>

## and following source-code file
- <code>matrixOperators.R</code> (in the <a href="https://github.com/cjabradshaw/DensityFeedbackSims/tree/main/source"><em>source</em></a> folder)




