# Density-Feedback Simulations
Stochastic simulations of population abundance with known component density feedback on survival to test for ability to return ensemble feedback signal

R code accompanying paper:

<a href="http://scholar.google.com.au/citations?sortby=pubdate&hl=en&user=1sO0O3wAAAAJ&view_op=list_works">BRADSHAW, CJA</a>, <a href="https://scholar.google.com.au/citations?hl=en&user=-BSGg1MAAAAJ">S HERRANDO-PÉREZ</a>. In review. Component and ensemble density feedbacks decoupled by density-independent processes. 

Prof <a href="http://scholar.google.com.au/citations?sortby=pubdate&hl=en&user=1sO0O3wAAAAJ&view_op=list_works">Corey J. A. Bradshaw</a> <br>
<a href="http://globalecologyflinders.com" target="_blank">Global Ecology</a>, <a href="http://flinders.edu.au" target="_blank">Flinders University</a>, Adelaide, Australia <br>
September 2021 <br>
<a href=mailto:corey.bradshaw@flinders.edu.au>e-mail</a> <br>

## INSTRUCTIONS

1. Run initial base population models (<code>baseModels.R</code>). More details and justification of the model components can be found in:
    - Bradshaw, CJA, CN Johnson, J Llewelyn, V Weisbecker, G Strona, F Saltré. 2021. <a href="http://doi.org/10.7554/eLife.63870">Relative demographic susceptibility does not explain the extinction chronology of Sahul’s megafauna</a>. <em>eLife</em> 10: e63870. doi:10.7554/eLife.63870
    - and its associated Github repository <a href="https://github.com/cjabradshaw/MegafaunaSusceptibility">here</a>

3. Run projection scenarios using outputs from Step 1:
    - <code>stableCatastrophe.R</code>: stable, fixed carrying capacity (<em>K</em>) with generationally scaled catastrophe
    - <code>pulseCatastrophe.R</code>: stable, fixed <em>K</em> with 90% mortality pulse disturbance at 20<em>G</em> (generations); generationally scaled catastrophe
    - <code>r001-01Catastrophe.R</code>: stable, fixed <em>K</em> with density-independent mortality to cause mean <em>r</em> = -0.001 or -0.01 (user choice); generationally scaled catastrophe
    - <code>KstochCatastrophe.R</code>: stable, mean stochastic <em>K</em> (constant variance); generationally scaled catastrophe
    - <code>KstochVarIncCatastrophe.R</code>: stable, mean stochastic <em>K</em> (variance doubles over projection interval); generationally scaled catastrophe
    - <code>KstochDeclCatastrophe.R</code>: declining (-0.001), mean stochastic <em>K</em> (constant variance); generationally scaled catastrophe
    - <code>NoDFstableCatastrophe.R</code>: stable, mean trajectory with no density feedback (no <em>K</em>)
4. Import resulting .RData files from Step 2 to estimate evidence and strength of ensemble density feedback (<code>logisticGrowthFit.R</code>
5. Import resulting .RData files from Step 2 to estimate stationarity metrics
    - <code>StatStableCatastrophe.R</code>: stable, fixed <em>K</em> with generationally scaled catastrophe
    - <code>StatPulseCatastrophe.R</code>: stable, fixed <em>K</em> with 90% mortality pulse disturbance at 20<em>G</em>; generationally scaled catastrophe
    - <code>Statr001Catastrophe.R</code>: stable, fixed <em>K</em> with density-independent mortality to cause mean <em>r</em> = -0.001; generationally scaled catastrophe
    - <code>Statr01Catastrophe.R</code>: stable, fixed <em>K</em> with density-independent mortality to cause mean <em>r</em> = -0.01; generationally scaled catastrophe
    - <code>StatKstochCatastrophe.R</code>: stable, mean stochastic <em>K</em> (constant variance); generationally scaled catastrophe
    - <code>StatKstochVarIncCatastrophe.R</code>: stable, mean stochastic <em>K</em> (variance doubles over projection interval); generationally scaled catastrophe
    - <code>StatKstochDeclCatastrophe.R</code>: declining, (-0.001) mean stochastic <em>K</em> (constant variance); generationally scaled catastrophe
    - <code>StatStableNoDFCatastrophe.R</code>: stable, mean trajectory with no density feedback (no <em>K</em>)

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
- <code>matrixOperators.R</code>




