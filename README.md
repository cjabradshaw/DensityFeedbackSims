# DensityFeedbackSims
Stochastic simulations of population abundance with known component density feedback on survival to test for ability to return ensemble feedback signal

R code accompanying paper:

BRADSHAW, CJA, S HERRANDO-PÉREZ. In review. Component and ensemble density feedbacks decoupled by density-independent processes. 

## INSTRUCTIONS:

1. Run initial base population models (<code>baseModels.R</code>
2. Run projection scenarios using outputs from Step 1:
    - <code>stableCatastrophe.R</code>: stable fixed K with generationally scaled catastrophe
    - <code>pulseCatastrophe.R</code>: stable fixed K with 90% mortality pulse disturbance at 20G; generationally scaled catastrophe
    - <code>r001-01Catastrophe.R</code>: stable fixed K with density-independent mortality to cause mean r = -0.001; generationally scaled catastrophe
    - <code>KstochCatastrophe.R</code>: stable fixed K with density-independent mortality to cause mean r = -0.01; generationally scaled catastrophe
    - <code>KstochVarIncCatastrophe.R</code>: stable mean stochastic K (constant variance); generationally scaled catastrophe
    - <code>KstochDeclCatastrophe.R</code>: stable mean stochastic K (variance doubles over projection interval); generationally scaled catastrophe
    - <code>NoDFstableCatastrophe.R</code>: declining (-0.001) mean stochastic K (constant variance); generationally scaled catastrophe
4. Import resulting .RData files from Step 2 to estimate evidence and strength of ensemble density feedback (<code>logisticGrowthFit.R</code>
5. Import resulting .RData files from Step 2 to estimate stationarity metrics



