Boumis, G., Geist, E. L., & Lee, D. (2023). Bayesian hierarchical modeling for probabilistic estimation of tsunami amplitude from far-field earthquake sources. Journal of Geophysical Research: Oceans, 128, e2023JC02000. [https://doi.org/10.1029/2023JC020002](https://doi.org/10.1029/2023JC020002)
<p align="center">
  <img src="map.png"/>
</p>

## Summary
Assessing tsunami hazard for a coastal area, e.g., nearshore tsunami height, typically requires us to resort to physics-based models since historical data from individual monitoring stations are scarce. In this work, however, we develop a purely statistical model for tsunami data from far-field earthquake sources observed at tide gauges along the shorelines of California and Oregon. Our spatial Bayesian hierarchical model can artificially augment the data catalog of each tide gauge by capturing spatial dependence between stations, and thus allows for a more robust tsunami hazard analysis than analyzing empirical station data in isolation.

## Details
The Stan code for the Bayesian hierarchical model along with an R script to run the model are included in the ***"code"*** folder. Different data used for this analysis are provided inside the ***"data"*** folder. Additional scripts for other analyses described in the article are available upon request.
