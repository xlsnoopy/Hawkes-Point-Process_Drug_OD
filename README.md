# Hawkes Point Process

A research project that studies Hawkes point process, and how it estimates missing information using fused data.

Applicable to geospatial data, with some high dimensional data associated with one of the datasets.

First tested on synthetic data. Simulated dataset using four different settings of parameters. Fused, sampled and labelled, using the labelled points to estimate unlabelled ones. Parameters are well recovered, and unlabelled points are decently well assigned.

Then tested on real data: opioid drug overdose cases (not disclosed). 

# Cite the work
The work has been published [here](https://www.imstat.org/publications/aoas/aoas_15_1/aoas_15_1.pdf#page=15). If you find this helpful, please consider citing:
> @article{LiuHawkes2021,
> author = {Xueying Liu and Jeremy Carter and Brad Ray and George Mohler},
> title = {{Point process modeling of drug overdoses with heterogeneous and missing data}},
> volume = {15}, 
> journal = {The Annals of Applied Statistics}, 
> number = {1},
> publisher = {Institute of Mathematical Statistics},
> pages = {88 -- 101}, 
> year = {2021}, 
> doi = {10.1214/20-AOAS1384},
> URL = {https://doi.org/10.1214/20-AOAS1384} 
> }

