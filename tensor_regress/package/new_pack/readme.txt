Package ``tensorregress" Update

Date: 5/20/2021 Last edition: 5/20/2021 (Jiaxin)

1. Update the main function ``tensor_regress": 
	- Add an option to the initialization, @param initial. initial = "random" for random initialization, initial = "tucker" for deterministic initialization using tucker decomposition.
- Add an option to the algorithm, @param alg. alg = "alter" for the original alternative optimization, alg = "unsup" for the one step estimation using the unsupervised tucker decomposition, which is only suitable for normal data.


The folder has not been complied to a package yet. Just changed the function file in tensor regress/R. 
