Package ``tensorregress" Update

Date: 5/20/2021 Last edition: 5/25/2021 (Jiaxin)

1. Update the main function ``tensor_regress": 
- Add an option to the initialization, @param initial. initial = "random" for random initialization, initial = "tucker" for deterministic initialization using naive tucker decomposition, and initial = "de_tucker" for QR-based initialization. 
- Add an option to the algorithm, @param alg. alg = "alter" for continuing the iterative algorithm until the converge, alg = "unsup" for returning the QR-based initialization results without the iterations, which is only suitable for normal data with initial = "de_tucker".


The folder has not been complied to a package yet. Just changed the function file in tensorregress/R/tensor_regress.R. The new main function is temporarily named as ``tensorregress1".
