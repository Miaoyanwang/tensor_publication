Package ``tensorregress" Update

Date: 5/20/2021 Last edition: 06/03/2021 (Jiaxin)

1. Update the main function ``tensor_regress": 

- Add an option to the initialization, @param initial. initial = "random" for random initialization, initial = "QR_tucker" for deterministic initialization using QR-based tucker initialization. 

- If dist == "normal", initial ==  "QR_tucker", then the function returns the results after initialization.

- Option @Nsim is changed to @niter

- Add warnings if updated W1, W2, W3 have a NA. Change the NA to 0 and continue the iteration. 

- With initial == "QR_tucker", add a warning if obtained core tensor G is not full rank. Note that numbers smaller than 1e-10 are set to 0. 

- Change glm_modify for normal case, in line 77 of bricks.R. If fit1 = speedlm() returns negative RSS, try fit1 = lm(). 

The folder has not been complied to a package yet. Just changed the function file in tensorregress/R/tensor_regress.R. 