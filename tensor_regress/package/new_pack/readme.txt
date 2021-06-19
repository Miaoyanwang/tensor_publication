Package ``tensorregress" Update

Date: 5/20/2021 Last edition: 06/19/2021 (Jiaxin)

1. Update the main function ``tensor_regress": 

- Add an option to the initialization, @param initial. initial = "random" for random initialization, initial = "QR_tucker" for deterministic initialization using QR-based tucker initialization. 

- If dist == "normal", initial ==  "QR_tucker", then the function returns the results after initialization.

- Option @Nsim is changed to @niter

- Add warnings if rank for the input data is smaller than the input core_shape. Change the NA to 0 when updating W1,W2,W3 and continue the iteration. 

- Change glm_modify for normal case, in line 77 of bricks.R. If fit1 = speedlm() returns negative RSS, try fit1 = lm(). 

- Add orthogonalization of the factor matrices with initial == "QR_tucker".

2. Update the function ``sim_Data":

- Add an option @param ortho. ortho == TRUE generates the side information matrices with orthogonal columns; ortho == FALSE generates the side information matrices with gaussian entries. 

3. Update the function ``sele_rank":

- Add an option @param initial. Update the lines using ``tensor_regress" with new option initial.

The folder has not been complied to a package yet. Just changed the function file in tensorregress/R/tensor_regress.R. 