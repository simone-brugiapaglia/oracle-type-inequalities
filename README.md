# On oracle-type local recovery guarantees in compressed sensing

This repository complements the paper 

*On oracle-type local recovery guarantees in compressed sensing*. Information and Inference: A Journal of the IMA (Accepted, 2020) by Ben Adcock, Claire Boyer, and Simone Brugiapaglia

A preprint of the paper can be found at [https://arxiv.org/abs/1806.03789](https://arxiv.org/abs/1806.03789).

## Content

To generate all plots, run the script main_script.m 

This is the role of each individual script:

- plotting_probabilities_1D_2D generates the plots in Figure 2;
- Haar_Fourier_test_1D.m generates the plots in Figure 3;
- function_approximation_test_1.m generates the plots in Figures 6 and 7;
- function_approximation_test_1.m generates the plots in Figure 8 and the data in Table 1.


## Dependencies

To run this script, you need to add to your Matlab path:

1. Numerical tour of Gabriel Peyré [https://www.numerical-tours.com/matlab/](https://www.numerical-tours.com/matlab/)
2. export_fig [https://github.com/altmany/export_fig](https://github.com/altmany/export_fig)
3. SPGL1 [https://github.com/mpf/spgl1](https://github.com/mpf/spgl1)

## Disclaimer

The following experiments involve randomization. Hence, the plots obtained by running this script might slightly differ from those in the paper.
