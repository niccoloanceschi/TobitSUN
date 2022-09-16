# Bayesian Conjugacy for Tobit Regression via Unified Skew-Normal Random Variables: Exact and Approximate Inference

This repository is associated with the article [Anceschi, Fasano, Durante and Zanella (2022+). *Bayesian conjugacy in probit, tobit, multinomial probit and extensions: A review and new results*](https://arxiv.org/abs/2206.08118). The **key contribution of the paper is outlined below**.

> In this article we review, unify and extend recent advances in Bayesian inference and computation for a core class of statistical models, crucially relying on partially or fully discretized latent Gaussian utilities. 
This includes pivotal statistical constructions for which conjugacy has been recelty derived, such as probit regression ([Durante, 2019](https://doi.org/10.1093/biomet/asz034)), multinomial probit ([Fasano and Durante, 2022](https://jmlr.org/papers/v23/20-735.html)), dynamic multivariate probit ([Fasano et al., 2021](https://doi.org/10.1007/s11222-021-10022-w)), Gaussian processes ([Cao et al., (2022)](https://doi.org/10.1080/10618600.2022.2036614)), skewed Gaussian processes ([Benavoli (2020)](https://doi.org/10.1007/s10994-020-05906-3)), and many others for which conjugacy results are still missing in the literature, above all tobit models ([Tobin, 1958](https://doi.org/10.2307/1907382)) and several exentesions of the above constructions to multivariate, skewed, non--linear and dynamic contexts.
We prove that the likelihoods induced by all these formulations share a common analytical structure that implies conjugacy with a broad class of distributions, namely the unified skew--normals ([SUN](https://doi.org/10.1017/CBO9781139248891)), that generalize Gaussians to skewed contexts.
We address avenues for improved posterior inference via novel closed–form expressions, i.i.d. samplers from the exact sun posteriors, and more accurate and scalable approximations from variational Bayes and expectation–propagation.
The performances of both exact and approximate routines are illustrated through extensive simulation studies in different scenarios, focusing in particular on tobit regression.

This repository provides **codes and tutorials to implement the inference methods discussed in the review**, in the specific case of **tobit regression**. 
In particular, the code from the present repository allows the replicate the results of the simulation studies, reported in Section 5 of the paper. 
The complete tutorial can be found in the file [`ApplicationTutorial.md`](https://github.com/niccoloanceschi/TobitSUN/blob/main/ApplicationTutorial.md) where we also provide details how to generate the synthetic data analyzed in the paper.

- The goal of the ***first part*** of the analysis is to compare the performance of an i.i.d. sampler from the exact posterior (`i.i.d.`), exploiting the properties of SUN random variables, with those of state-of-the-art competitor sampling schemes, in particular with the no-u-turn Hamiltonian sampler ([`NUTS`](http://jmlr.org/papers/v15/hoffman14a.html)).
- The ***second part*** focuses instead on the comparison of different approximate approaches to posterior inference, including mean-field variational Bayes ([`MF-VB`](https://doi.org/10.1080/01621459.2017.1285773)), partially-factorized variational Bayes ([`PFM-VB`](https://doi.org/10.1093/biomet/asac026)), and expectation-propagation ([`EP`](https://doi.org/10.1214/16-STS581)), using the results of the `i.i.d.` sampler as a ground-truth to validate empirically the accuracy of the different approximate methods.

We highlight that the `EP` routine reported in the present repository leverage on a novel efficient implementations proposed in the paper, crucially improving scalability to high-dimensional scenarios, where alternative implementations in the literature result non-competitively slow. </br>
We refer to Section 4 of the paper for further details on all the methods considered here.

The **functions to implement the above methods** can be found in the `R` source file [`functionsTobit.R`](https://github.com/niccoloanceschi/TobitSUN/blob/main/functionsTobit.R), and a **tutorial** explaining in detail the usage of these functions is available in the file [`functionsTutorial.md`](https://github.com/niccoloanceschi/TobitSUN/blob/main/functionsTutorial.md). 

All the analyses are performed with a **MacBook Pro (OS Big Sur, version 11.6.8, Processor 2,7 GHz Intel Core i7, RAM 16 GB)**, using an `R` version **4.1.0**.

**IMPORTANT**: Although a seed is set at the beginning of each routine, the outputs reported in [`ApplicationTutorial.md`](https://github.com/niccoloanceschi/TobitSUN/blob/main/ApplicationTutorial.md) may be subject to slight variations depending on which version of the `R` packages has been used in the code implementation. This is due to possible internal changes of certain functions when the package version has been updated. **However, the magnitude of these minor variations is negligible and does not affect the conclusions**.
