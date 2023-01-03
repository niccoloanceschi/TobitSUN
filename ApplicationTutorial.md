Introduction
============

As described in the [`README.md`](https://github.com/niccoloanceschi/TobitSUN/blob/main/README.md) file, this tutorial contains general guidelines and code to **replicate the  simulation studies reported in Section 5** of the paper, which focus on **tobit regression**.
For ease of presentation, we hereby consider single combinations of number of covariates *p* $\in \mathbb{N}$ and censoring percentage *k* = $\lfloor 100 \cdot \kappa \rfloor \in$ { $0,1,\dots,100$ }.
Original Figures and Tables from the paper gather the outcomes for multiple values of *p* and *k*, and can be thus reproduced by running the code below for all the settings analyzed therein.

Tobit simulations
=================

The code below shows how to **generate and load the data** corresponding to the different settings analyzed in Section 5 and explains in detail the `R` code to **implement the different methods for posterior inference** discussed in Section 4, including both sampling-based schemes and fast deterministic approximations. 

In the empirical studies from Section 5 we simulate a total of $n$ $\in \mathbb{N}$ observations from a standard tobit model under censoring proportions $\kappa \in [0,1]$, dividing the data in $n_0=\kappa\cdot n$ censored observations and $n_1=(1-\kappa)\cdot n$ uncensored ones.
Exploiting the latent utility interpretation of tobit regression, the responses $y_i$, $i=1, \ldots, n$, are obtained by first simulating the associated  utilities $z_i$, $i=1, \ldots, n$, from a univariate normal distribution $\cal{N}$(**x**<sub>i</sub>$^{\intercal}$**β**, 1), and then setting $y_i=z_i \mathbb{1} (z_i>z_T)$, for each $i=1, \ldots, n$ where $z_T$ is a pre-specified truncation threshold to obtain the desired proportion of censored observations under different values of $\kappa$.
The $p$ unit-specific predictors in **x**<sub>i</sub>, $i=1, \ldots, n$, are instead simulated from standard Gaussians, except for the intercept term, whereas the regression coefficients in **β** are generated from a uniform distribution in the range $[-5,5]$.
The response vector **y** $=  (y_1,\dots,y_n)^\intercal$ and design matrix **X** $=$ $($**x**$_1,\dots,$**x**$_n)^\intercal$ are divided into two sets **y**<sub>0</sub>,**X**<sub>0</sub> and **y**<sub>1</sub>,**X**<sub>1</sub>, associated respectively with censored and uncensored observations.
Other than the $n$ in-sample observations, we generate analogous $n$<sub>Test</sub> out-sample statistical units, later used for assessing predictive performances.
In doing so, we abide to the recommended practice of standardizing the predictors to have $0$ mean and standard deviation $0.5$ (see [Gelman et al., 2008](https://doi.org/10.1214/08-AOAS191) and [Chopin and Ridgway, 2017](https://doi.org/10.1214/16-STS581)).
Posterior inference relies on  spherical Gaussian priors $\cal{N}_p$(**0**$,\omega_p^2$**Ι**$_p$), with $\omega_p^2=25 \cdot 10 \, /p$, inducing increasing shrinkage in high dimensions.
Note that the varying threshold $z_T$ poses no difficulties in Bayesian inference since it directly enters the intercept term.

The results from Section 5 correspond to different combinations of the aforementioned settings, in particular:

-   *p* $\in$ { $10,20,50,100,200,400,800,1200$ }
-   *k* $\in$ { $15,50,85$ }
-   $n = 200$
-   $n$<sub>Test</sub> $= 200$

**For illustrative purposes, we focus here on the setting *p* $=$ 50 and *k* $=$ 50**. The other scenarios analyzed in Section 5 can be readily implemented by changing *p* and *k* in the code below.

Setting the hyperparameters, loading data and required packages
===============================================================

We recommend beginning by cleaning the workspace and setting the working directory to the folder where the source code is located.
Then, we load the dedicated source code [`functionsTobit.R`](https://github.com/niccoloanceschi/TobitSUN/blob/main/functionsTobit.R) and other required libraries, assumed to have been previously downloaded.

``` r
rm(list = ls())

#TODO: set custom working directory
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library("TruncatedNormal")
library("truncnorm")
library("emulator")
library("rlang")
library("rstan")
```

``` r
source("functionsTobit.R")
```

Next, we set the desired values of *p* and *k*, and the other hyperparameters.

``` r
# number of covariates
p = 50         
# percentange of censored observations
k = 50            

n = 200             # total in-sample observations
nTest = 200         # total out-sample observations

om2p = 25*10/p      # prior variances

seed = 123          # seed for random number generation 
    
nBurnIn = 1e3       # number of burn-in samples
nSample = 5e3       # number of MCMC samples

tol = 1e-3          # tolerance for stopping iterative updates
```

Then, we conclude the preliminary operations by loading (or generating) corresponding simulated data, and dividing them between censored and uncensored units.

``` r
input_dir = 'dataset'

if(!dir.exists(input_dir)){dir.create(input_dir)}

if(!file.exists(paste0(input_dir,"/simulatedData_n",n,"_p",p,"_k",k,".RData",sep=""))){
  generateData(n,p,k,nTest,seed=seed,input_dir=input_dir)
}
load(paste0(input_dir,"/simulatedData_n",n,"_p",p,"_k",k,".RData",sep=""))
    
# censored observations    
X0 = X[y==0,]
y0 = y[y==0]

# uncensored observations
X1 = X[y>0,]
y1 = y[y>0]

n1 = length(y1)
n0 = n-n1
```

Sampling-based posterior inference
==================================

We begin our analysis by considering sampling-based approaches to posterior inference, including both an i.i.d. sampler from the exact SUN posterior and Hamiltonian Monte Carlo, specifically the no-u-turn sampler by [Hoffman and Gelman (2014)](http://jmlr.org/papers/v15/hoffman14a.html).

Exact sampling from the unified skew-normal posterior (i.i.d.)
--------------------------------------------------------------

Samples from the exact unified skew-normal posterior are drawn via the algorithm by [Durante  (2019)](https://doi.org/10.1093/biomet/asz034), which relies on the method by [Botev (2017)](https://doi.org/10.1111/rssb.12162) for sampling from the multivariate truncated normal component.

``` r
startTime = Sys.time()
betaSUN = rSUNpost(X1,y1,X0,y0,om2p,zT,nSample,seed=seed)
timeSUN = difftime(Sys.time(), startTime, units=("secs"))[[1]]
```

We leverage the i.i.d. samples obtained above to evaluate empirical **posterior moments** and **predictive functionals**. Later on, we are going to use these quantities as a ground-truth to validate and compare the outcomes of the different approximation schemes. As for predictive functionals, we assess performance by computing the expected value and censoring probability for every statistical unit *y*<sub>NEW</sub>$,$**x**<sub>NEW</sub> within the $n$<sub>Test</sub> observations (refer to the article for the analytical expression of such quantities). Note that, as mentioned before, in evaluating  such quantities we implicitly absorb the censoring threshold *z*<sub>T</sub> within the intercept term.

``` r
# First 2 marginal moments
meanBetaMCMC = apply(betaSUN$sampleSUN,1,mean)
logDevStdBetaMCMC = log(apply(betaSUN$sampleSUN,1,sd))
    
# Predictives functionals 
locPredMCMC = Xtest%*%betaSUN$sampleSUN
predMeanMCMC = dnorm(locPredMCMC) + locPredMCMC * pnorm(locPredMCMC)
predProbMCMC = pnorm(-locPredMCMC)

# MC average
locPredMCMC = apply(locPredMCMC,1,mean)
predMeanMCMC = apply(predMeanMCMC,1,mean)
predProbMCMC = apply(predProbMCMC,1,mean)
```

Hamiltonian no u-turn sampler (NUTS)
------------------------------------

We now turn our attention to the Hamiltonian no u-turn sampler by [Hoffman and Gelman (2014)](http://jmlr.org/papers/v15/hoffman14a.html), resorting in particular to its [`rstan`](https://mc-stan.org/) implementation.

As a first step, `rstan` requires defining a probabilistic model compatible with the structure of tobit regression.
To do so, we recall that the tobit regression can be mapped into a totally unbalanced probit regression, where the Gaussian prior has been updated by taking into account the uncensored data, while the binary response vector is made only of zeros corresponding to censored observations.

Accordingly, we define below the familiar structure of probit regression under a general Gaussian prior, further representing the parameters of interest **β** as a suitable transformation of an underlying vector **β**<sub>0</sub> with standard normal entries.  

``` r
probitModel = " data{
                    int<lower=0> K;
                    int<lower=0> N;
                    int<lower=0,upper=1> Y[N];
                    matrix[N,K] X;
                    matrix[K,K] Chol;
                    vector[K] xi;
                }
                parameters {
                    vector[K] beta0;
                }
                transformed parameters {
                    vector[K] beta;
                    beta = xi + Chol * beta0;
                }
                model {
                    for(j in 1:K)
                        target += normal_lpdf(beta0[j] | 0, 1);
                    for(n in 1:N)
                        target += normal_lcdf((2*Y[n]-1)*(X[n]*beta) | 0, 1);
                } "
```

Secondly, we evaluate the parameters of the updated Gaussian prior $\cal{N}_p$(**ξ**<sub>1</sub>,**Ω**<sub>1</sub>).
Furthermore, we select a good initialization point for the sampler, by exploiting the results of an approximation of the posterior (EP in this case), with the goal of accelerating convergence of the MCMC chain.

``` r
startTime = Sys.time()
if(p<n1){
  Omega_post = solve(diag(1./om2p,p,p) + crossprod(X1))
}else{
  Omega_post = diag(om2p,p,p) - 
    om2p*crossprod(X1,solve(diag(1.,n1,n1) + om2p*tcrossprod(X1))%*%X1)*om2p
}
Omega_post = 0.5*(Omega_post+t(Omega_post))
Chol_post = t(chol(Omega_post))
xi_post = Omega_post %*% (crossprod(X1,as.matrix(y1)) + c(-zT/om2p,rep(0,p-1)))

data_stan = list(N=n0,K=p,Y=as.vector(y0),X=X0,Chol=Chol_post,xi=as.vector(xi_post))

paramsEP = getParamsEP(X1,y1,X0,y0,om2p,zT)
beta0_init = as.vector(solve(Chol_post)%*%(paramsEP$meanBeta-xi_post))

timeHMC = difftime(Sys.time(), startTime, units=("secs"))[[1]]
```

We can now implement the sampler itself, and compare its running time with that of the i.i.d. sampler.

``` r
startTime = Sys.time()
HMC_Samples = stan(model_code = probitModel, iter = nSample+nBurnIn, chains=1,
                        data = data_stan, warmup=nBurnIn, algorithm="NUTS", 
                        seed=seed, init=list(list(beta0=beta0_init)))
     
timeHMC = timeHMC + difftime(Sys.time(), startTime, units=("secs"))[[1]]
```

Running times
-------------

We can now visualize the running times, which correspond to producing a zoom of Table 1 from Section 5 of the paper, for the selected values of *p* and *k* (or equivalently $\kappa$).

``` r
table_sampling <- matrix(c(timeHMC,timeSUN),2,1)
rownames(table_sampling) <- c("NUTS","i.i.d")
colnames(table_sampling) <- c("Running Times [sec]")
table_caption = paste("k = ",k," & p = ",p,sep="")
knitr::kable(table_sampling,caption = table_caption, align = "c", digits=2)
```

Fast deterministic approximations
=================================

We now turn our attention to the deterministic approximation procedures described in Section 4.3 of the paper, including in particular mean-field VB, partially-factorized variational inference and expectation propagation.

Mean-field variational Bayes (MF-VB)
------------------------------------

We start by implementing mean-field variational Bayes to obtain the optimal Gaussian density *q*\*<sub>MF-VB</sub>(**β**) within the factorized family of joint densities
$\cal{Q}$<sub>MF-VB</sub> = { *q*(**β**,**z**<sub>0</sub>) :  *q*(**β**,**z**<sub>0</sub>) = *q*(**β**) *q*(**z**<sub>0</sub>) }, where **z**<sub>0</sub> represents the latent Gaussian utilities associated with censored observations **y**<sub>0</sub>.

``` r
# Get parameters
paramsMF = getParamsMF(X1,y1,X0,y0,om2p,zT)
```

As mentioned above, we validate empirically the quality of the different approximation schemes by comparing the associated approximate posterior moments and predictive functionals with that obtained via i.i.d. sampling.
As for predictive functionals of interest, the Gaussianity of MF approximation leads to simple closed-form expressions easily evaluated.

``` r
# First 2 marginal moments
nIterMF = paramsMF$nIter
meanBetaMF = paramsMF$meanBeta
logDevStdBetaMF = log(sqrt(paramsMF$diagV))

# Predictive functionals
if(p>n){
  invXXt = solve(diag(1.,n,n)+om2p*tcrossprod(X))
  XXtT = tcrossprod(X,Xtest)
  sdPredMF = as.vector(sqrt(1.+om2p*(rowSums(Xtest^2)-om2p*quad.diag(invXXt,XXtT))))
} else {
  V = solve(diag(1.,p,p)/om2p + crossprod(X))
  sdPredMF = as.vector(sqrt(1+quad.diag(V,t(Xtest))))
}

locPredMF = as.vector(Xtest%*%meanBetaMF)
predMeanMF = locPredMF*pnorm(locPredMF/sdPredMF) + sdPredMF*dnorm(locPredMF/sdPredMF)
predProbMF = pnorm(-locPredMF/sdPredMF)
```

Partially-factorized mean-field variational Bayes (PFM-VB)
----------------------------------------------------------

Secondly, we implement partially-factorized mean-field variational Bayes, returning the optimal approximate joint posterior *q*\*<sub>PFM-VB</sub>(**β**,**z**<sub>0</sub>) within the class $\cal{Q}$<sub>PFM-VB</sub> = {*q*(**β**,**z**<sub>0</sub>) :  *q*(**β**,**z**<sub>0</sub>) = *q*(**β** |  **z**<sub>0</sub>) $\prod~_{i=1}^{n_0}$ *q*(*z*<sub>0[i]</sub>) }.

As detailed in Section 4.3.1 of the paper, the optimal density *q*\*<sub>PFM-VB</sub>(**β**) = $\int$ *q*\*<sub>PFM-VB</sub>(**β**, **z**<sub>0</sub>) d**z**<sub>0</sub> = $\int$ *p*(**β** | **z**<sub>0</sub>) $\prod~_{i=1}^{n_0}$ *q*\*<sub>PFM-VB</sub>(*z*<sub>0[i]</sub>) d*z*<sub>0[i]</sub> corresponds to a SUN random variable with a specific structure.
Indeed, the covariance matrix appearing within the Gaussian cdf term of the SUN density *q*\*<sub>PFM-VB</sub>(**β**) has non-zero elements only on the main diagonal, which allow for straightforward calculations of the corresponding moments.
We refer to [Fasano, Durante and Zanella (2022)](https://doi.org/10.1093/biomet/asac026) for further explanations.

``` r
# Get parameters
paramsPFM = getParamsPFM(X1,y1,X0,y0,om2p,zT,predictive=TRUE)
```

Conversely, predictive functionals of interest can be evaluated by sampling from the resulting approximate posterior, still benefitting from the simplified structure of the latter.

``` r
# First 2 marginal moments
nIterPFM = paramsPFM$nIter
meanBetaPFM = paramsPFM$postMoments.meanBeta
logDevStdBetaPFM = log(sqrt(paramsPFM$postMoments.varBeta))

# Predictive functionals
if(p>n){
  invXXt = solve(diag(1.,n,n)+om2p*tcrossprod(X))
  XXtT = tcrossprod(X,Xtest)
  sdPredPFM = as.vector(sqrt(1.+om2p*(rowSums(Xtest^2)-om2p*quad.diag(invXXt,XXtT))))
} else {
  V = solve(diag(1.,p,p)/om2p + crossprod(X))
  sdPredPFM = as.vector(sqrt(1+quad.diag(V,t(Xtest))))
}

set.seed(seed)
muTN = paramsPFM$mu
muTN[y0==0] = -muTN[y0==0]
sampleTruncNorm = matrix(rtruncnorm(n0*nSample, a = 0, b = Inf, mean = muTN,
    sd = sqrt(paramsPFM$sigma2)), nrow = n0, ncol = nSample, byrow = F )
sampleTruncNorm[y0==0,] = -sampleTruncNorm[y0==0,]
locationVector = paramsPFM$postPredictive.VXt%*%sampleTruncNorm + 
  matrix(rep(paramsPFM$postPredictive.VinvOmega0beta0,nSample), ncol = nSample)

locPredPFM = Xtest%*%locationVector
predMeanPFM = locPredPFM*pnorm(locPredPFM/sdPredPFM) + sdPredPFM*dnorm(locPredPFM/sdPredPFM)
predProbPFM = pnorm(-locPredPFM/sdPredPFM)

# MC average
locPredPFM = apply(locPredPFM,1,mean)
predMeanPFM = apply(predMeanPFM,1,mean)
predProbPFM = apply(predProbPFM,1,mean)
```

Expectation-propagation (EP)
----------------------------

We conclude by producing the code for obtaining the optimal EP approximation *q*\*<sub>EP</sub>(**β**), as described in Section 4.3.2 of the paper. As highlighted therein, we introduce a novel and efficient implementation of the EP updates, achieving linear cost per iteration in the number of covariates *p*. 
Without such a feature, the quadratic cost in *p* of previous EP implementations would provide prohibitively slow routines in high-dimensional scenarios, compared to the variational procedures implemented above.

``` r
# Get parameters
paramsEP = getParamsEP(X1,y1,X0,y0,om2p,zT,predictive=TRUE)
```

As before, we evaluate marginal moments and predictive functionals, benefitting once more from the Gaussianity of the approximation.

``` r
# First 2 moments
nIterEP = paramsEP$nIter
meanBetaEP = paramsEP$meanBeta
logDevStdBetaEP = log(sqrt(paramsEP$diagOmega))

# Predictive functionals
if(p>=n0){
  XXtT = X1%*%t(Xtest)
  invXXt = solve(diag(1.,n1,n1)+om2p*tcrossprod(X1))
  XtestU = Xtest%*%paramsEP$postPredictive.U
  XtestU0 = Xtest%*%paramsEP$postPredictive.U0
  sdPredEP = as.vector(sqrt(1.+om2p*(rowSums(Xtest^2) -
                om2p*quad.diag(invXXt,tcrossprod(X1,Xtest))) -
                rowSums(XtestU*t(paramsEP$kEP*t(XtestU0))) ))
} else {
  OmegaEP = solve(diag(1.,p,p)/om2p + crossprod(X1) + crossprod(X0,paramsEP$kEP*X0))
  sdPredEP = as.vector(sqrt(1+quad.diag(OmegaEP,t(Xtest))))
}

locPredEP = as.vector(Xtest%*%meanBetaEP)
predMeanEP = locPredEP*pnorm(locPredEP/sdPredEP) + sdPredEP*dnorm(locPredEP/sdPredEP)
predProbEP = pnorm(-locPredEP/sdPredEP)
```

Running times and number of iterations
--------------------------------------

We can now visualize the running times, as medians over 20 independent repetitions, which correspond to producing a zoom of Table 2 from Section 5 of the paper, for the selected values of *p* and *k* (or equivalently $\kappa$).

``` r
library("microbenchmark")

run_approx <- summary(microbenchmark(getParamsMF(X1,y1,X0,y0,om2p,zT),
                                     getParamsPFM(X1,y1,X0,y0,om2p,zT),
                                     getParamsEP(X1,y1,X0,y0,om2p,zT),
                                     unit = "s", times = 20))

table_approx_time <- matrix(c(run_approx$median,nIterMF,nIterPFM,nIterEP),3,2)
rownames(table_approx_time) <- c("MF-VB","PFM-VB","EP")
colnames(table_approx_time) <- c("Running Times [sec]","Number of Iterations")
table_caption = paste("k = ",k," & p = ",p,sep="")
knitr::kable(table_approx_time,caption = table_caption, align = "cc", digits=c(3,0))
```

Assessing quality of different approximations
---------------------------------------------

Finally, we can visualize the results of different approximation schemes, comparing their outcomes with the ground-truth provided by the i.i.d. sampler from the exact posterior.
Since we are hereby concentrating on a single combination of $\kappa$ and $p$, we cannot directly reproduce Figure 1 of the paper, which nonetheless arises by gathering analogous results for different values of $\kappa$ and $p$.
Instead, we report here scatter plots for the same posterior functionals analyzed therein, unraveling the behavior behind the median deviations reported in the corresponding points of Figure 1.

We begin by loading the required library to generate the desired data visualization.

``` r
library(latex2exp)
library(reshape2)
library(ggplot2)
library(plyr)
```

Then, we structure the data to construct a single plot.

``` r
# Method's column
method_list = c("MF-VB","PFM-VB","EP")
methodMom = matrix(mapply(rep,method_list,MoreArgs=list(times=p)),ncol=1)
methodPred = matrix(mapply(rep,method_list,MoreArgs=list(times=nTest)),ncol=1)

# Posterior Mean
meanData = data.frame("Method"=methodMom,
                     "value"=c(meanBetaMF,meanBetaPFM,meanBetaEP),
                     "groundtruth"=rep(meanBetaMCMC+c(zT,rep(0,p-1)),length(method_list)))
meanData$value = meanData$value + c(zT,rep(0,p-1))
meanData$obj = "mean"

# Posterior Standard Deviations
varData = data.frame("Method"=methodMom,
                     "value"=c(exp(logDevStdBetaMF),exp(logDevStdBetaPFM),exp(logDevStdBetaEP)),
                     "groundtruth"=rep(exp(logDevStdBetaMCMC),length(method_list)))
varData$obj = "stdDev"

# Predictive Probabilities
predData = data.frame("Method"=methodPred,
                     "value"=c(predProbMF,predProbPFM,predProbEP),
                     "groundtruth"=rep(predProbMCMC,length(method_list)))
predData$obj = "predProb"

# Predictive Means
predmeanData = data.frame("Method"=methodPred,
                     "value"=c(predMeanMF,predMeanPFM,predMeanEP),
                     "groundtruth"=rep(predMeanMCMC,length(method_list)))
predmeanData$obj = "predMean"

# Combining the data
dataPoints = rbind(meanData,varData,predData,predmeanData)
dataPoints$obj = factor(dataPoints$obj, levels=c("mean", "stdDev","predProb", "predMean"))
dataPoints$obj = mapvalues(dataPoints$obj,
        from = c("mean", "stdDev", "predProb", "predMean"),
        to = c(TeX("$E \\[ \\beta_j | y \\]$"),
        TeX("$Var \\[ \\beta_j | y \\,\\]^{1/2}$"),
        TeX("$P \\[ y_{new} = 0 | y \\]$"),
        TeX("$E \\[ y_{new} | y \\]$")))
```

Finally, we produce the Figure with the desired scatter plots.

``` r
myColors <- c("EP" = "#FFB000", "MF-VB" = "#DC267F", "PFM-VB"="#785EF0")
myShapes <- c("EP" = 0, "MF-VB"=1, "PFM-VB"=2)

ggplot(dataPoints, aes(y=value, x=groundtruth,color=Method)) +
  geom_point(aes(shape=Method),alpha=0.7,size=2) + 
  facet_wrap(obj~., scales="free", ncol=2,labeller=label_parsed) + 
  geom_abline(slope = 1, intercept = 0, alpha = 0.7, lty = 2) +
  labs(x="Quantities computed from the exact posterior",
       y = "Quantities computed from the approximate posterior") + 
  theme_bw() + theme(axis.title.x = element_text(size=9), axis.title.y = element_text(size=9), 
        plot.margin = margin(0.1, 0.1, 0.05, 0.2, "cm"), legend.position = "bottom",
        aspect.ratio = 0.7) +  
  scale_color_manual(values=myColors) + scale_shape_manual(values=myShapes)
```







