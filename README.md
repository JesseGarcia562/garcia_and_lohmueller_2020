# garcia_and_lohmueller_2020
Scripts for Garcia and Lohmueller 2020


# Computing Hrj

## For slim

To compute Hrj for slim we used the R Based script titled:

"computing_hrj_slim.R"

To compute Hrj and other measures of association based on entropy check out the C++ implementation use the script titled:

"measures_of_association_based_on_entropy.cpp"
 
Here we implement this computation in C and C++. Given a two-dimensional (2D) contingency table in the form of an integer array ‘nn[i][j]’, where i labels the x variable and ranges from 1 to ni, j labels the y variable and ranges from 1 to nj, this routine returns the entropy ‘h’ of the whole table, the entropy hx of the x distribution, the entropy hy of the y distribution, the entropy hygx of y given x, the entropy hxgy of x ‘given’ y, the dependency of uygx of y on x (equations above), the ‘dependency’ uxgy of x on y and the symmetrical dependency ‘uxy’. This code is heavily adapted off The Art of Scientific Computing Second Edition Page 633.

# Computing Significance of Supplementary Figure 5 (KS Test)

> Effect of background selection on LD between S variants: this is interesting, and I think it merits further investigation, but Figure S5 does not convince me that the effect is real. Is the slight difference that is observed at d<500bp statistically significant?

This is a great question! Given these two sets of data we asked a similar question: "Can we disprove, to a certain required level of significance, the null hypothesis that two data sets are drawn from the same population distribution function?" Let's assume that disproving the null hypothesis in effect proves that the data sets are from different distributions. However, failing to disprove the null hypothesis only shows that the data sets can be expected under a single distribution function. According to the scientific method, a scientist can never prove that two data sets come from a single distribution since no practical amount of data can distinguish between two distributions wich differ only by one part in $10^{10}$ [1].

For continuous data as a function of a single variable the most generally accepted test is the Kolmogorov Smirnov test (KS test). The KS test is applicable to unbinned distributions that are functions of a single independent variable, that is, to data sets where each data point can be associated with a single number.This is important because we can covert this list of data points to an unbiased estimator. In such cases, the list of data points can be easily converted to an unbiased estimator $S_{N}(x)$ of the cumulative distribution function of the probability distribution from which it was drawn: If the $N$ events are located at values $x_{i}, i=1, \ldots, N,$ then $S_{N}(x)$ is the function giving the fraction of data points to the left of a given value $x .$ This function is obviously constant between consecutive (i.e., sorted into ascending order) $x_{i}$ 's, and jumps by the same constant $\left.1 / N \text { at each } x_{i} . \text   .\right)$

Different distribution functions, or sets of data, give different cumulative distribution function estimates by the above procedure. However, all cumulative distribution functions agree at the smallest allowable value of $x$ (where they are zero), and at the largest allowable value of $x$ (where they are unity). (The smallest and largest values might of course be $\pm \infty .$ So it is the behavior between the largest and smallest values that distinguishes distributions.

One can think of any number of statistics to measure the overall difference between two cumulative distribution functions: the absolute value of the area between them, for example. Or their integrated mean square difference. The KolmogorovSmirnov $D$ is a particularly simple measure: It is defined as the maximum value of the absolute difference between two cumulative distribution functions. Thus, for comparing one data set's $S_{N}(x)$ to a known cumulative distribution function $P(x),$ the $\mathrm{K}-\mathrm{S}$ statistic is
\[
D=\max _{-\infty<x<\infty}\left|S_{N}(x)-P(x)\right|
\]
while for comparing two different cumulative distribution functions $S_{N_{1}}(x)$ and $S_{N_{2}}(x),$ the $\mathrm{K}-\mathrm{S}$ statistic is
\[
D=\max _{-\infty<x<\infty}\left|S_{N_{1}}(x)-S_{N_{2}}(x)\right|
\]


One cool thing of the KS statistic is that its **distribution** in the case of the null hypothesis (which, lets remember, is 'data sets drawn from the same distribution') can be computed to an approximation. With this, we can get the significance that this reviewer is referring to given any observed nonzero value of this KS statistics (we call it D here). Another unique feature of the $\mathrm{K}-\mathrm{S}$ test is that it is invariant under reparametrization of $x$. For example, you will get the same significance using $x$ as using $\log x$.

The function that enters into the calculation of the significance can be written as the following sum:
\[
Q_{K S}(\lambda)=2 \sum_{j=1}^{\infty}(-1)^{j-1} e^{-2 j^{2} \lambda^{2}}
\]
which is a monotonic function with the limiting values
\[
Q_{K S}(0)=1 \quad Q_{K S}(\infty)=0
\]
In terms of this function, the significance level of an observed value of $D$ (as a disproof of the null hypothesis that the distributions are the same) is given approximately by the formula
\[
\text { Probability }(D>\text { observed })=Q_{K S}([\sqrt{N_{e}}+0.12+0.11 / \sqrt{N_{e}}] D)
\]

where $N_{e}$ is the effective number of data points, $N_{e}=N$ for the case of one distribution, and
\[
N_{e}=\frac{N_{1} N_{2}}{N_{1}+N_{2}}
\]
for the case of two distributions, where $N_{1}$ is the number of data points in the first distribution, $N_{2}$ the number in the second.

The nature of the approximation involved in is that it becomes asymptotically accurate as the $N_{e}$ becomes large, but is already quite good for $N_{e} \geq 4,$ as small a number as one might ever actually use [1]. von Mises discusses this in more detail but has a great review on current theory of the KS test.

For the case of two distributions we implemented two C++ functions kstwo() and probks() as defined in Press et. al. In base R this function is called ks.test() and in python3 this function can be called with scipy.stats.kstest()


Here we use base R's implementation and include a C++ implementation for readers with inexpensive GPU's who would like to replicate this statistical test: 

> https://github.com/JesseGarcia562/garcia_and_lohmueller_2020/blob/master/ks_test.cpp

## References
1. von Mises, R. 1964, Mathematical Theory of Probability and Statistics (New York: Academic Press), Chapters IX(C) and IX(E)

