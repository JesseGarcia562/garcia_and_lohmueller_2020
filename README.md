# garcia_and_lohmueller_2020
Scripts for Garcia and Lohmueller 2020


# Computing Hrj



Consider the game of **twenty questsions**, 'where by repeated yes/no questions you try to eliminate all except one correct possibility for an unknown object. Better yet, consider a generalization of the game, where you are allowed to ask multiple choice questions as well as binary (yes/no) ones'. The categories in these multilpe choice questions are **assumed** or **supposed** to be **mutually exclusive** and **exhaustive** ("as are 'yes' and 'no" <3 )



The value to you of an 'answer' increases with the number of possibilities that it eliminates. MOre specifically, an answer that eliminates all except a fraction **p** of the remaining possibilites can be assigned a value  ** - ln p** (**this is going to be a positive number because p < 1**). The purpose of the logarithm is to make the value **additive** since, for one reason, one question that eliminates all but 1/6 of the possibilities considered as good as two questions, that in sequence reduce the number by factors .5 and .3333. This **toy example** is important for our future forays into genomics, linkage disequilibrium, and haplotype structure. However, instead of focusing on biallelic Single Nucleotide Polymorphisms that might be physically linked due to the linear/circular structure of chromosomes on Earth, we will continue to study this example of **twenty questions**


In the above text we talk about the **value of an answer**. What is the **value of a question** ? If there are **I** possible answers to the question (i = 1, ..., I) and the fraction of possibilities consistent with the 'ith' are $I$ possible answers to the question $(i=1, \ldots, I)$ and the fraction of possibilities consistent with the $i$ th answer is $p_{i}$ (with the sum of the $p_{i}$ 's cqual to one), then the value of the question is the expectation value of the value of the answer, denoted $H$
\[
H=-\sum_{i=1}^{I} p_{i} \ln p_{i}
\]
In evaluating $(14.4 .6),$ note that
\[
\lim _{p \rightarrow 0} p \ln p=0
\]

The value $H$ lies between 0 and $\ln I .$ It is zero only when one of the $p_{i}$ 's is one, all the others zero: In this case, the question is valueless, since its answer is preordained. $H$ takes on its maximum value when all the $p_{i}$ 's are equal, in which case the question is sure to eliminate all but a fraction $1 / I$ of the remaining possibilities.

The value $H$ is conventionally termed the entropy of the distribution given by the $p_{i}$ 's, a terminology borrowed from statistical physics.

So far we have said nothing about the association of two variables; but suppose we are deciding what question to ask next in the game and have to choose between two candidates, or possibly want to ask both in one order or another. Suppose that one question, $x,$ has $I$ possible answers, labeled by $i,$ and that the other question,
$y,$ as $J$ possible answers, labeled by $j .$ Then the possible outcomes of asking both questions form a contingency table whose entries $N_{i j},$ when normalized by dividing by the total number of remaining possibilities $N,$ give all the information about the
$p^{\prime}$ s. In particular, we can make contact with the notation $(14.4 .1)$ by identifying
\[
\begin{array}{l}
p_{i j}=\frac{N_{i j}}{N} \\
\left.p_{i}=\frac{N_{i}}{N} \quad \text { (outcomes of question } x \text { alone }\right)
\end{array}
\]
$p_{\cdot j}=\frac{N_{\cdot j}}{N} \quad$ (outcomes of question $y$ alone)
The entropies of the questions $x$ and $y$ are, respectively,
\[
H(x)=-\sum_{i} p_{i} . \ln p_{i} . \quad H(y)=-\sum_{j} p_{\cdot j} \ln p_{\cdot j}
\]
The entropy of the two questions together is
\[
H(x, y)=-\sum_{i, j} p_{i j} \ln p_{i j}
\]
Now what is the entropy of the question $y$ given $x$ (that is, if $x$ is asked first)? It is the expectation value over the answers to $x$ of the entropy of the restricted $y$ distribution that lies in a single column of the contingency table (corresponding to the $x$ answer
\[
H(y \mid x)=-\sum_{i} p_{i} \cdot \sum_{j} \frac{p_{i j}}{p_{i}} \ln \frac{p_{i j}}{p_{i}}=-\sum_{i, j} p_{i j} \ln \frac{p_{i j}}{p_{i}}
\]
Correspondingly, the entropy of $x$ given $y$ is
\[
H(x \mid y)=-\sum_{j} p_{\cdot j} \sum_{i} \frac{p_{i j}}{p_{i j}} \ln \frac{p_{i j}}{p_{\cdot j}}=-\sum_{i, j} p_{i j} \ln \frac{p_{i j}}{p_{i j}}
\]


We can readily prove that the entropy of $y$ given $x$ is never more than the cntropy of $y$ alone, i.e., that asking $x$ first can only reduce the usefulness of asking
14. Statistical Description of Data
$634 \quad$ Chapter
$y$ (in which case the two variables are associated!):
\[
\begin{aligned}
H(y \mid x)-H(y) &=-\sum_{i, j} p_{i j} \ln \frac{p_{i j} / p_{i}}{p \cdot j} \\
&=\sum_{i, j} p_{i j} \ln \frac{p_{\cdot j} p_{i}}{p_{i j}} \\
& \leq \sum_{i, j} p_{i j}\left(\frac{p_{i j} p_{i}}{p_{i j}}-1\right) \\
&=\sum_{i, j} p_{i} \cdot p_{\cdot j}-\sum_{i, j} p_{i j} \\
&=1-1=0
\end{aligned}
\]
where the inequality follows from the fact
\[
\ln w \leq w-1
\]
We now have everything we need to define a measure of the "dependency" of $y$ on $x,$ that is to say a measure of association. This measure is sometimes called the uncertainty coefficient of $y .$ We will denote it as $U(y \mid x)$
\[
U(y \mid x) \equiv \frac{H(y)-H(y \mid x)}{H(y)}
\]
This measure lies between zero and one, with the value 0 indicating that $x$ and $y$ have no association, the value 1 indicating that knowledge of $x$ completely predicts
$y .$ For in-between values, $U(y \mid x)$ gives the fraction of $y$ 's entropy $H(y)$ that is lost if $x$ is already known (i.e., that is redundant with the information in $x$ ). In our game of "twenty questions," $U(y \mid x)$ is the fractional loss in the utility of question $y$ if question $x$ is to be asked first.

If we wish to view $x$ as the dependent variable, $y$ as the independent one, then interchanging $x$ and $y$ we can of course define the dependency of $x$ on $y$
\[
U(x \mid y) \equiv \frac{H(x)-H(x \mid y)}{H(x)}
\]



If we want to treat $x$ and $y$ symmetrically, then the useful combination turns out to be
\[
U(x, y) \equiv 2\left[\frac{H(y)+H(x)-H(x, y)}{H(x)+H(y)}\right]
\]
If the two variables are completely independent, then $H(x, y)=H(x)+H(y),$ so (14.4.17) vanishes. If the two variables are completely dependent, then $H(x)=$ $H(y)=H(x, y),$ so $(14.4 .16)$ equals unity. In fact, you can use the identities (easily proved from equations $14.4 .9-14.4 .12$ )
\[
H(x, y)=H(x)+H(y \mid x)=H(y)+H(x \mid y)
\]
to show that
\[
U(x, y)=\frac{H(x) U(x \mid y)+H(y) U(y \mid x)}{H(x)+H(y)}
\]
i.e., that the symmetrical measure is just a weighted average of the two asymmetrical measures $(14.4 .15)$ and $(14.4 .16),$ weighted by the entropy of each variable separately. Here is a program for computing all the quantities discussed, $H(x), H(y)$ $H(x \mid y), H(y \mid x), H(x, y), U(x \mid y), U(y \mid x),$ and $U(x, y)$
