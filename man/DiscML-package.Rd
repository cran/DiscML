\name{DiscML-package}
\alias{DiscML-package}
\alias{DiscML}
\alias{print.DiscML}
\docType{package}
\title{
DiscML: An R package for estimating evolutionary rates
of discrete characters using maximum likelihood
}
\description{
DiscML was developed as a unified R program for estimating the evolutionary 
rates of discrete characters with no restriction on the number of character
states and having a great flexibility on transition models. DiscML performs 
maximum likelihood estimation with the options to correct for unobservable 
data, to implement a gamma distribution for rate variation, and to estimate 
the prior root probabilities from the empirical data. It gives users the 
ability to customize the instantaneous rate transition matrices, and to 
choose a variety of pre-determined matrices. DiscML is ideal for the 
analysis of binary (1s/0s) patterns, gene families, and multistate discrete 
morphological characteristics.
}

\details{
\tabular{ll}{
Package: \tab DiscML\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-04-28\cr
License: \tab GPL (>= 2) \cr
}

DiscML is flexible on both the size and type of the rate transition matrix.
The argument, 'model', specifies the rate transition matrix, which can be
customized by the user or chosen from pre-determined matrices. The 
pre-determined matrices in DiscML are:

1. 'ER': an equal-rate matrix, in which all non-diagonal entries are equal, 
e.g., matrix(c(0,1,1,1,0,1,1,1,0),ncol = 3, nrow = 3).

2. 'SYM': a symmetric matrix, which is identical with its transpose, e.g., 
matrix(c(0,1,2,1,0,3,2,3,0), ncol= 3, nrow = 3).

3. 'ARD': an all-rates-different matrix, in which all non-diagonal entries are
 free to vary, e.g., matrix(c(0,1,2,3,0,4,5,6,0), ncol = 3, nrow = 3).

4. 'GTR': a general time reversible matrix, which is the same as a combination
of arguments, model ="SYM", reversible = TRUE, and rootprobability = TRUE.

5. 'BDER': a birth-and-death matrix, in which all non-zero entries are equal, 
e.g., matrix( c(0,1,0,0,1,0,1,0,0,1,0,1,0,0,1,0), ncol= 4, nrow = 4)

6. 'BDSYM': a birth-and-death matrix with symmetric entries, e.g.,
 matrix( c(0,1,0,0,1,0,2,0,0,2,0,3,0,0,3,0), ncol= 4, nrow = 4).

7. 'BDARD': a birth-and-death matrix with all non-zero entries free to vary, 
e.g., matrix( c(0,1,0,0,2,0,3,0,0,4,0,5,0,0,6,0), ncol= 4, nrow = 4).

8. 'BDIER': a birth-death-and-innovation matrix with equal entries, ultimately
the same as 'BDER'.

8. 'BDISYM': a birth-death-and-innovation matrix with symmetric entries, e.g.,
 matrix( c(0,1,0,0,1,0,2,0,0,2,0,2,0,0,2,0), ncol= 4, nrow = 4).

9. 'BDIARD': a birth-death-and-innovation matrix with variable entries, e.g.,
 matrix( c(0,1,0,0,3,0,2,0,0,4,0,2,0,0,4,0), ncol= 4, nrow = 4).

10. 'BDBI' : a birth-and-death matrix with birth entries and death entries 
being equal respectively.  e.g.,
 matrix( c(0,1,0,0,2,0,1,0,0,2,0,1,0,0,2,0), ncol= 4, nrow = 4).

When the prior root probabilities are estimated, the argument 'reversible' 
allows reversibility for the pre-determined symmetric matrices, namely,
ER, SYM, BDER, BDIER, BDSYM, and BDISYM, by multiplying the entries by the 
corresponding root probabilities.

The user can customize entries in the rate matrix by assigning numerical
characters. The user can also assign constant numbers to particular 
entries. To do so, the user needs to set a matrix containing alphabetical 
characters, such that all numerical entries will be optimized, and all 
alphabetical characters or words will be assigned will assigned by constant 
values using argument 'const'.

For instance, model = matrix(c("0", "Amy", "2", "1", "0", "Jacob", "3", "Amy",
 "0"), ncol = 3, nrow =3), constant values can be assigned by,
const = list("Amy" = 0.4, "Jacob" = 0.9), 
const = data.frame("Amy" = 0.4, "Jacob" = 0.9), 
or const = matrix(c("Amy", 0.4, "Jacob", 0.9), nrow =2, ncol= 2). 

DiscML allows the optimization of the prior root probabilities by setting
rootprobability =TRUE. To avoid being trapped at a local maximum, 
DiscML initiates the optimization from multiple initial prior root 
probabilities. The user can customize a delta parameter, where 1/delta is the 
increment of different root probability, to specify the thoroughness of the 
optimization. The user can also specify a single set of initial prior root 
probabilities using, e.g., irootprobability = c(1/5,2/5,2/5).

DiscML considers rate variation among the character sites by implementing
 a discrete gamma distribution using alpha = TRUE.  
}

\arguments{
\item{x}{a vector, a matrix (in which each row is a vector for a gene family), 
a data frame, or a list.}

\item{phy}{phylogenetic information; an object of class "phylo" in "ape".}

\item{CI}{a logical specifying whether to return the 95\% confidence 
intervals the likelihood of the different states.}

\item{model}{a customized numeric matrix, or one of the pre-determined ones,
 {"ER", "ARD", "SYM", "BDER", "BDARD","BDSYM", "BDBI", "BDIER", "BDIARD",
 "BDISYM"}}

\item{const}{a matrix, a list, or a data frame}

\item{p}{a numeric vector that indicates the values according to the index of 
rate matrix. This is used when you want to use a constant rate matrix.}

\item{ip}{the initial values of the entries in the rate matrix used in the 
maximum likelihood estimation.}

\item{alpha}{a logical specifying whether to estimate the alpha parameter of 
a discrete gamma distribution in the maximum likelihood estimation, or a 
numeric value specifying an alpha value.}

\item{ialpha}{the initial alpha value in a gamma distribution used in the 
maximum likelihood estimation.}

\item{rootprobability}{a logical specifying whether to estimate the prior root 
probabilities in the maximum likelihood estimation.}

\item{irootprobability}{an initial numeric vector for the prior root
probabilities used in the maximum likelihood optimization, or randomly
generated initial vectors used in the optimization, by setting
irootprobability = "RANDOM".}

\item{ivnum}{a positive integer specifying the number of random initial
  vectors to be tested on for prior root probabilities, when
  irootprobability = "RANDOM".}

\item{zerocorrection}{a logical specifying whether to correct for unobservable 
data; see details in Felsenstein (1992).}

\item{mu}{a logical specifying whether to display the standardized rate estimate.}

\item{simplify}{a logical specifying whether to convert all nonzero character 
states to character state "1" and perform binary analysis.}

\item{delta}{a positive integer, such that 1/delta specifies the smallest increment 
an element of initial prior root probability can be increased for optimization.}

}

\author{
Tane Kim, Weilong Hao 

Maintainer: Weilong Hao <haow@wayne.edu>
}

\references{
Felsenstein, J. (1992). Phylogenies from restriction sites: A 
maximum-likelihood approach. Evolution, 46, 159--173.

Paradis, E., Claude, J., and Strimmer, K. (2004). APE: Analyses of 
Phylogenetics and Evolution in R language. Bioinformatics, 20, 289--290.

Yang, Z. (1994), Maximum likelihood phylogenetic estimation from DNA sequences
 with variable rates over sites: approximate methods. J Mol Evol, 39, 306--314.
}

\examples{
  
# The default arguments in DiscML are:
# model = "ER"
# mu = FALSE
# reversible = FALSE
# alpha = FALSE
# rootprobability = FALSE

x<- c(1,2,1,0,1)
phy <- rtree(length(x))
# x is a vector with 5 elements, and phy is a randomly generated
# 5-taxon tree (using rtree from the 'ape' package).
DiscML(x, phy)

x <- matrix(c(1,2,0,1,0,1),2,3)
# x here is a matrix
phy <- rtree(3) 
# phy is a randomly generated tree containing 3 tips.

DiscML(x, phy, model = "SYM")
# a symmetric rate transition matrix is used in the estimation.

DiscML(x, phy, rootprobability = TRUE)
# the prior root probabilities will be estimated.

DiscML(x, phy, rootprobability = c(1/16,5/16,10/16))
# the prior root probabilities are fixed to be 1/16, 5/16, and 10/16.

DiscML(x, phy, alpha = TRUE)
# the alpha value in a gamma distribution will be estimated.

DiscML(x, phy, alpha = 0.1)
# the alpha value in a gamma distribution is fixed to be 0.1.

DiscML(x, phy, rootprobability = TRUE, reversible = TRUE)
# DiscML allows the reversibility for the symmetric matrices, e.g., ER, SYM..

DiscML(x, phy, simplify = TRUE)
# DiscML can convert all non-zero character states to be '1's to perform
# binary analysis. 

DiscML(x, phy, mu =TRUE)
# During the likelihood calculation, the matrices are standardized,
# the rate parameter mu will be displayed.
}
