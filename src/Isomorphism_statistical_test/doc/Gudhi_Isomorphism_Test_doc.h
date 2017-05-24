/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Claire Brécheteau
 *
 *    Copyright (C) 2017  INRIA Sophia-Saclay (France)
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


#ifndef DOC_GUDHI_ISOMORPHISM_TEST_H_
#define DOC_GUDHI_ISOMORPHISM_TEST_H_

namespace Gudhi {

namespace isomorphism_test {

/**  \defgroup isomorphism_test Isomorphism Test
 *
 * \author   Claire Brécheteau
 *
 * @{
 *

 * This module aims at comparing two samples of points from two general metric spaces. (Recall that a metric, or distance \f$d\f$ on a space \f$\mathcal{X}\f$ is 
 * a non-negative map defined on \f$\mathcal{X}\times\mathcal{X}\f$, symetrical \f$\forall x,y\in\mathcal{X}, d(x,y)=d(y,x)\f$, satisfying the triangle inequality 
 * \f$\forall x,y,z\in\mathcal{X}, d(x,z)\leq d(x,y)+d(y,z)\f$ and such that \f$\forall x,y\in\mathcal{X}, (d(x,y)=0)\Rightarrow (x=y)\f$).
 * Thus, the input should be, either a sample \f$P=\{X_1,X_2,\ldots,X_N\}\f$, a distance \f$d\f$ defined on all pairs \f$(X_i,X_j)\f$, a sample 
 * \f$Q=\{{X'}_1,{X'}_2,\ldots,{X'}_{N'}\}\f$ and a distance \f$d'\f$ defined on all pairs \f$({X'}_i,{X'}_j)\f$, either two distance matrices
 * (a distance matrix \f$M=(m_{i,j})_{i,j}\f$ satisfies : \f$m_{i,j}=d(X_i,X_j)\f$ for some sample \f$P=\{X_1,X_2,\ldots,X_N\}\f$ and some distance \f$d\f$.).

 * Note that we can use as an input a matrix  \f$M=(m_{i,j})_{i,j}\f$ for which \f$m_{i,j}=d(X_i,X_j)\f$ if the point \f$X_j\f$ is one of the \f$k\f$ nearest neighbours of \f$X_i\f$, and \f$+\infty\f$ if not. Then the function test can be used for parameters m smaller than \f$\frac{k}{N}\f$, we can also plot signatures. Nonetheless, the function check_test_error cannot be used, since subsampling presumes the knowledge of all pairwise distances.

\section theory Some Theory

 * In a more statistical setting, the points \f$X_1,\ X_2,\ldots\ X_N,\ {X'}_1,\ {X'}_2,\ldots\ X'_{N'}\f$ are independent random variables, the \f$X_i\f$s are from 
 * a borel probability measure \f$\mu\f$ defined on a metric space \f$(\mathcal{X},d)\f$. We say that the \f$X_i\f$s are from the metric-measure space, or 
 * <em>mm-space</em> \f$(\mathcal X,d,\mu)\f$. As well, the \f$X'_i\f$s are independent random variables from the mm-space \f$(\mathcal Y,d',\nu)\f$.

 * The isomorphism statistical test proposed in this module aims at answering the following question, with low probability of mistake :
 * Are the mm-spaces \f$(\mathcal X,d,\mu)\f$ and \f$(\mathcal{Y},d',\nu)\f$ isomorphic ? That is, can we find a map 
 * \f$\phi:(\mathcal X,d,\mu)\rightarrow(\mathcal Y,d',\nu)\f$ which is one-to-one and onto, which is an isometry, that is 
 * \f$\forall\ x,x'\in\mathcal{X}, d'(\phi(x),\phi(x'))=d(x,x')\f$ and which preserves the measure, that is \f$\nu(\phi(B))=\mu(B)\f$ for all Borel set 
 * \f$B\f$ of \f$\mathcal X\f$.
 * We can set the hypotheses

 * \f$H_0\f$ : \f$(\mathcal{X},d,\mu)\f$ and \f$(\mathcal{Y},d',\nu)\f$ are isomorphic,

 * and the hypothesis 

 * \f$H_1\f$ : \f$(\mathcal X,d,\mu)\f$ and \f$(\mathcal Y,d',\nu)\f$ are not isomorphic".

 * The statistical test comes with theoretical guarantees, see the paper \cite brecheteau2017isomorphism_test. Although these guarantees are asymptotic, 
 * which means that if the numbers of points \f$N\f$ and \f$N'\f$ are not big enough, then the test may not "work". 
 * Fortunately, a function is proposed in the module to try to determine whether the test has a chance to work or not, based on a heuristic.

 * What do we mean by "the test works" ? We fix a parameter \f$\alpha\f$, that we may choose equal to 0.05 for instance, and we say that the test works
 * if the <em>error of type I</em> is smaller than \f$\alpha\f$, or at least very close to \f$\alpha\f$. By error of type I, we mean the probability that the
 * Hypothesis \f$H_1\f$ is retained when \f$H_0\f$ is actually true. The test is designed to fix this error of type I at a level smaller than \f$\alpha\f$, at least when 
 * \f$N\f$ and \f$N'\f$ go to \f$\infty\f$.

 * Nonetheless, the reader should be aware of the fact that fixing the error of type I does not allow to fix the <em>error of type II</em>, that is the probability
 * that the Hypothesis \f$H_0\f$ is retained when \f$H_1\f$ is actually true. Indeed, if the hypothesis \f$H_0\f$ is retained, it does not mean that \f$H_0\f$ is true
 * (that the mm-spaces are isomorphic), it actually means that there is no obvious evidence of non isomorphism.

 * The first reason for this impossibility to fix the error of type II is that we do not use a distance between mm-spaces up to an isomorphism, only pseudo-distances.
 * Indeed, it is possible that all the pseudo-distances we consider are zero for two mm-spaces that are not isomorphic. In this situation, the test will retain the
 * hypothesis \f$H_0\f$ with probability around \f$1-\alpha\f$, whereas \f$H_0\f$ is not true.
 * The second reason is that for general statistical tests, fixing the error of type I does not allow to fix the error of type II. But sometimes, we can derive 
 * upper-bounds for this error, this is linked to the notion of power. Some theoretical bounds depending on the pseudo-distances are available 
 * in \cite brecheteau2017isomorphism_test.

 * We also provide a function to help choosing the "best" pseudo-distance to use for the test, that is the pseudo-distance which offers a test of error of type II
 * as small as possible, while keeping the error of type I close to the parameter \f$\alpha\f$.


\section the_method The method used for the isomorphism statistical test
 * The statistical test of the hypothesis
 * \f$H_0\f$ : \f$(\mathcal{X},d,\mu)\f$ and \f$(\mathcal{Y},d',\nu)\f$ are isomorphic
 * versus
 * \f$H_1\f$ : \f$(\mathcal{X},d,\mu)\f$ and \f$(\mathcal{Y},d',\nu)\f$ are not isomorphic
 * is based on a pseudo-metric built as follows.
 * We associate to any mm-space \f$(\mathcal{X},d,\mu)\f$ a signature \f$S(\mu)\f$, that is, a Borel probability measure on \f$\mathbf{R}_+\f$ such that \f$S(\mu)=S(\nu)\f$ whenever \f$(\mathcal{X},d,\mu)\f$ and \f$(\mathcal{Y},d',\nu)\f$ are isomorphic.
 * \image html "images_doc/cap2.png" \f$\mu\f$ uniform measure on a cap shaped open subset of \f$\mathbf{R}^2\f$.
 * \image html "images_doc/hat2.png" \f$\nu\f$ uniform measure on a cap shaped open subset of \f$\mathbf{R}^2\f$.
 * \image html "images_doc/capsignature.png" cumulative distribution function of \f$S(\mu)\f$.
 * \image html "images_doc/hatsignature.png" cumulative distribution function of \f$S(\nu)\f$.
 *The pseudo-metric \f$D(\mu,\nu)\f$ considered is then equal to the \f$L_1\f$-norm between the cumulative distribution functions of the signatures, also knows as the \f$L_1\f$-Wasserstein distance:
 * \f[D(\mu,\nu)=W_1(S(\mu),S(\nu)).\f]
 * \image html "images_doc/Comparaisonhatcap.png" \f$W_1(S(\mu),S(\nu))\f$.

 * The true mm-spaces are unknown. Indeed, we do only have an access to two samples of points. Thus, the statistical test is based on the signatures associated to the samples :
 * \f[D(\mathbf{1}_P,\mathbf{1}_Q)=W_1(S(\mathbf{1}_P),S(\mathbf{1}_Q)),\f]
 * where \f$\mathbf{1}_P\f$ corresponds to the uniform measure on the sample \f$P\f$.

 * \image html "images_doc/cap2echant70pt.png" Sample P
 * \image html "images_doc/hat2echant70pt.png" Sample Q
 * \image html "images_doc/Cap70pt.png" cumulative distribution function of \f$S(\mathbf{1}_P)\f$.
 * \image html "images_doc/Hat70pt.png" cumulative distribution function of \f$S(\mathbf{1}_Q)\f$.

 * For each point \f$X_i\f$ in \f$P\f$, we use the notation
 * \f[d_{\mathbf{1}_P,m}(X_i)=\frac{1}{k}\sum_{j=1}^{k}d(X_i,X^{j}),\f]
 * where the \f$X^{j}\f$ are the \f$k\f$ nearest-neighbours of \f$X_i\f$ in \f$P\f$.

 * The family of signatures considered in this test depends on two parameters : m and n, and is defined as :
 * \f[S_{n,m}(P) = \frac{1}{n}\sum_{i=1}^{n}\delta_{d_{\mathbf{1}_P,m}(X_i)}.\f]
 * Here \f$\delta_x\f$ is the Dirac mass on the point \f$x\f$.
 * We also denote
 * \f[D_{n,m}(P,Q)=W_1(S_{n,m}(P),S_{n,m}(Q)).\f]
 * For m not too small and n not too big, but not too small neither, it is proved in \cite brecheteau2017isomorphism_test that we can approximate the distribution of \f$D_{n,m}(P,Q)\f$ under the hypothesis \f$H_0\f$, by subsampling.
 * Concretely, subsampling means repeating the following operation \f$\frac{N_{boot}}{2}\f$ times :
 * Compute \f$D^*_{n,m}(P)=W_1({S^*}_{n,m}(P),{{S^*}'}_{n,m}(P))\f$ and \f$D^*_{n,m}(Q)\f$,
 * where \f${S^*}_{n,m}(P)=\frac{1}{n}\sum_{i=1}^{n}\delta_{d_{\mathbf{1}_P,m}(X^*_i)}\f$ with the \f$X^*_i\f$ uniformly sampled from \f$\mathbf{1}_P\f$, the same for \f${{S^*}'}_{n,m}(P)\f$, independently.

 * The p-value is then equal to the proportion of elements \f$D^*_{n,m}(P)\f$ and \f$D^*_{n,m}(Q)\f$ bigger than the statistic observed \f$D_{n,m}(P,Q)\f$.

 * The hypothesis \f$H_0\f$ will be rejected whenever this p-value is smaller than \f$\alpha\f$.
 * Note that when \f$N\f$ and \f$n\f$ are large enough, it happens with probability \f$\alpha\f$ when \f$(\mathcal{X},d,\mu)\f$ and \f$(\mathcal{Y},d',\nu)\f$ are isomorphic, for some relevant choice of parameters m and n.

\section{Examples of usage}

 * This first example presents the method to get a p-value. Beware of the fact that the number of points in each sample in this example, \f$N=4\f$, is not big enough for the p-value to be meaningful. 
 * \include example_simple.cpp

 * This second example proposes an overview of all of the possible methods to read or compute a distance matrix for the test.
 * \include example_read_distance_matrix.cpp

 * This last example is a more complete version of the test since it first presents the method to choose the best parameter m, then plot the signatures associated to this m, helps choosing the biggest parameter n for which the approximate type I error is \f$\alpha\f$, and finally makes the test for the chosen parameters m and n, returning a p-value and the hypothesis retained.
 * \include example_step_by_step.cpp


\verbatim
{m,cost} = {0.05,0.000518096}.
{m,cost} = {0.1,0.000788196}.
{m,cost} = {0.2,0.00301515}.
{m,cost} = {0.3,0.00314847}.
{m,cost} = {0.4,0.00255544}.
{m,cost} = {0.5,0.00187678}.
{m,cost} = {0.6,0.00141873}.
{m,cost} = {0.7,0.00124559}.
{m,cost} = {0.8,0.000853334}.
Which m do you want to choose ?
0.3
{n,err type I} = {5,0.05}
{n,err type I} = {10,0.06}
{n,err type I} = {20,0.15}
{n,err type I} = {30,0.16}
Which n do you want to choose ?
10
The p-value is equal to : 0.49
The hypothesis retained is H0
\endverbatim



 * \section IsomorphismTestExamples Examples
 * End user programs are available in example/Isomorphism_test folder.
 * 
 * \copyright GNU General Public License v3.
 */
/** @} */  // end defgroup isomorphism_test

}  // namespace isomorphism_test

namespace Isomorphism_test = isomorphism_test;

}  // namespace Gudhi

#endif  // DOC_GUDHI_ISOMORPHISM_TEST_H_
