# garcia_and_lohmueller_2020
Scripts for Garcia and Lohmueller 2020


# Computing Hrj

To compute Hrj and other measures of association based on entropy check out the C++ implementation use the script titled:

"measures_of_association_based_on_entropy.cpp"
 
Here we implement this computation in C and C++. Given a two-dimensional (2D) contingency table in the form of an integer array ‘nn[i][j]’, where i labels the x variable and ranges from 1 to ni, j labels the y variable and ranges from 1 to nj, this routine returns the entropy ‘h’ of the whole table, the entropy hx of the x distribution, the entropy hy of the y distribution, the entropy hygx of y given x, the entropy hxgy of x ‘given’ y, the dependency of uygx of y on x (equations above), the ‘dependency’ uxgy of x on y and the symmetrical dependency ‘uxy’. This code is heavily adapted off The Art of Scientific Computing Second Edition Page 633.

