# Inversion of a Dirac Matrix

In lattice calculations, often the following problem appears: A ensemble of configurations of gauge field \( U_\mu(x) \) have been generated. To calculate an observable of interest, a certain matrix, which depends on the gauge fields, needs to be inverted and applied to a given vector. Here, we will look at the following simplified setting:

Consider a 2-dimensional square lattice of size \( N_1 \times N_2 \) as a discretization of a 1+1 dimensional space-time. A fermion field is discretized in such a way that the field values are specified at each lattice point. Furthermore, an abelian gauge field is introduced at the links between neighboring lattice sites.

This abelian gauge field consists of complex numbers with absolute square 1. That is, the link between the lattice site \( x \) and the neighboring site in the direction \( \mu \), \( x + \hat{\mu} \), is called \( U_\mu(x) \) and has the form \( U_\mu(x) = \exp(ia(x)) \) where \( a(x) \) is a real number.


