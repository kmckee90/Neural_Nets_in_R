# Neural_Nets_in_R
Some neural nets in R, mainly for didactic purposes. 
Just a bit of fun and education.

The scripts were written in order of Boltzman, Helmholtz, then RecurrentHelmholtz.

The Boltzmann machine can learn discrete binary network, but with multiple layers the network is unlikely to reflect the data-generating structure.

Helmholtz demonstrates the recovery of data-generating node structures over multiple hidden layers.
The Helmholz script is the most refined and is successful at recovering very large data-generating neural networks.

The Recurrent Helmholz script does not use simulation, but instead learns from part of the wikipedia article on the American Civil War, with words mapped to binary sequences.
It definitely learns, but the best speech product I have seen so far was not convincing. 
The method is crude however, and a faster language would help allow it to learn more.
