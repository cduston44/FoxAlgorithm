# FoxAlgorithm
An Implementation of the Fox Algorithm for finding the fundamental groups of branched covering spaces over graphs.

This code is writen in Sage (https://www.sagemath.org/) and makes extensive use of GAP (https://www.gap-system.org/) functionality. The idea is to calculate the fundamental group of a branched coverging space of the 3-sphere, branched over a graph. The original motivation was to use the code to study the topological state space of Loop Qauntum Gravity (see https://arxiv.org/abs/1005.1057, https://arxiv.org/abs/1111.1252, https://arxiv.org/abs/1308.2934, https://arxiv.org/abs/2106.14188), but there is no explicit reference to LQG needed for any of this. Purely math fun!

Included here are:

FoxAlg-Wn.sage: Finds the fundamental group of all covers over the wheel graphs Wn for n=4, 5, 6, and 7.

FoxAlg-Bn.sage: Finds the fundamental group of all covers over the "bubble" graphs Bn for n=3 and 4.

FoxAlg-MC-ByGraph: Uses Monte Carlo to find the fundamental groups of a random set of covers over an order N wheel graph.

FoxAlgorithm-SingleGraph: This is a jupyter notebook that finds the fundamental group over a particular graph with a particular set of permutations labels, set by the user.

MWE: This is a minimal working example, demonstrating the memory / exception errors generated by some of...many of...the above codes (see below).

Although all the included scripts are guaranteed to find the fundamental group of a particular cover (in the form of a finitely presented group), they are not guaranteed to be able to make any further identification of those groups. There is a GAP routine for this (structure_description), which is included in various ways in most of the codes, but since it amounts to solving the word problem, it sometimes fails to find a structure. To accommidate this, the codes have built-in exceptions so that the whole thing doesn't crash while GAP struggles to identify a group, but these built-in exceptions also sometimes fail, leading to crashes. The MWE demonstrates this behavior, and by trial and error it can be observed that the codes fail when running over a few hundred covers. If anyone has any thoughts about how to fix this, please let me know!

A work around for this is to watch when they break, restart them from that point, and combine the data files at the end.

Further caution: I refer to "the number of generators" of a finitely presented group as "the order" of the group basically everywhere here - even if "rank" is closer to what I mean, I am still only talking about the number of generators!
