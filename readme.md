# Global LocoMotif Clustering

- newWrapper.py : used to call Fast Unfolding algorithm of Blondel et al. 
- renumber.py : used to renumber graphs
- newEvaluation.py : used to evaluate clusterings
- localmotifclustermain.cpp : used to run our extension of MAPPR.


Note that localmotifclustermain.cpp and newWrapper.py require the libraries Louvain modularity clustering (https://sites.google.com/site/findcommunities/)
and MAPPR (http://snap.stanford.edu/mappr/code.html).

Then place the newWrapper in the Community_latest folder (of Louvain modularity clustering), and localmotifclustermain.cpp in snap/examples/localmotifcluster/
