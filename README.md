# Distributed Optimization of Average Consensus Containment with Multiple Stationary Leaders
This repository is the official implementation of the paper [Distributed Optimization of Average Consensus Containment with Multiple Stationary Leaders.](https://arxiv.org/abs/2203.16451)

## Abstract :

In this paper, we consider the problem of containment control of multi-agent systems with multiple stationary leaders, interacting over a directed network. While, containment control refers to just ensuring that the follower agents reach the convex hull of the leaders states, we focus on the problem where the followers achieve a consensus to the average values of the leaders states. We propose an algorithm that can be implemented in a distributed manner to achieve the above consensus among followers. Next we optimize the convergence rate of the followers to the average consensus by proper choice of weights for the interaction graph. This optimization is also performed in a distributed manner using Alternating Direction Method of Multipliers (ADMM). Finally, we complement our results by illustrating them with numerical examples.

## Usage : 

This implementation provides a setup as demonstrated in the paper having 24 agents (10 leaders and 14 followers), connected over a directed network topology. The topology has been chosen randomly such that it satisfies the assumptions listed in the paper.

## Requirements :

1. Code is good to go with MATLAB v2015b and above.
2. CVX package is required to be pre-installed with the MATLAB.

## Citation : 


```bibtex
@article{chatterjee2022distributed,
  title={Distributed Optimization of Average Consensus Containment with Multiple Stationary Leaders},
  author={Chatterjee, Sushobhan and Kalaimani, Rachel Kalpana},
  journal={arXiv preprint arXiv:2203.16451},
  year={2022}
}
```
