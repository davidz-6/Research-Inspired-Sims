
# MATLAB Simulations

A collection of MATLAB simulations I use for research, primarily exploring networked dynamical systems, synchronization, and distributed control.

Most simulations consider networks of coupled agents of the form

$$
\dot{x}_i = f(x_i) + c\sum_{j=1}^{N} a_{ij}\Gamma(x_j-x_i) + u_i
$$

where $a_{ij}$ describes the network topology, $c$ is the coupling strength, and $u_i$ represents the control input.

Some simulations explore event-triggered control, where updates occur only when a triggering condition is satisfied, rather than at fixed time intervals:

$$
\|e_i(t)\|^2 \geq \sigma_i \|x_i(t)\|^2
$$

The code is mainly intended for numerical experimentation and visualization of theoretical results.

## References

- References cited within the individual simulations.
- P. Tabuada, *Event-Triggered Real-Time Scheduling of Stabilizing Control Tasks*, IEEE Transactions on Automatic Control, 2007.
- W. Yu et al., research on pinning control and synchronization of complex dynamical networks.
