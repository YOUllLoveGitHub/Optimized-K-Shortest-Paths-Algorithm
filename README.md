# Optimized-K-Shortest-Paths

The OKSP algorithm efficiently solves the three-tier caching problem, involving users, nearby server cluster, and cloud server. The time complexity of the algorithm is $O((K\cdot M + N) \cdot M \cdot \log N))$, where $K$ is the total number of tiles that all $M$ MEC servers can cache in collaboration, and $N$ is the number of different tile files. 

And we added the **output_lp_file** function to the code. This function outputs the linear programming input file of the original problem, and the lp solver can be called to verify the correctness of the solution given by the OKSP algorithm. 

To use the lp solver, we downloaded the [winglpk-4.65](https://jaist.dl.sourceforge.net/project/winglpk/winglpk/GLPK-4.65/winglpk-4.65.zip) package and followed these steps:

1. Download glpk.
2. Open a command prompt and navigate to the `winglpk-4.65\glpk-4.65\w32` directory.
3. Run the command `glpsol --lp 'lp.input' -o 'test.txt'`.

In this command, "lp.input" is the output file of the **output_lp_file** function. The solver will give the exact solution using the lp algorithm.

We are currently in the process of submitting a paper that provides a detailed proof of the algorithm's correctness. We plan to add more information about the proof in the future.
