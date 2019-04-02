# Cannon-s-algorithm-project
Project of parallel computation of Matrix multiplication. This was part of the course Parallel Computations for Large-Scale Problems at KTH.

## Planning

27/3: Discussed practicalities + DLs + strategies to deal with non square matrices and non square processor grids.

3/4: Further and final discussions on non square cases.

12/4: Itermediate meeting. The implementation of the parallel algorithm should be done.

25/4: The Timings should be done and we should have most of the plots ready

9/5: DL report

16-17/5: Oral presentation

## Main tasks

* Write a serial algorithm (Usual cache oblivious implementation and Strassen algorithm)
* Write the parallel algorithm with MPI
* Compute the theoretical complexity of our algorithm
* Time the computation by taking into account the broadcast/gathering parts and without.
* Analyse speedup and efficiency for different sizes of matrices and processor grid (interesting case when data transfer or computation in each processor becomes too important: bottle neck problem)
* Write the report ...

## Discussions on non square examples and matrices that can't be equally divided by the number of processors.

We thought about two way to deal with this: complete with zeros or truncate the matrices and do the remaining calculus another way.

The last approach we thought of was:
* Decide to split the number R of processors in R=PxQ where P and Q satisfy the following requirements.
* If we multiply A and B we can have different heights for the partitions of A and different wiwth for the patrtition of B. However P should divide the #column of A and Q should divide #row of B.

NB: We should check to what extend this solution works. Can we always find a P and Q that works? Will the final partition be clever?
