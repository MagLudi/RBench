RBench
======

Benchmark hardware using R.

The benchmark() function generates a random dataset, then performs multiple density-based clusterings
using the dbscan implementation from the fpc package. By default, the number of datapoints being
analyzed is 1800. Use the 'size' parameter to scale this up or down (e.g. benchmark(size=2) would
use a dataset with 3600 points). 

Here is the output from a very old PC with a 3.2 GHz Pentium and 3 GB ram running XP: 

> benchmark(size=1.5)

> user  system elapsed

> 121.65    0.42  114.24 
 
> benchmark(size=0.5)

>  user  system elapsed
   
> 13.69    0.00   12.91 
  
> benchmark(size=0.5)

> user  system elapsed 
   
> 13.71    0.00   12.95  
