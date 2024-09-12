# Welcome!!

## Population Genetics

## Notes by Miles Garvin

## Week 3: Sep 10th / 12th

Where to find class data?:

```         
cd /gpfs1/cl/pbio3990
```

How to read zipped fq.gz files:

```         
zcat \<file.gz\> \| head
```

also:

```         
zcat \<file.gz\> \| head -n <number>
```

trying to get this to work... not the same in Rstudio and unix :(

```         
cat "/gpfs1/cl/pbio3990/PopulationGenomics/example_data/WV_9.R1.fq.gz"
```

library(tidyverse)

rnorm(n=1, mean = 10.2, sd = 3.8) rnorm(n=15, mean = 10.2, sd = 3.8)

how to load packages?????????:

```         
spack load <name>
```
