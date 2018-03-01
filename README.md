# mBed Pipeline for single distance

[![https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg](https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg)](https://singularity-hub.org/collections/624)

To run this mBed distance pipeline, you need to have the input data in the correct format.

You will need 2 files. foo.fa and foo.seed ,the extensions could be whatever, but the file name has to be the same for both of it (same family)

In the fasta file, you will have all the sequences (seed+non seed) and in the foo.seed file you will need to have the name of the seeds, one name per line.

Then, you just need to run nextflow with 
```
nextflow run main.nf --seqs=<fasta file> --refs=<seed file> -with-singularity
```

By default it is set to the location ```./tutorial/*.fa``` and ```./tutorial/*.seed```

The output will be 3 files:

**./results/\<seqID\>_seed.out** where you can find the distance between the seeds against himself

**./results/\<seqID\>_clean.out** where you can find the distance of your family against the seeds (seeds are not printed here)

**./results/\<seqID\>_coordinates.out** where you can find all the distances

