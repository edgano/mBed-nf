# mBed Pipeline for single distance

To run this mBed distance pipeline, you need to have the input data in the correct format.

You will need 2 files. foo.fa and foo.seed ,the extensions could be whatever, but the file name has to be the same for both of it (same family)

In the fasta file, you will have all the sequences (seed+non seed) and in the foo.seed file you will need to have the name of the seeds, one name per line.

Then, you just need to run nextflow with 
```
nextflow run main.nf --seqs=<fasta file> --refs=<seed file> -with-singularity
```

By default it is set to the location ```./tutorial/*.fa``` and ```./tutorial/*.seed```

The output file will be in ```./results/<seqID>_coordinates.out```

