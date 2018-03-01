// Name
params.name = "mBed_Distance"

// input sequences to align [FASTA]
params.seqs = "$baseDir/tutorial/*.fa"
//params.seqs = "$baseDir/tutorial/seatoxin.fa"

// input reference sequences aligned [Aligned FASTA]
params.refs = "$baseDir/tutorial/*.seed"
//params.refs = "$baseDir/tutorial/seatoxin.seed"


// output directory [DIRECTORY]
params.output = "$baseDir/results"


log.info """\
         m B e d   A n a l y s i s  ~  version 0.1"
         ======================================="
         Name                                                  : ${params.name}
         Input sequences (FASTA)                               : ${params.seqs}
         Input references (Aligned FASTA)                      : ${params.refs}
         Output directory (DIRECTORY)                          : ${params.output}
         """
         .stripIndent()

// Channels for sequences [REQUIRED]
Channel
  .fromPath(params.seqs)
  .ifEmpty{ error "No files found in ${params.seqs}"}
  .map { item -> [ item.baseName, item] }
  .into { seqs; seqs2}

// Channels for sequences [REQUIRED]
Channel
  .fromPath(params.refs)
  .ifEmpty{ error "File NOT found in ${params.seqs}"}
  .map { item -> [ item.baseName, item] }
  .into { seeds; seeds2  }

seqs
    .combine( seeds, by: 0 )
    .set { seqsAndRefs }

process fillRefs {
  tag "${id}"
  
  input:
  set val(id), \
      file(sequences), \
      file(references) from seqsAndRefs

  output:
  set val(id),\
      file(sequences), \
      file (references), \
      file("${id}_populated.ref") \
    into filledRefs
 
  script:
    """
    ${baseDir}/bin/populate_fasta.py ${sequences} ${references} ${id}_populated.ref
    """
}

process removeSeeds{
  tag "${id}"
  
  input:
  set val(id),\
      file(sequences), \
      file(references), \
      file (ref_populated) \
    from filledRefs
  
  output:
  set val(id), file("${id}_clean.fa"), \
      file(references), \
      file(ref_populated) into filledAndRemoved
  
  script:
    """
    /mBed-nf/bin/SeqFilter/bin/SeqFilter ${sequences} --ids ${ref_populated} --ids-exclude --out ${id}_clean.fa -q >>/dev/null
    """
}

process mergeFiles{
  tag "${id}"

  input:
  set val(id),\
      file(sequences), \
      file(references), \
      file(ref_populated)\
    from filledAndRemoved
  
  output:
  set val(id), file("${id}_complete"), \
      file(references), \
      file(ref_populated)\
    into mergedFiles

  script:
    """
    cat ${ref_populated} ${sequences} > ${id}_complete
    """
}

process mBed {
  tag "${id}"

  publishDir "${params.output}", mode: 'copy', overwrite: true

  input:
  set val(id), file(sequences), \
      file(references), \
      file(ref_populated)\
    from mergedFiles
  
  output:
  set file("${id}_seed.out"), \
      file("${id}_coordinates.out"), \
      file("${id}_clean.out") \
    into finalResult
  
  script:
    """
    set +e
    numRef=\$(cat ${references} | wc -l)
    /mBed-nf/bin/mBed/mBed -infile ${sequences} -seedfile ${ref_populated} -method SeedMap -numInputSeeds \$numRef >>/dev/null
    foo=\$?
    
    mv coordinates.out ${id}_coordinates.out
    awk 'FNR==NR {a[\$1]=\$0; next}; \$1 in a {print a[\$1]}' ${id}_coordinates.out ${references} > ${id}_seed.out
    awk 'NR==FNR{a[\$0]=1;next}!a[\$0]' ${id}_seed.out ${id}_coordinates.out > ${id}_clean.out
    """
}

