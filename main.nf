// Name
params.name = "mBed_Distance"

// input sequences to align [FASTA]
params.seqs = "$baseDir/tutorial/PF00004.fa"
//params.seqs = "$baseDir/tutorial/seatoxin.fa"

// input reference sequences aligned [Aligned FASTA]
params.refs = "$baseDir/tutorial/PF00004.seed"
//params.refs = "$baseDir/tutorial/seatoxin.ref"

// output directory [DIRECTORY]
params.output = "$baseDir/results"


log.info """\
         D P A   A n a l y s i s  ~  version 0.1"
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
  set val(id),file(sequences),file(references) from seqsAndRefs

  output:
  set val(id), file(sequences), \
      file("${id}_populated.ref") \
      into filledRefs 
  script:
    """
    $baseDir/bin/populate_fasta.py ${sequences} ${references} ${id}_populated.ref
    """
}

process removeSeeds{
  tag "${id}"
  
  input:
  set val(id), file(sequences), \
      file(references) \
      from filledRefs
  
  output:
  set val(id), file("${id}_clean.fa"), \
      file(references) \
      into filledAndRemoved
  script:
    """
    $baseDir/bin/SeqFilter/bin/SeqFilter ${sequences} --ids $params.refs --ids-exclude --out ${id}_clean.fa -q >>/dev/null

    """
}

process mergeFiles{
  tag "${id}"

  input:
  set val(id), file(sequences), \
      file(references) \
      from filledAndRemoved

  output:
  set val(id), file("${id}_complete"), \
      file(references) \
      into mergedFiles

  script:
    """
    cat ${references} ${sequences} > ${id}_complete
    """
}

process mBed {
  tag "${id}"
  publishDir "${params.output}", mode: 'copy', overwrite: true

  input:
  set val(id), file(sequences), \
      file(references) \
      from mergedFiles

  output:
  file("${id}_coordinates.out") \
      into finalResult

  script:
    """
    set +e
    numRef=\$(cat $params.refs | wc -l)
    $baseDir/bin/mBed/mBed -infile ${sequences} -seedfile ${references} -method SeedMap -numInputSeeds \$numRef >>/dev/null

    foo=\$?
    mv coordinates.out ${id}_coordinates.out
    """
}









