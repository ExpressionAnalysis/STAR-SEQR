#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

requirements:
  - class: InlineJavascriptRequirement
  - class: InitialWorkDirRequirement
    listing:
      - entry: $(inputs.transcript_gtf)
        entryname: $(inputs.transcript_gtf.basename)
      - entry: $(inputs.genome_fasta)
        entryname: $(inputs.genome_fasta.basename)

hints:
  - class: DockerRequirement
    dockerPull: eagenomics/starseqr:0.6.1
  - class: ResourceRequirement
    coresMin: 8
    ramMin: 50000

baseCommand: ['starseqr.py']

inputs:
  star_index_dir:
    type: Directory?

  fq1:
    type: File?
    inputBinding:
      prefix: '-1'
      position: 1

  fq2:
    type: File?
    inputBinding:
      prefix: '-2'
      position: 2

  star_sam:
    type: File?
    inputBinding:
      prefix: '-ss'
      position: 2

  star_bam:
    type: File?
    inputBinding:
      prefix: '-sb'
      position: 2

  star_juncs:
    type: File?
    inputBinding:
      prefix: '-sj'
      position: 2

  genome_fasta:
    type: File?

  transcript_gtf:
    type: File?

  name_prefix:
    type: string
    inputBinding:
      prefix: -p
      position: 3

  worker_threads:
    type: int?
    default: 8
    inputBinding:
      prefix: -t
      position: 3

  mode:
    type: int?
    default: 1
    inputBinding:
      prefix: -m
      position: 3

outputs:
  starseqr_bedpe:
    type: File
    outputBinding:
      glob: $(inputs.name_prefix + '_STAR-SEQR' + '/' + inputs.name_prefix + '_STAR-SEQR_breakpoints.bedpe')
  starseqr_text:
    type: File
    outputBinding:
      glob: $(inputs.name_prefix + '_STAR-SEQR' + '/' + inputs.name_prefix + '_STAR-SEQR_breakpoints.txt')

arguments:
  - valueFrom: $(inputs.transcript_gtf.basename)
    position: 3
    prefix: '-g'

  - valueFrom: $(inputs.genome_fasta.basename)
    position: 3
    prefix: '-r'