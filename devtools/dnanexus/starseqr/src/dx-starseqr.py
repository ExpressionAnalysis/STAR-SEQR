#!/usr/bin/env python
# STAR-SEQR 0.0.1

import os
import shlex
import subprocess
import tempfile
from glob import glob
import magic
import traceback
import time
import dxpy


def download_link(dxlink):
	"""Download a DNAnexus object reference to the original filename.
	If the object is actually None, return None.
	"""
	if dxlink is not None:
		start = time.time()
		dxf = dxpy.DXFile(dxlink)
		fname = dxf.describe()['name']
		print("Downloading " + str(fname))
		if ' ' in fname:
			fname_nospace = fname.replace(' ', '_')
			print("Filename", repr(fname), "contains spaces; renaming to",
				  repr(fname_nospace))
			fname = fname_nospace
		dxpy.download_dxfile(dxf.get_id(), fname)
		print("Download took %g seconds" % (time.time() - start))
		return fname


def maybe_gunzip(fname, base, ext):
	if fname and 'gzip' in magic.from_file(fname):
		start = time.time()
		print("Gunzip file " + str(fname))
		newf = safe_fname(base, ext)
		sh("gunzip", fname, "-c >", newf)
		fname = newf
		print("Gunzip took %g seconds" % (time.time() - start))
	return fname


def safe_fname(base, ext):
	"""Ensure the chosen filename will not overwrite an existing file.
	If it would, insert some garbage in the middle of the chosen filename to
	avoid the naming conflict.
	"""
	fname = base + "." + ext
	if os.path.exists(fname):
		fname = tempfile.mktemp(prefix=base + "-", suffix="." + ext, dir=".")
	return fname


def sh(*command):
	"""Run a shell command."""
	cmd = " ".join(map(str, command))
	print("$", cmd)
	try:
		subprocess.check_call(cmd, shell=True)
	except subprocess.CalledProcessError as exc:
		check_files(command[1:])
		raise exc
	print()


def check_files(maybe_filenames):
	fnames = []
	for fname in maybe_filenames:
		if isinstance(fname, basestring):
			fnames.extend(shlex.split(fname))
	for fname in fnames:
		if '.' in fname:
			# It might be a filename
			if os.path.exists(fname):
				print("File:", os.path.abspath(fname))
			else:
				print("Not file:", os.path.abspath(fname))


def starseqr_docker(*args):
	"""Run docker"""
	docker_prefix = ["dx-docker", "run",
					 "-v", "/home/dnanexus:/data", "-w", "/data",
					 "eagenomics/starseqr:0.6.5", "starseqr.py"]
	docker_cmd = docker_prefix + list(args)
	print(" ".join(map(str, docker_cmd)))
	sh(*(docker_cmd))


def starseqr_docker_index(*args):
	"""Run docker"""
	docker_prefix = ["dx-docker", "run",
					 "-v", "/home/dnanexus:/data", "-w", "/data",
					 "eagenomics/starseqr:0.6.6", "samtools", "faidx"]
	docker_cmd = docker_prefix + list(args)
	print(" ".join(map(str, docker_cmd)))
	sh(*(docker_cmd))


@dxpy.entry_point('main')
def main(mode, prefix, genome_fasta, transcript_gtf, threads, star_index=None, fastq1=None, fastq2=None, chimeric_juncs=None, chimeric_sam=None , chimeric_bam=None):
	print("Starting python script")

	# Check necessary params and combinations exist
	# align = [star_index]
	# call = [chimeric_juncs, chimeric_sam, chimeric_bam]
	# fqs = [fastq1, fastq2]
	# star_reads = [chimeric_sam, chimeric_bam]
	# if any(align) and any(call):
	#     print("Error: Please specify either a STAR Index or STAR existing files as input!")
	#     sys.exit(1)
	# if any(align) and None in fqs:
	#     print("Error: Fastq1, Fastq2, and the STAR index must be specified if doing alignment")
	#     sys.exit(1)
	# if any(call) and not any(star_reads):
	#     print("Error: The STAR .junction and .sam|.bam file must be specified if using existing alignment")
	#     sys.exit(1)
	# if all(star_reads):
	#     print("Error: Please specify only one of .sam|.bam .")
	#     sys.exit(1)

	# Download data and convert variables to filenames
	## Required data
	fastq1 = download_link(fastq1)
	fastq2 = download_link(fastq2)
	genome_fasta = download_link(genome_fasta)
	genome_fasta = maybe_gunzip(genome_fasta, "genome", "fa")
	transcript_gtf = download_link(transcript_gtf)
	transcript_gtf = maybe_gunzip(transcript_gtf, "transcripts", "gtf")

	## Optional data
	chimeric_juncs = download_link(chimeric_juncs)
	chimeric_juncs = maybe_gunzip(chimeric_juncs, "this.Chimeric.out", "junction")
	chimeric_sam = download_link(chimeric_sam)
	chimeric_sam = maybe_gunzip(chimeric_sam, "this.Chimeric.out", "sam")
	chimeric_bam = download_link(chimeric_bam)
	star_index = download_link(star_index)
	if star_index:
		star_index_dir = star_index[:-len(".tar.gz")]
		try:
			print("Preparing STAR Index directory")
			sh("mkdir -p " + str(star_index_dir))
			subprocess.check_call(["tar", "-zxvf", str(star_index), "-C", str(star_index_dir)])
		except Exception:
			traceback.print_exc()
			raise dxpy.AppError("Error extracting reference tarball")

	# Build cmd
	ss_cmd = ["-p", prefix, "-t", threads, "-m", mode, "-vvv"]
	ss_cmd.extend(["-1", fastq1])
	ss_cmd.extend(["-2", fastq2])
	ss_cmd.extend(["-r", genome_fasta])
	ss_cmd.extend(["-g", transcript_gtf])

	if chimeric_juncs:
		ss_cmd.extend(["-sj", chimeric_juncs])
	if chimeric_sam:
		ss_cmd.extend(["-ss", chimeric_sam])
	if chimeric_bam:
		ss_cmd.extend(["-sb", chimeric_bam])
	if star_index:
		ss_cmd.extend(["-i", star_index_dir])

	## Run STAR-SEQR with Docker
	sh("ls -Altr")  # Show the generated files in the DNAnexus log
	print("Running STARSEQR")
	starseqr_docker_index(*[genome_fasta])
	starseqr_docker(*ss_cmd)
	sh("ls -Altr")  # Show the generated files in the DNAnexus log

	# Upload the output
	print("Uploading output")
	ss_breakpoint_bedpe = os.path.join(prefix + "_STAR-SEQR", prefix + "_STAR-SEQR_breakpoints.bedpe")
	ss_breakpoint_text = os.path.join(prefix + "_STAR-SEQR", prefix + "_STAR-SEQR_breakpoints.txt")
	breakpoint_bedpe = dxpy.upload_local_file(ss_breakpoint_bedpe)
	breakpoint_text = dxpy.upload_local_file(ss_breakpoint_text)

	# Link the output to app
	output = {}
	output["breakpoint_bedpe"] = dxpy.dxlink(breakpoint_bedpe)
	output["breakpoint_text"] = dxpy.dxlink(breakpoint_text)

	return output

start = time.time()
dxpy.run()
print("Program took  %g seconds" % (time.time() - start))
