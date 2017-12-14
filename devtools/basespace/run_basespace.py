#!/usr/bin/env python
# encoding: utf-8

from __future__ import (absolute_import, division, print_function)
from collections import defaultdict
import subprocess as sp
import os
import sys
import errno
import argparse
import time
import json


def parse_args():
    class FullPaths(argparse.Action):
        """Expand user- and relative-paths"""
        def __call__(self, parser, namespace, values, option_string=None):
            setattr(namespace, self.dest, os.path.abspath(os.path.expanduser(values)))

    usage = " "
    parser = argparse.ArgumentParser(description="STAR-SEQR Parameters:", epilog=usage)
    parser.add_argument('-j', '--json', type=str, required=False, action=FullPaths,
                        default='/data/input/AppSession.json',
                        help='json input file from basespace app')
    parser.add_argument('-s', '--sample_path', type=str, required=False, action=FullPaths,
                        default='/data/input/samples/',
                        help='default path for samples')
    parser.add_argument('-r', '--results_path', type=str, required=False, action=FullPaths,
                        default='/data/output/appresults/',
                        help='default path for results')
    parser.add_argument('-t', '--scratch_path', type=str, required=False, action=FullPaths,
                        default='/data/scratch/',
                        help='default path for scratch')
    args = parser.parse_args()
    return args


def make_new_dir(newdir):
    try:
        os.mkdir(newdir)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise e
        pass


def check_file_exists(path):
    if (os.stat(os.path.realpath(path)).st_size == 0):
        print("Exiting. Cannot find file: " + os.path.abspath(path))
        sys.exit(1)


def read_json(in_file):
    with open(in_file) as data_file:
        data = json.load(data_file)
    return data


def extract_projectid(raw_cfg):
    for item in raw_cfg['Properties']['Items']:
            if item['Name'] == 'Input.Projects':
                project_id = item['Items'][0]['Id']
    return project_id


def extract_params(raw_cfg):
    '''This method is used to extract parameters specific to the method.'''
    avoid_items = ['app-session-name', 'basespace-labs-disclaimer',
                   'project-id.attributes', 'project-id',
                   'sample-id.attributes', 'sample-id',
                   'Projects', 'Samples']
    parameters = defaultdict(dict)
    for item in raw_cfg['Properties']['Items']:
        # add parameters to parameters list
        item_name = item['Name']
        if item_name.startswith('Input.'):
            item_key = str(item_name.replace('Input.', ''))
            if item_key not in avoid_items:
                if 'Content' in item.keys():
                    item_val = item['Content']
                    parameters[item_key] = str(item_val)
                else:  # used for flags/radio buttons
                    parameters[item_key] = item_key
    return parameters


def extract_samples(raw_cfg, args):
    '''This funtion extracts the samples to be run from the json file'''
    sinfo = defaultdict(dict)
    for item in raw_cfg['Properties']['Items']:
        if item['Name'] == 'Input.Samples':
            for sample in item['Items']:
                sname = str(sample['Id'])
                sinfo[sname]['Id'] = str(sample['Id'])
                sinfo[sname]['folder'] = os.path.join(args.sample_path, str(sample['Id']))
    return sinfo


def create_results_folders(args, res):
    '''
    This function create the sample-level folders for writing results to later.results_path.
    The folder structure is /data/output/appresults/<project-id>/<sample-id>/
    '''
    base = os.path.join(args.results_path, res['projid'])
    make_new_dir(base)
    for sample in res['samples']:
        output_dir = os.path.join(base, sample)
        make_new_dir(output_dir)


def get_pairfq_paths(res, extension=[".fastq", ".fastq.gz", ".fq", ".fq.gz"]):
    '''
    This function finds the fastq files in the sample-level folder.
    Currently only works for one set of paired fastqs.
    '''
    for sample in res['samples']:
        fq1 = ''
        fq2 = ''
        for root, dirs, files in os.walk(res['samples'][sample]['folder']):
            for file in files:
                if file.endswith(tuple(extension)):
                    if "_R1_" in file:
                        fq1 = os.path.join(root, file)
                    elif "_R2_" in file:
                        fq2 = os.path.join(root, file)
        check_file_exists(fq1)
        check_file_exists(fq2)
        res['samples'][sample]['fq1'] = fq1
        res['samples'][sample]['fq2'] = fq2
    return res


def get_indexes(res, args):
    '''
    Downloads the index files from AWS s3 bucket and places in scratch folder
    The scratch folder is recommended for analysis as it has no size restrictions
    '''
    if res['params']['species'] == "human":
        cmd = ['aws', 's3', 'sync', 's3://starseqr-inputs/human/gencodev25lift37', args.scratch_path, '--no-sign-request']
        run_cmd = list(map(str, cmd))
        print("Command: " + " ".join(run_cmd))
        try:
            p = sp.Popen(run_cmd, stdout=sp.PIPE, stderr=sp.PIPE)
            stdout, stderr = p.communicate()
            if stdout:
                print(stdout)
            if stderr:
                print(stderr)
        except OSError as e:
            print('downloading indexes failed')
            print('Exception: ' + str(e))
            traceback = sys.exc_info()[2]
            raise(ValueError, e, traceback)
        # set paths
        res['params']['gtf'] = os.path.join(args.scratch_path, 'gencode.v25lift37.annotation.gtf.gz')
        res['params']['ref'] = os.path.join(args.scratch_path, 'GRCh37.primary_assembly.genome.fa.gz')
        res['params']['idx'] = os.path.join(args.scratch_path, 'STAR_INDEX')
        check_file_exists(res['params']['gtf'])
        check_file_exists(res['params']['ref'])
        check_file_exists(res['params']['idx'])
    else:
        print("No species found!")
        sys.exit(1)
    return res


def get_file_paths(res, fkey, extension):
    '''General function that will find a file type based on extensions'''
    for sample in res['samples']:
        query_file = ''
        for root, dirs, files in os.walk(res['samples'][sample]['folder']):
            for file in files:
                if file.endswith(tuple(extension)):
                    query_file = str(os.path.join(root, file))
        check_file_exists(query_file)
        res['samples'][sample][fkey] = query_file
    return res


def run_starseqr(res):
    '''This runs starseqr'''
    for sample in res['samples']:
        fq1 = res['samples'][sample]['fq1']
        fq2 = res['samples'][sample]['fq2']
        gtf = res['params']['gtf']
        ref = res['params']['ref']
        idx = res['params']['idx']
        msens = res['params']['sensitivity']
        cmd = ['starseqr.py', '-1', fq1, '-2', fq2,
               '-g', gtf, '-r', ref, '-i', idx,
               '-p', sample, '-m', msens]
        run_cmd = list(map(str, cmd))
        print("Command: " + " ".join(run_cmd))
        try:
            p = sp.Popen(run_cmd, stdout=sp.PIPE, stderr=sp.PIPE)
            stdout, stderr = p.communicate()
            if stdout:
                print(stdout)
            if stderr:
                print(stderr)
        except OSError as e:
            print('starseqr Failed')
            print('Exception: ' + str(e))
            traceback = sys.exc_info()[2]
            raise(ValueError, e, traceback)
    return


def copy_results(res):
    '''
    Used to copy the final results to output directory for upload
    '''
    pass


def main():
    start = time.time()
    args = parse_args()
    os.chdir(args.scratch_path)
    with open('changed_dir.txt', 'w') as out1:
        out1.write("Done")
    raw_cfg = read_json(args.json)
    with open('read_json.txt', 'w') as out1:
        out1.write("Done")
    # store all necessary items in the res dict. Update along the way
    res = {}
    res['projid'] = extract_projectid(raw_cfg)
    res['params'] = extract_params(raw_cfg)
    res['samples'] = extract_samples(raw_cfg, args)
    with open('ex_samples.txt', 'w') as out1:
        out1.write("Done")
    res = get_pairfq_paths(res)
    with open('get_fq.txt', 'w') as out1:
        out1.write("Done")
    # if 'bed' in res['params'].keys():
    #     res = get_file_paths(res, 'bed', ['.bed'])
    res = get_indexes(res, args)
    with open('parsed_params.txt', 'w') as outfile2:
        json.dump(res, outfile2)
    run_starseqr(res)
    create_results_folders(args, res)
    copy_results(res)
    print("Program took  %g seconds" % (time.time() - start))


if __name__ == "__main__":
    sys.exit(main())
