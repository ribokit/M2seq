#!/usr/bin/env python
###############################################################################
# Analysis of 2D signal in mutational profiling sequencing data (M2seq)
###############################################################################
#
# M2seq
#
# This script takes raw FASTQ files from a sequencing run and performs the following:
#
# 1. Demultiplexes the raw FASTQs using user-provided barcodes using Novobarcode
# 2. Uses the demultiplexed FASTQs and the WT sequence to generate 2D mutational profiling data
#       Use ShapeMapper for read alignment to reference sequence and generation of mutation strings
# 3. Calculates 2D datasets and outputs RDAT files using simple_to_rdat.py
#
# (C) Clarence Cheng, 2015-2016
# (C) Joseph Yesselman, Rhiju Das, 2017

import os, sys, time
import shutil
import argparse
import glob
import yaml
from string import Template

from utils import which, timeStamp, make_dir, get_sequence

M2SEQ_FILEPATH = os.path.realpath(__file__)
M2SEQ_FOLDER = os.path.dirname(M2SEQ_FILEPATH)
CFG_TEMPLATE_PATH = os.path.join(M2SEQ_FOLDER, 'ShapeMapper_config_template.txt')


def check_required_programs():
    print "checking required external programs"
    if which('novobarcode') is None:
        raise ValueError("novobarcode program is not installed but required!, "
                         "please install at: http://www.novocraft.com/support/download/")
    else:
        print "novobarcode is detected ..."


    if which("ShapeMapper.py") is None:
        raise ValueError("ShapeMapper 1.2 is required, software developed by the Weeks lab at UNC Chapel Hill for 1D "
                         "analysis ofmutational profiling data. Available at http://www.chem.unc.edu/rna/software.html "
                         "(Make sure you go into that directory and run make.)")
    else:
        print "ShapeMapper is detected ..."


    if which("bowtie2") is None:
        raise ValueError("BowTie2 is needed for ShapeMapper. Available here:"
                         " https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.9/. Version 2.2.9 works.")
    else:
        print "BowTie2 is detected ..."


def parse_commandline_args():
    parser = argparse.ArgumentParser()

    # parser.add_argument('sequencefile', type=argparse.FileType('r'))
    # parser.add_argument('barcodes', type=argparse.FileType('r'))
    # parser.add_argument('read1fastq', type=argparse.FileType('r'))
    # parser.add_argument('read2fastq', type=argparse.FileType('r'))
    # parser.add_argument('--config', type=argparse.FileType('r'))
    # parser.add_argument('--name', type=str, default='PLACEHOLDER')
    # parser.add_argument('--offset', type=int, default=0)
    # parser.add_argument('--outprefix', type=str, default='out')
    # parser.add_argument('--force_demultiplex', action='store_true')

    # parser.add_argument('--only_demultiplex', type=bool)

    parser.add_argument('--manifest')

    args = parser.parse_args()
    # if args.name == 'PLACEHOLDER':
    #     seq_file = args.sequencefile.name.split("/")[-1]
    #     name = seq_file.split(".")[0]
    #     args.name = name

    return args

def parse_manifest(manifest_path):
    """
    Parse the manifest YAML file at the specified path
    """
    with open(manifest_path, 'r') as f:
        manifest_data = yaml.load(f)
    
    print 'Loaded manifest file, it was read in as the following object data:'
    print manifest_data
    
    # Validate that the manifest is in the correct format
    validate_manifest(manifest_data)

    return manifest_data


    

def validate_manifest(manifest_data):
    pass


def get_barcode_sequences(args, f_log):
    print 'Reading barcodes from: ' + args.barcodes.name
    f_log.write('Primers:\n')
    lines = open(args.barcodes.name).readlines()
    primer_tags = []
    barcodes= {}
    print 'Current barcodes found:'
    for line in lines:
        if len(line) < 2: continue
        col = line.rstrip('\n')
        cols = col.split('\t')
        if len(cols) != 2: continue
        if len(cols[1]) < 2: continue
        primer_tags.append(cols[0])
        barcodes[cols[0]] = cols[1]
        f_log.write(line)
        print cols[0] + '\t' + cols[1]

    return barcodes


def demultiplex_fastq_files(read1_path, read2_path, RTB_file_path, output_folder):
    # f_log.write('Starting Novobarcode demultiplexing at: ' + timeStamp())
    print 'Starting Novobarcode demultiplexing'
    novobarcode_outpath = os.path.join(output_folder, '1_Demultiplex')
    novobarcode_log_path = os.path.join(novobarcode_outpath, 'novobarcode_log_Distance4.txt')
    novobarcode_command = 'novobarcode -b {} -f {} {} -d {} > {}'.format(RTB_file_path, \
    read1_path, read2_path, novobarcode_outpath, novobarcode_log_path)
    print novobarcode_command
    os.system(
        'novobarcode -b ' + RTB_file_path + ' -f ' + read1_path + ' ' + read1_path +  \
        ' -d 1_Demultiplex > 1_Demultiplex/novobarcode_log_Distance4.txt')
    # f_log.write('\nFinished demultiplexing at: ' + timeStamp() + '\n')


def valid_demultiplex_output(args, barcodes):
    read1fastq_name = args.read1fastq.name.split("/")[-1]
    read2fastq_name = args.read2fastq.name.split("/")[-1]
    for name, seq in barcodes.iteritems():
        fname_1 = '1_Demultiplex/' + seq + '/' + read1fastq_name
        fname_2 = '1_Demultiplex/' + seq + '/' + read2fastq_name
        if not os.path.isfile(fname_1) or not os.path.isfile(fname_2):
            return 0
    return 1


def create_sym_link_to_demultiplex_files(read1_path, read2_path, construct_barcodes, output_folder):
    read1fastq_name = read1_path.split("/")[-1]
    read2fastq_name = read2_path.split("/")[-1]
    print read1fastq_name
    print read1_path
    for tag, seq in construct_barcodes:
        fname_1 = '1_Demultiplex/' + seq + '/' + read1fastq_name
        fname_2 = '1_Demultiplex/' + seq + '/' + read2fastq_name
        new_fastq_names = [os.path.join(output_folder, tag + '_S1_L001_R1_001.fastq'),
                           os.path.join(output_folder, tag + '_S1_L001_R2_001.fastq')]
        os.system('ln -s ' + os.path.abspath(fname_1) + " " + os.path.abspath(new_fastq_names[0]))
        os.system('ln -s ' + os.path.abspath(fname_2) + " " + os.path.abspath(new_fastq_names[1]))


def run_shapemapper(shapemapper_folder, shapemapper_config_path):
    starting_folder = os.getcwd()
    os.chdir(shapemapper_folder)
    print 'Starting ShapeMapper analysis'
    os.system("ShapeMapper.py " +  shapemapper_config_path.split('/')[-1])

    # f_log.write('\nGenerating simple files at: ' + timeStamp())
    currdir = os.getcwd()
    outdir = os.path.join(currdir, 'output', 'mutation_strings_oldstyle')
    for file in os.listdir(outdir):
        if file.endswith('.txt'):
            muts_to_simple( os.path.join(outdir, file) )
    # f_log.write('\nFinished generating simple files at: ' + timeStamp())

    os.chdir(starting_folder)
    print os.getcwd()

def generate_sequence_fasta(name, sequence, save_folder):
    save_path = os.path.join(save_folder, name + '.fa')
    with open(save_path, 'w') as f:
        f.write('> {}\n'.format(name))
        f.write('{}'.format(sequence))
    return save_path

def generate_shapemapper_config(construct_name, construct_barcodes, config_template):
    pass

def valid_shapemapper_output(args):
    pass

def muts_to_simple(mutsfile):
    simplename = mutsfile + '.simple'
    f_simple = open(simplename, 'w')

    start_pos = []
    count = 0
    f = open(mutsfile)
    lines = f.readlines()
    f.close()

    for line in lines:
        fields = line.strip().split('\t')
        if len(fields) < 3: continue

        start_pos = int(fields[0])
        end_pos = int(fields[1])
        mut_string = fields[2]
        assert (len(mut_string) == (int(fields[1]) - int(fields[0]) + 1))

        count += 1  # record total sequences
        if count % 50000 == 0: print 'Reading line number ', count

        # ignore non-overlapping reads
        temp = mut_string.strip('s~')
        non_overlap = 0
        for i in xrange(len(temp)):
            if temp[i] == 's' or temp[i] == '~':
                if temp[i + 1] == '|':
                    non_overlap = 1
                    break

        # adjust start and end positions
        if non_overlap == 0:
            for i in xrange(len(mut_string)):
                if mut_string[i] == 's' or mut_string[i] == '~':
                    start_pos += 1
                else:
                    break

            for i in reversed(xrange(len(mut_string))):
                if mut_string[i] == 's' or mut_string[i] == '~':
                    end_pos -= 1
                else:
                    break

            simple_line = mut_string.strip('s~').replace('|', '0').replace('~', '0').replace('A', '1').replace('T','1').replace('G', '1').replace('C', '1').replace('-', '1')

            f_simple.write(str(start_pos) + '\t' + str(end_pos) + '\t' + simple_line + '\n')

    print '\nSimple file created: ' + f_simple.name
    print '\nTotal number of sequences: ' + str(count)

    f_simple.flush()
    f_simple.close()


def m2_seq_final_analysis(construct_name, offset=0):
    print 'Starting M2seq analysis'
    # f_log.write('\nStarting M2seq analysis at: ' + timeStamp() + '\n')
    currdir = os.getcwd()
    make_dir(os.path.join('3_M2seq', construct_name, 'simple_files'))

    print os.getcwd()
    simple_files = glob.glob(os.path.join('2_ShapeMapper', construct_name, 'output/mutation_strings_oldstyle/*.simple'))
    for sf in simple_files:
        print sf
        f_name = sf.split("/")[-1]
        os.system('ln -s ' + os.path.abspath(sf) + " " + '3_M2seq/' + construct_name + '/simple_files/' + f_name)

    
    # Copy over our fasta file from the shapemapper folder to the m2seq folder
    shutil.copy(os.path.join('2_ShapeMapper', construct_name, construct_name + '.fa'), os.path.join('3_M2seq', construct_name))
    fasta_path = os.path.join('3_M2seq', construct_name, construct_name + '.fa')

    os.chdir(currdir + '/3_M2seq/' + construct_name + '/simple_files')
    simple_files = glob.glob('*.simple')
    for sf in simple_files:
        f_name = sf.split('.')[0]
        print sys.executable
        os.system(sys.executable + ' ' + M2SEQ_FOLDER + "/simple_to_rdat.py " + fasta_path + " --simplefile " + sf + " --name " +
                  construct_name + " --offset " + str(offset) + " --outprefix " + f_name)
    os.chdir(currdir)

def single_sample_analysis(args):
    currdir = os.getcwd()
    with open(os.path.join(currdir, 'AnalysisLog.txt'), 'w') as f_log:
        # make directories now
        make_dir(currdir + '/1_Demultiplex')
        make_dir(currdir + '/2_ShapeMapper')
        make_dir(currdir + '/3_M2seq')
        make_dir(currdir + '/3_M2seq/simple_files')

        # TODO move this into settings.py
        file_path = os.path.realpath(__file__)
        spl = file_path.split("/")
        base_dir = "/".join(spl[:-1])

        sequence = get_sequence(args.sequencefile)
        barcodes = get_barcode_sequences(args, f_log)

        # should we demultiplex?
        if not valid_demultiplex_output(args, barcodes) or args.force_demultiplex:
            demultiplex_fastq_files(args, f_log)
            # create links instead of actually moving files, much cleaner
            create_sym_link_to_demultiplex_files(args, barcodes)
        else:
            print "not demultiplexing it's done already, use --force_demultiplex to redo it"

        # run shaper mapper
        run_shapemapper(args)

        # run m2seq analysis
        m2_seq_final_analysis(args, f_log)


def generate_novobarcode_file(manifest_data, outpath):
    distance = 4
    format = 5
    barcodes = []
    for construct in manifest_data['constructs']:
        for sample in construct['samples']:
            barcodes.append([sample['barcode_tag'], sample['barcode']])

    with open(outpath, 'w') as outfile:
        outfile.write('Distance:\t{}\n'.format(str(distance)))
        outfile.write('Format:\t{}\n'.format(str(format)))  
        for barcode in barcodes:
            outfile.write('{}\n'.format('\t'.join(barcode)))

    return barcodes


def manifest_analysis(manifest_path):
    currdir = os.getcwd()
    make_dir(currdir + '/1_Demultiplex')
    make_dir(currdir + '/2_ShapeMapper')
    make_dir(currdir + '/3_M2seq')
    make_dir(currdir + '/3_M2seq/simple_files')

    # Load in the manifest data
    manifest_data = parse_manifest(manifest_path)
    read1_path = manifest_data['undemultiplexed']['read1']
    read2_path = manifest_data['undemultiplexed']['read2']
    output_folder = manifest_data['output_folder']

    # Generate barcode file
    novobarcode_file_path = os.path.join(manifest_data['output_folder'], 'RTBbarcodes.fa')
    barcodes = generate_novobarcode_file(manifest_data, novobarcode_file_path)
    print 'Generated novobarcode RTB file and stored it at ' + novobarcode_file_path

    # Do the demultiplexing
    # demultiplex_path = os.path.join(manifest_data['output_folder'], '1_Demultiplex')
    # demultiplex_fastq_files(read1_path, \
    #                         read2_path, \
    #                         novobarcode_file_path, \
    #                         manifest_data['output_folder'])

    # Run SHAPEmapper for each construct
    for construct in manifest_data['constructs']:
        # Start by creating the necessary folder and symlinks into the 2_ShapeMapper folder
        construct_sequence = construct['sequence']
        construct_name = construct['name']
        construct_barcodes = []
        for sample in construct['samples']:
            construct_barcodes.append([sample['barcode_tag'], sample['barcode']])
        
        shapemapper_output_folder = os.path.join('2_ShapeMapper', construct_name)
        make_dir(shapemapper_output_folder)

        create_sym_link_to_demultiplex_files(read1_path, read2_path, construct_barcodes, shapemapper_output_folder)

        # Create ShapeMapper .cfg file for the construct
        # We do this by importing a string template from an external file, then replacing the
        # specified config lines with our own, as generated for this given construct
        with open(CFG_TEMPLATE_PATH, 'r') as f:
            cfg_template = Template(f.read())

        construct_cfg_lines_template = '{0}: {0}_S1_L001_R1_001.fastq, {0}_S1_L001_R2_001.fastq = {1}\n'
        construct_cfg_lines = ''
        for sample in construct['samples']:
            construct_cfg_lines += construct_cfg_lines_template.format(sample['barcode_tag'], construct['name'])

        cfg_contents = cfg_template.substitute({'cfg_lines':construct_cfg_lines})

        cfg_path = os.path.join(shapemapper_output_folder, 'ShapeMapper_config.cfg')
        with open(cfg_path, 'w') as f:
            f.write(cfg_contents)

        # Generate sequence fasta in the construct folder
        generate_sequence_fasta(construct_name, construct_sequence, shapemapper_output_folder)

        # Run shapemapper
        run_shapemapper(shapemapper_output_folder, cfg_path)

        # Run the final M2-seq analysis
        # m2_seq_final_analysis(construct_name)



if __name__ == "__main__":
    check_required_programs()

    args = parse_commandline_args()

    if args.manifest is not None:
        manifest_analysis(args.manifest)
    else:
        single_sample_analysis(args)


