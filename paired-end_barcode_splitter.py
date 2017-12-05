#!/usr/bin/env python

'''
USAGE
    (file pair)
    paired-end_barcode_splitter.py --bcfile=FILE --prefix=STR 
     --position=bol|eol FILE_1 FILE_2
    
    (batch)
    paired-end_barcode_splitter.py -b --position=bol|eol FILE

DESCRIPTION
    Multiplexed paired end reads extracted with SRAtoolkit's fastq-dump
    (with the option --split-3) are initially separated in two files with
    suffixes *_1.fastq(.gz) (FILE_1) *_2.fastq(.gz) (FILE_2), and the 
    barcodes are only present in the reads from the FILE_1. Therefore,
    this script demultiplexes FILE_1 in the first place, and then sort
    FILE_2 reads according to FILE_1 demultiplexing.
    
    You can run it for one file pair, specifying the FILE_1 and the
    FILE_2 as terminal arguments. FILE_1 is the one that contains the 
    barcoded reads. Alternatively, specifying the option '-b', you can 
    list the input in a batch file that replaces FILE_1 and FILE_2 and 
    will be handled as standard input.
    
    The script requires FASTX-toolkit:
    http://hannonlab.cshl.edu/fastx_toolkit/download.html

OPTIONS (maybe not the best term, most of them are mandatory...)
    -b
        Read a list of organized inputs from the standard input rather 
        than handling a pair of fastq files. The prefix is extracted
        from the file names.
        
        This list is a 3 columns tab separated files organized as:
         column 1: barcode file name
         column 2: mate 1 fastq (barcoded reads)
         column 3: mate 2 fastq
        
        Overrides --bcfile and --prefix values
        
    --bcfile=FILE
        FILE contains the barcode list.
    
    --mismatch=INT
        Number of tolerated mismatches for barcodes.
    
    --prefix=STR
        STR is append to output file names.
        
    --position=bol|eol
        Set up as bol if the adapter is at the 5' end, eol if it is at
        the 3' end. [MANDATORY]
        
    --help
        Display this message

'''

# Require 

import getopt, sys, fileinput, os, subprocess, gzip, shutil
from os import path

class Options(dict):
    '''
    Parse options from the command line arguments and store values.
    '''
    
    def __init__(self, argv):
        
        # set default
        self.set_default()
        
        # handle options with getopt
        try:
            opts, args = getopt.getopt(argv[1:], "b", 
                                       ['help',
                                        'bcfile=',
                                        'mismatch=',
                                        'prefix=',
                                        'position='])
        except getopt.GetoptError, e:
            sys.stderr.write(str(e) + '\n\n' + __doc__)
            sys.exit(1)

        for o, a in opts:
            if o == '--help':
                sys.stdout.write(__doc__)
                sys.exit(0)
            elif o == '-b':
                self['batch'] = True
            elif o == '--bcfile':
                self['bcfile'] = a
            elif o == '--mismatch':
                self['mismatch'] = int(a)
            elif o == '--prefix':
                self['prefix'] = a
            elif o == '--position':
                self['position'] = a
                self['position'].lower()
            else:
                raise getopt.GetoptError("unhandled option: {}".format(o))

        self.args = args
        self.opts = opts
        
    def set_default(self):
    
        # default parameter value
        self['batch'] = False
        self['bcfile'] = None
        self['mismatch'] = 1
        self['prefix'] = None
        self['position'] = None

class BarcodeSplitter(object):
    '''
    Split a multiplexed pair of FASTQ files into corresponding 
    demultiplexed FASTQ files.
    
    Handle the following parameters:
        bcfile      Barcode file 
        mismatch    Number of tolerated mismatches with the barcodes
        prefix      Output prefix 
        position    Take value "eol" (3') or "bol" (5')
    '''
    
    def __init__(self, file_1, file_2, **kwargs): 
        self.args = ["fastx_barcode_splitter.pl", "--suffix", "_demux1.fastq"]
        self.log = { "not matched"       : 0,
                     "demultiplexed"     : 0 }
        self.file_1, self.file_2 = file_1, file_2 
        for k in kwargs:
            
            # barcode file name
            if k == "bcfile":
                self.bcfile = kwargs[k]
                self.args.extend(["--bcfile", self.bcfile])
            
            # number of mismatches
            elif k == "mismatch":
                self.mismatch = kwargs[k]
                self.args.extend(["--mismatch", str(kwargs[k])])
            
            # prefix that will be added to the output
            elif k == "prefix":
                self.prefix = kwargs[k]
                self.args.extend(["--prefix", self.prefix])
                
            # either bol (5') or eol (3')
            elif k == "position":
                self.position = kwargs[k]
                self.args.append("--" + self.position)
            
            # raise an error if the option is not recognized
            else:
                raise ValueError("unknown or unhandled option: {}".format(k))
        
    def run(self):
        
        # demultiplex file_1
        openf1 = gzip.open if self.file_1.endswith(".gz") else open
        with openf1(self.file_1, "rb") as f:
            p = subprocess.Popen(self.args, stdin=subprocess.PIPE)
            shutil.copyfileobj(f, p.stdin)
            p.stdin.close()
        p.wait()
            
        # get output file names
        fnames = [ fname for fname in os.listdir(".") 
                   if fname.startswith(self.prefix) and 
                      fname.endswith("_demux1.fastq") ]
                   
        # read all outputs and the file_2
        openf2 = gzip.open if self.file_2.endswith(".gz") else open
        with openf2(self.file_2, "rb") as f_2:
            
            # to finally close the open file arrays
            try:
            
                # open the mate #1 demultiplexed files for reading
                demultiplexed_files_1 = [ open(fname) 
                                          for fname in fnames ]
                
                # open the mate #2 demultiplexed files for writing
                demultiplexed_files_2 = [ open(fname.replace("_demux1", "_demux2"), "w")
                                          for fname in fnames ]
                                          
                # read 4 lines (1 read info) of each of the mate #1 
                # demultiplexed files
                slots = [ [ f_1.readline().strip() for j in xrange(4) ] 
                          for f_1 in demultiplexed_files_1 ]
                
                # read 4 lines of the file to be demultiplexed
                read_2_lines = [ f_2.readline() for j in xrange(4) ]
                
                # continue while the lines are not null
                while all(read_2_lines):
                    
                    # get the read name before any whitspace
                    read_name = read_2_lines[0].split()[0]

                    # search the read in the demultiplex files, if found,
                    # write read_2_lines to the corresponding mate #2 
                    # demultiplexed file then read the next 4 lines of the 
                    # corresponding mate #1 demultiplexed file.
                    found = False
                    for i in xrange(len(fnames)):
                        read_1_lines = slots[i]
                        
                        # pass when it has reach the end of the file
                        if not read_1_lines[0]: continue
                        
                        if read_1_lines[0].split()[0] == read_name: 
                            demultiplexed_files_2[i].writelines(read_2_lines)
                            slots[i] = [ demultiplexed_files_1[i].readline()
                                         for j in xrange(4) ]
                            self.log["demultiplexed"] += 1
                            found = True
                            break
                    
                    # report unmatched mates                    
                    if not found:
                        sys.stderr.write(
                            "mate missing in file_2: {}\n".format(read_name))
                        self.log["not matched"] += 1
                    
                    # read the next 4 lines of file_2
                    read_2_lines = [ f_2.readline() for j in xrange(4) ]
            finally:
                map(lambda x: x.close(), demultiplexed_files_1)
                map(lambda x: x.close(), demultiplexed_files_2)

def barcode_splitter(file_1, file_2, bcfile, prefix, position, mismatch=1):
    '''
    Run a BarcodeSplitter instance.
    '''
    
    # handle options
    b = BarcodeSplitter(file_1, file_2,
                        bcfile=bcfile, prefix=prefix,
                        position=position, mismatch=mismatch)    
    
    # run the process
    b.run()

def find_prefix(a, b):
    '''
    Return the rightward common part of strings a and b.
    '''
    
    c = ""
    for x, y in zip(a, b):
        if x != y: break
        c += x
    return c

def main(argv=sys.argv):
    '''
    Main function.
    '''
    
    # read options and remove options strings from argv (avoid option 
    # names and arguments to be handled as file names by
    # fileinput.input().
    options = Options(argv)
    sys.argv[1:] = options.args

    # batch input
    if options['batch']:
        for line in fileinput.input():
            if not line.strip(): continue
            bcfile, file_1, file_2 = line.split()
            prefix = find_prefix(file_1, file_2)
            barcode_splitter(file_1, file_2, bcfile, prefix,
                             options['position'], options['mismatch'])
        return 0

    # single input
    else:
        file_1, file_2 = options.args[0], options.args[1]
        barcode_splitter(file_1, file_2, options['bcfile'], options['prefix'],
                         options['position'], options['mismatch'])                

    # return 0 if everything succeeded
    return 0

# does not execute main if the script is imported as a module
if __name__ == '__main__': sys.exit(main())
