#!/usr/bin/env python

#//////////////////////////////////#
#///   Script to set up homer   ///#
#//////////////////////////////////#
import argparse
import subprocess
import os

parser = argparse.ArgumentParser(description = 'This script checks HOMER configuration')
parser.add_argument("-g", "--genome", help = "Genome which to check. For example 'mm10' .")
args = parser.parse_args()

# Get path of the homer installation
out = subprocess.Popen(['which','homer'], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
stdout,stderr = out.communicate()
ENV_path = stdout.decode().rstrip().replace('bin/homer', 'share')

HOMER_path = subprocess.Popen(
            ['find', ENV_path, '-type', 'd', '-name', 'homer*'], 
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT).communicate()[0].decode().rstrip()

# check if required genome exists and if not,
# Install it using homer perl script

if not os.path.exists(HOMER_path + '/data/genomes/' + args.genome):
  os.system(HOMER_path + '/configureHomer.pl -install ' + args.genome)
else:
    print("\n#############\n\nGenome "
     + args.genome + 
     ' is already installed.\nNothing to be done\nClossing...' + 
     "\n\n#############\n")
