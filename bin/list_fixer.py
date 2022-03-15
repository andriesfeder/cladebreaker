#!/usr/bin/env python

import os
import sys
import gzip
import argparse

def parse_args(args=None):
    Description = "Reformat nf-core/whatsgnu topgenomes samplesheet file"
    Epilog = "Example usage: python list_fixer.py <FILE_IN> <FILE_OUT>"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("FILE_IN", help="Input topgenomes text file.")
    parser.add_argument("FILE_OUT", help="Output file.")
    return parser.parse_args(args)

def clean_list(file_in, file_out):
    gca = []
    file=open(file_in, "r")
    newFile=open(file_out,"w")
    file.readline()
    for raw_num in file:
        raw_num = raw_num.split("GCA")
        gca_num = raw_num[1]
        gca_num = gca_num.split("_"[0])
        final_num = "GCA_"+gca_num[1]
        # if final_num in gca:
        #    continue
        #else:
        if len(final_num) > 6:
            newFile.write(final_num+"\n")
            gca.append(final_num)

    newFile.close()
    file.close()

def main(args=None):
    args = parse_args(args)
    clean_list(args.FILE_IN, args.FILE_OUT)

if __name__ == "__main__":
    sys.exit(main())
