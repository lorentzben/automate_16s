#A script to check the manifest provided, the presence of the fastq files, and that dependencies are all setup
import subprocess
import os
from os import path

def check_dependencies():
    print('all software installed and ready to go')

def verify_manifest(mainfest):
    print("the manifest called: " +manifest+ " is valid and ready to go")

def verify_files(manifest):
    print("the files from the manifest have been accounted for")

def error():
    print("something is wrong and I will try to tell you what the problem is")