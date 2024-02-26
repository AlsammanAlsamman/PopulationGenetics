#!/usr/bin/python

############################################## About Author #########################################
# Created by: Alsamman M. Alsamman                                                                  #
# Emails: smahmoud [at] ageri.sci.eg or A.Alsamman [at] cgiar.org or SammanMohammed [at] gmail.com  #
# License: MIT License - https://opensource.org/licenses/MIT                                        #
# Disclaimer: The script comes with no warranty, use at your own risk                               #
# This script is not intended for commercial use                                                    #
#####################################################################################################

import sys
import threading
import time
import os

# performing some tasks using threading

def task1(firstInput):
    print("Task 1 is running")
    # run the bash command
    os.system("sh command.sh")
    time.sleep(2)
    print("Task 1 is done")


# creating 10 threads
threads = []
targetFolder = "path/to/folder"
files2process = os.listdir(targetFolder)
# run 10 each time to avoid memory issues
for i in range(0, len(files2process), 10):
    for j in range(i, i+10):
        if j < len(files2process):
            t = threading.Thread(target=task1, args=(files2process[j],))
            threads.append(t)
            t.start()
    for t in threads:
        t.join()
    threads = []