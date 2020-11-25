#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
NUM_PROC  = 5
 
for process in range(NUM_PROC):
    pid = os.fork()
    if pid &gt; 0: 
        print("This is the parent process {}".format(os.getpid()))
    else: 
        print("This is the child process {}".format(os.getpid()))) 
        os._exit(0)
print("Parent process is closing")

