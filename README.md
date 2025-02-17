## Installation 
This software requires the CMIG core utilities, which should be cloned alongside this repo:
```
git clone --depth 1 https://github.com/cmig-research-group/cmig_core_utils.git && \
git clone --depth 1 https://github.com/cmig-research-group/rsi_pelvis.git
```

after checking this out, cd to the rsi_prostate dir and run `make install` to copy files into that data directory.  

Then you can run `make` to build the program.  After that you can test it with some test files with the `make run` command.

You will need to modify the makefile to point to the test files you have put in the data directory.  
