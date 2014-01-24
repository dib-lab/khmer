#! /usr/bin/env python
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2013. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt. Contact: ctb@msu.edu
#

import os

# Enum to list possible File statuses
class FileStatus:
    FileExistsNonEmpty, FileNonExistent, FileEmpty = range(3)

def check_file_status(filePath):
    ''' Check status of file - return if file exists,
    is empty, or neither '''
    if not os.path.exists(filePath):
        return FileStatus.FileNonExistent
    else:
        if os.stat(filePath).st_size==0:
            return FileStatus.FileEmpty
        else:
            return FileStatus.FileExistsNonEmpty
        
def check_space(obj,inFile):
    ''' Estimate size of obj passed, then calculate
    disk space and return extra space needed to store
    file or 0 if sufficient space'''
    
    # Get disk free space in Bytes assuming non SU
    dirPath = os.dirname(inFile)
    target=os.statvfs(dirPath)
    freeSpace = target.f_bsize * target.f_bavail 
    #<TODO>: If SU, use target.f_bfree
    
    # Get input file size as worst case estimate of 
    # output file size
    fileSize = os.stat(inFile).st_size
    
    sizeDiff = fileSize-freeSpace
    if sizeDiff > 0:
        return sizeDiff
    else:
        return 0