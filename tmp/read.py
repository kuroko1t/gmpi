#!/usr/bin/env python3

import re

def readapi(path):
    lines = []
    with open(path) as f:
        lines = f.readlines()

    api_key0 = "OMPI_DECLSPEC  .* (MPI_\S*)\("
    api_hash = {}
    env = True
    api_line = []
    apiname = ""
    for line in lines:
        m = re.search(api_key0,line)
        if m:
            apiname = m.group(1)
            if re.search(";",line):
                api_hash[apiname] = [line]
            else:
                env = True
                api_line.append(line)
        elif env:
            api_line.append(line)
            if re.search(";",line):
                api_hash[apiname] = api_line
                api_line = []
                env = False
    print(api_hash)

readapi("/usr/local/include/mpi.h")
