"""
Greg McInes
Altman Lab
gmcinnes@stanford.edu

Some custom commands I find useful
"""

# Detect if a file is gzipped before opening
def g_open(file):
    if file.endswith('gz'):
        import gzip
        return gzip.open(file)
    return open(file)

#def g_line(file_object):


# Check if a function exists on the os
# Source: https://stackoverflow.com/a/377028/3570907
def which(program):
    import os
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None

# Run a system command and return the output from stdout
# Optionally, check if it succeeded (exit code = 0)
def run_system_cmd(command, check_exit=True):
    import subprocess
    #cmd_list = command.split()
    p = subprocess.Popen(command, stdout=subprocess.PIPE, stderr = subprocess.PIPE, shell=True)
    out, err = p.communicate()
    if check_exit:
        exit_code = p.returncode
        if exit_code != 0:
            print("Process failed!\nCOMMAND: %s\nEXIT CODE: %s\nSTDERR: %s" % (command, exit_code, err))
            exit(1)
    return out

def check_file_paths(paths, terminate=True):
    import os
    import sys
    for p in paths:
        if not os.path.exists(p):
            print("Invalid file path: %s" % p)#, file=sys.stderr) # shits broken and i don't know why
            #if terminate:
            #    print("Exiting", file=sys.stderr)
            #    exit(1)

def timestamp(sep="", pretty=False):
    from datetime import datetime
    if pretty:
        return datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    return datetime.now().strftime(sep.join(["%Y", "%m", "%d", "%H", "%M"]))

def glog(message, type, end=False):
    import sys
    print("%s %s: %s" % (timestamp(pretty=True), type, message))#, file=sys.stderr)
    if end:
        print("Exiting")#, file=sys.stderr)
        exit()

def intersection(a, b):
    # Return the intersection of both lists
    return list(set(a) & set(b))

def difference(list_a, list_b):
    # Return a list of the elements in a that are not in b
    return list(set(list_a) - set(list_b))

def mean(numbers):
    return float(sum(numbers)) / max(len(numbers), 1)

def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i + n]

def byte_decoder(a):
    #print(a)
    return a.decode("utf-8")
