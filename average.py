from contextlib import contextmanager
from itertools import izip
from glob import iglob
from math import sqrt
from sys import exit

@contextmanager
def multi_file_manager(files, mode='rt'):
    files = [open(file, mode) for file in files]
    yield files
    for file in files:
        file.close()

# generator function to read, convert, and yield each value from a text file
def read_values(file, datatype=float):
    for line in file:
        for value in (datatype(word) for word in line.split()):
            yield value

# enumerate multiple egual length iterables simultaneously as (i, n0, n1, ...)
def multi_enumerate(*iterables, **kwds):
    start = kwds.get('start', 0)
    return ((n,)+t for n, t in enumerate(izip(*iterables), start))

DATA_FILE_PATTERN = '*.dat'
MIN_DATA_FILES = 2

with multi_file_manager(iglob(DATA_FILE_PATTERN)) as datfiles:
    num_files = len(datfiles)
    if num_files < MIN_DATA_FILES:
        print('Less than {} .dat files were found to process, '
              'terminating.'.format(MIN_DATA_FILES))
        exit(1)

    # determine number of rows and cols from first file
    temp = [line.split() for line in datfiles[0]]
    num_rows = len(temp)
    num_cols = len(temp[0])
    datfiles[0].seek(0)  # rewind first file
    del temp  # no longer needed
    print '{} .dat files found, each must have {} rows x {} cols\n'.format(
           num_files, num_rows, num_cols)

    means = []
    std_devs = []
    divisor = float(num_files-1)  # Bessel's correction for sample standard dev
    generators = [read_values(file) for file in datfiles]
    for _ in xrange(num_rows):  # main processing loop
        for _ in xrange(num_cols):
            # create a sequence of next cell values from each file
            values = tuple(next(g) for g in generators)
            mean = float(sum(values)) / num_files
            means.append(mean)
            means_diff_sq = ((value-mean)**2 for value in values)
            std_dev = sqrt(sum(means_diff_sq) / divisor)
            std_devs.append(std_dev)

print 'Average and (standard deviation) of values:'
with open('means.dat', 'wt') as averages:
    for i, mean, std_dev in multi_enumerate(means, std_devs):
        print '{:.9f} ({:.9f})'.format(mean, std_dev),
        averages.write('{:.9f}'.format(mean))  # note std dev not written
        if i % num_cols != num_cols-1:  # not last column?
             averages.write(' ')  # delimiter between values on line
        else:
            print  # newline
            averages.write('\n')
