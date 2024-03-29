#!/usr/bin/env python3

# DISCLAIMER
# the code below was taken from stackoverflow.com/questions/489861/locking-a-file-in-python
# by 'Thomas Lux'; no license stated; accessed 23-12-2022, at 15:00

import os

try:
    # Posix based file locking (Linux, Ubuntu, MacOS, etc.)
    #   Only allows locking on writable files, might cause
    #   strange results for reading.
    import fcntl

    def lock_file(f):
        if f.writable():
            fcntl.lockf(f, fcntl.LOCK_EX)

    def unlock_file(f):
        if f.writable():
            fcntl.lockf(f, fcntl.LOCK_UN)

except ModuleNotFoundError:
    # Windows file locking
    import msvcrt

    def file_size(f):
        return os.path.getsize( os.path.realpath(f.name))

    def lock_file(f):
        msvcrt.locking(f.fileno(), msvcrt.LK_RLCK, file_size(f))

    def unlock_file(f):
        msvcrt.locking(f.fileno(), msvcrt.LK_UNLCK, file_size(f))


# Class for ensuring that all file operations are atomic, treat
# initialization like a standard call to 'open' that happens to be atomic.
# This file opener *must* be used in a "with" block.
class AtomicOpen:
    # Open the file with arguments provided by user. Then acquire
    # a lock on that file object (WARNING: Advisory locking).
    def __init__(self, path, *args, **kwargs):
        # Open the file and acquire a lock on the file before operating
        self.file = open(path, *args, **kwargs)
        # Lock the opened file
        lock_file(self.file)

    # Return the opened file object (knowing a lock has been obtained).
    def __enter__(self, *args, **kwargs): return self.file

    # Unlock the file and close the file object.
    def __exit__(self, exc_type=None, exc_value=None, traceback=None):
        # Flush to make sure all buffered contents are written to file.
        self.file.flush()
        os.fsync(self.file.fileno())
        # Release the lock on the file.
        unlock_file(self.file)
        self.file.close()
        # Handle exceptions that may have come up during execution, by
        # default any exceptions are raised to the user.
        if exc_type is not None:
            return False
        else:
            return True

# Now, AtomicOpen can be used in a with block where one would normally use an open statement.
#
# WARNINGS:
#     If running on Windows and Python crashes before exit is called, I'm not sure what the lock behavior would be.
#     The locking provided here is advisory, not absolute. All potentially competing processes must use the
#     "AtomicOpen" class.
#     As of (Nov 9th, 2020) this code only locks writable files on Posix systems. At some point
#     after the posting and before this date, it became illegal to use the fcntl.lock on read-only files.
