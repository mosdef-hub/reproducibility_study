"""Delete folders in ws which are empty."""

import os
from glob import glob

os.chdir("workspace")
hashes = glob("./*", recursive=True)

print(len(hashes))

for hash in hashes:
    os.chdir(hash)

    if len(os.listdir()) < 2:
        print(os.listdir())

        os.chdir("..")
        os.system("rm -rf {}".format(hash))
    else:
        os.chdir("..")
