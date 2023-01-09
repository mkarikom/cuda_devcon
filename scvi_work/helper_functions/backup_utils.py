import os

def makepath(mypath):
    if not os.path.exists(mypath): 
        os.makedirs(mypath)