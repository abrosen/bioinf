import os
import time


QUEUE =  []

def getListOfCurrentFiles():
    return os.listdir('uploads')

def runJob(filename):
    pass




while True:
    time.sleep(10)  # wait for 10 seconds 
    current = getListOfCurrentFiles()
    for x in current:
        if x not in  QUEUE:
            QUEUE.append(x)
    if len(QUEUE):
        task = QUEUE.pop()
        runJob(task)