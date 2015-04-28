import os
import time 
import threading #upgrade this to multiprocessing 
from Andres_optimizing_v2 import have_a_blast as foo




QUEUE =  []

def getListOfCurrentFiles():
    return os.listdir('uploads')

def dummy(filename):
    print(filename, "dum dum")

#given the file and the first script
def runTheJob(filename, script=foo):
    t = threading.Thread(target= script)#, args = (filename,))
    t.daemon = True
    print "starting", script
    t.start()
    t.join()


def runScheduler():
    while True:
        time.sleep(10)  # wait for 10 seconds 
        current = getListOfCurrentFiles()
        for x in current:
            if x not in  QUEUE:
                QUEUE.append(x)
        if len(QUEUE):
            task = QUEUE.pop()
            runTheJob(task,foo)

def doTheThing():
    runTheJob("andres",foo)

if __name__ == '__main__':
    doTheThing()
