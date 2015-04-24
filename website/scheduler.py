import os
import time 
import threading #upgrade this to multiprocessing





QUEUE =  []

def getListOfCurrentFiles():
    return os.listdir('uploads')

def dummy(filename):
    print(filename, "dum dum")

#given the file and the first script
def runTheJob(filename, script=dummy):
    t = threading.Thread(target= script, args = (filename,))
    t.daemon = True
    try:
        t.start()
    except e:
        print("Something went wrong and someone isn't getting their data.", e)



def runScheduler():
    while True:
        time.sleep(10)  # wait for 10 seconds 
        current = getListOfCurrentFiles()
        for x in current:
            if x not in  QUEUE:
                QUEUE.append(x)
        if len(QUEUE):
            task = QUEUE.pop()
            runTheJob(task)

def doTheThing():
    runTheJob("someName")

if __name__ == '__main__':
    doTheThing()