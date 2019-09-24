from multiprocessing import freeze_support

from convertuploaddocument import *

if __name__ == '__main__':
    osSleep = WindowsInhibitor()
    osSleep.inhibit()
    freeze_support()
    cud = convertuploaddocument(searchpath = '\\\\data.wexac.weizmann.ac.il\\rivlinlab-arc\\raw\\', 
                                WexacH5Path = '\\\\data.wexac.weizmann.ac.il\\rivlinlab-arc\\h5s\\',
                                startfresh=False, cloudflag = "None", localcopyflag = False)
    cud.scanforexperiments()
    cud.parsefields(getvalchar=['=', '>'], nextfieldchar=[',', '_'])
    cud.wexac_copy()
    osSleep.uninhibit()
