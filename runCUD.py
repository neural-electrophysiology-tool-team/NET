from multiprocessing import freeze_support

from convertuploaddocument import *

if __name__ == '__main__':
    osSleep = WindowsInhibitor()
    osSleep.inhibit()

    cud = convertuploaddocument(startfresh=True, cloudflag='aws',localcopyflag = True, clearlogflag = True, overwriteh5 = True)
    freeze_support()
    cud.scanforexperiments()
    cud.converttoh5()
    _  = cud.createspreadsheet()
    
    osSleep.uninhibit()
