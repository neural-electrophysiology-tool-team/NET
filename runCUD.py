from multiprocessing import freeze_support

from convertuploaddocument import *

if __name__ == '__main__':
    osSleep = WindowsInhibitor()
    osSleep.inhibit()

    cud = convertuploaddocument(startfresh=True, cloudflag='aws',localcopyflag = True, clearlogflag = True, overwriteh5 = True)
    freeze_support()
#     cud.scanforexperiments()
#     cud.converttoh5()
#     _  = cud.createspreadsheet()
    cud = convertuploaddocument(searchpath = '\\\\data.wexac.weizmann.ac.il\\rivlinlab-arc\\raw\\', 
                                WexacH5Path = '\\\\data.wexac.weizmann.ac.il\\rivlinlab-arc\\h5s\\',
                                startfresh=True, cloudflag = "None", localcopyflag = False)
    cud.scanforexperiments()
    cud.parsefields(getvalchar=['=', '>'], nextfieldchar=[',', '_'])
    cud.wexac_copy()
    osSleep.uninhibit()
