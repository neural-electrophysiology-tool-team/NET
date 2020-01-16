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
    
    #Fills out entire spreadsheet again (work around until I rewrite excel creation to append)
#     cud = convertuploaddocument(searchpath = '\\\\data.wexac.weizmann.ac.il\\rivlinlab-arc\\raw\\', 
#                             WexacH5Path = '\\\\data.wexac.weizmann.ac.il\\rivlinlab-arc\\h5s\\',
#                             startfresh=True, cloudflag = "None", localcopyflag = False)
#     cud.scanforexperiments()
#     cud.parsefields(getvalchar=['=', '>'], nextfieldchar=[',', '_'])
#     cud.field.to_excel(os.path.join("\\\\data.wexac.weizmann.ac.il\\rivlinlab-arc\\h5s\\", 'experiments.xlsx'), index=False)
    osSleep.uninhibit()
