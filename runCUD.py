from multiprocessing import freeze_support

from convertuploaddocument import convertuploaddocument

if __name__ == '__main__':
    cud = convertuploaddocument(startfresh=True, cloudflag='aws')
    freeze_support()
    cud.scanforexperiments()
    cud.converttoh5()
    _  = cud.createspreadsheet()
