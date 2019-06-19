import _pickle as pickle
import logging

import numpy as np
import pandas as pd

import concurrent.futures as cf

import os
import subprocess

import re
import hashlib

import boto3
from boto3.s3.transfer import TransferConfig
from google.cloud import storage


class convertuploaddocument:
    """Class which allows for upload and batch conversion
    of mcds files generated by MultiChannel Systems Experimenter
    The class is meant to work on the MCS recording computer in Windows 10.
    Version 1: 10-06-19 Noah Dolev
    Version 1.1: 11-06-19 Noah Dolev [added spreadsheet export]
    """

    def walk(self):
        """Get all files in path
        """
        for p, d, f in os.walk(self.searchpath):
            for ff in f:
                if ff.endswith(self.suffix):
                    self.files.append(os.path.join(p, ff))

    def scanforexperiments(self):
        """Takes a search path and looks recursively for all MSRD files
        """

        self.logger.info('Start: Scan for new files')
        self.walk()
        self.logger.info('Process: %d _total_ files found' % len(self.files))
        oldfilelist = []

        # Checks if this was run in the past, if so, only processes the new files
        oldfileslog = os.path.join(self.logpath, 'filelist.p')
        if self.startfresh:
            self.logger.info('Start: Erasing File History')
            pickle_out = open(oldfileslog, "wb")
            pickle.dump([], pickle_out)
            pickle_out.close()

        if os.path.getsize(oldfileslog) > 0:
            pickle_in = open(oldfileslog, "rb")
            oldfilelist = pickle.load(pickle_in)
            pickle_in.close()
            self.logger.info('Process: %d _old_ files found' % len(oldfilelist))

        pickle_out = open(oldfileslog, "wb")
        pickle.dump(self.files, pickle_out)
        pickle_out.close()
        self.logger.info('Process: Old files list updated')

        if len(oldfilelist) != 0:
            self.files = list(set(self.files) - set(oldfilelist))
        self.numfiles = len(self.files)
        self.logger.info('Process: %d _new_ files found' % self.numfiles)
        self.logger.info('End: Scan for files')

    def parsefields(self, getvalchar='<', nextfieldchar='_'):
        self.field['fullpathmsrd'] = self.files * 2
        self.field.sort_values(by=['fullpathmsrd'], inplace=True)
        self.field.reset_index(drop=True, inplace=True)
        self.field['MEAfiles'] = self.field['fullpathmsrd'].apply(
            lambda x: x.split('\\')[-1].replace(self.suffix, 'h5'))
        self.field['MEAfiles'].loc[::2] = self.field['MEAfiles'].loc[::2].apply(lambda x: x.replace('h5', 'bin'))
        self.field['folder'] = ['\\' + os.path.join(*word[:-1]) for word in
                                [f.split('\\') for f in self.field['fullpathmsrd']]]
        self.field['recNames'] = self.field['MEAfiles'].apply(lambda x: hashlib.md5(x.encode()).hexdigest())

        def setrecformat(x):
            if 'h5' in x:
                return ('MCH5Recording')
            else:
                return ('binaryRecording')

        self.field['recFormat'] = [setrecformat(x) for x in self.field['MEAfiles']]

        flatten = lambda l: [item for sublist in l for item in sublist]
        nvc = '|'.join(flatten([nextfieldchar]))
        gvc = '|'.join(flatten([getvalchar]))
        for ind in self.field.index:
            word = self.field['fullpathmsrd'].loc[ind].split('\\')
            keyvalpair = [re.split(nvc, w) for w in word]
            for kv in keyvalpair:
                for k in kv:
                    if '=' in k:
                        temp = re.split(gvc, k)[0]
                        temp = temp[-np.min([6, len(temp)]):]  # maximum field length is 6 characters
                        key = ''.join(j for j in temp if not j.isdigit())  #
                        temp = re.split(gvc, k)[1:]
                        val = ''.join(temp)
                        if key not in self.field.columns:
                            self.field[key] = 'unspecified'
                        self.field[key].loc[ind] = val.split('.' + self.suffix)[0]
        self.field.fillna('unspecified', inplace=True)

    def makeGSdir(self, directory='database/'):
        """Creates a google storage bucket for data
        """
        blob = self.bucket.blob(directory.replace('\\', '/'))
        blob.upload_from_string('', content_type='application/x-www-form-urlencoded;charset=UTF-8')

    def createDirStructure(self):
        """Parses pandas dataframe into a directory structure for easy reading by a database engine
        whereby each field and value are changed into "x=y" directories.
        _Should add code to automatically form the directory hierarchy so that the smallest number of files is in the lowest folder_
        """
        table = self.field
        table = table[table['MEAfiles'].str.contains('h5')].drop(
            ['fullpathmsrd', 'exclude', 'recFormat', 'recNames', 'folder', 'comments'], axis=1)
        table = table.reindex(columns=(['user'] + list([col for col in table.columns if col != 'user'])))
        dirs = [t for t in table.columns.values if t != 'MEAfiles']
        for i in range(0, table.shape[0]):
            dtemp = 'database'
            for d in dirs:
                if 'unspecified' not in table[d].iloc[i]:
                    dtemp = os.path.join(dtemp, d + '=%s' % table[d].iloc[i])
            self.directory.append(dtemp.replace('\\', '/'))

    def uploadFile(self, file, path, bucket):
        """Uploads files to Google Cloud Storage
        :param file: the file to be uploaded
        :param path: the generated path using the file's inferred experiment attributes
        :param bucket: the bucket to place the file
        """
        path = path + '\\'
        path = path.replace('\\', '/')
        testexists = bucket.blob(path + file.replace('.msrd', '.h5').split('\\')[-1])
        if not testexists.exists():
            self.makeGSdir(directory=path)
            bashCommand = 'gsutil -o GSUtil:parallel_composite_upload_threshold=150M -m cp "%s" "gs://%s"' % (
                file.replace('.msrd', '.h5'), path)
            # Multi-threaded composite upload. Still quite slow.
            subprocess.Popen(bashCommand, stdout=subprocess.PIPE, shell=True)

    def multi_part_upload_with_s3(self, file_path, key_path, bucket="meadata"):
        """Uploads files in parallel to Amazon s3 buckets
        :param file_path: path of file to upload
        :param key_path: generated path using the file's inferred experiment attributes
        :param bucket: bucket to place the file
        """
        # Multipart upload
        with open(self.awsaccesskey, 'rb') as csvfile:
            keys = pd.read_csv(csvfile)
        s3 = boto3.resource('s3', aws_access_key_id=keys["Access key ID"][0],
                            aws_secret_access_key=keys["Secret access key"][0])
        config = TransferConfig(multipart_threshold=1024 * 25, max_concurrency=8,
                                multipart_chunksize=1024 * 25, use_threads=True)
        b = s3.Bucket(bucket)
        objs = list(b.objects.filter(Prefix=key_path))
        if len(objs) > 0 and objs[0].key == key_path:
            print("%s already Exists - skipping upload!" % file_path.split('\\')[-1])
        else:
            s3.meta.client.upload_file(Filename=file_path, Bucket=bucket, Key=key_path,
                                       Config=config)

    def localnetworkcopy(self, f, d, archive=True):
        """
        :param f: file to copy
        :param d: directory to place file
        :param archive: Flag for whether to also copy file to local network archive
        """
        self.logger.info('Process: Copy h5 to local network location')
        subprocess.call(
            ["xcopy", f.replace("msrd", "h5"), os.path.join(self.localsharepath, d.replace('/', '\\')) + '\\',
             "/K/O/X/i"])
        self.logger.info('Process: %s copied succesfully' % f.replace("msrd", "h5"))
        if archive:
            self.logger.info('Process: Copy MSRD and MSRS to local network archive')
            subprocess.call(
                ["xcopy", f, os.path.join(self.localarchivepath, d.replace('/', '\\')) + '\\', "/K/O/X/i/Y"])
            subprocess.call(
                ["xcopy", f.replace("msrd", "msrs"), os.path.join(self.localarchivepath, d.replace('/', '\\')) + '\\',
                 "/K/O/X/i/Y"])
            self.logger.info('Process: %s msrd and msrs copied succesfully to archive' % f.split('.')[0])

    def createspreadsheet(self):
        """Takes parsed file dataframe and creates a spreadsheet based on the fields
            :return pandas dataframe with list of uploaded files and their metadata
        """
        self.field = self.field.replace("unspecified", "")
        spreadsheetpath = os.path.join(self.searchpath, 'experiments.xlsx')
        self.field.to_excel(spreadsheetpath, index=False)
        if self.cloudflag == "gcs":
            self.logger.info('Process: Starting upload of spreadsheet to GS')
            self.uploadFile(file=spreadsheetpath, path='/database',
                            bucket=self.bucket)
            self.logger.info('Process: Successfully uploaded spreadsheet to GS')
        elif self.cloudflag == "aws":
            with open(self.awsaccesskey, 'rb') as csvfile:
                keys = pd.read_csv(csvfile)
            self.multi_part_upload_with_s3(file_path=spreadsheetpath, key_path="database/experiments.xlsx")
            self.logger.info("Process: Spreadsheet uploaded to AWS s3")
        else:
            self.logger.info('Process: Not uploading spreadsheet to cloud storage')
        if self.localcopyflag == True:
            self.localnetworkcopy(spreadsheetpath, self.localsharepath)
        return (self.field)

    def converttoh5(self):
        """Takes files that were discovered and converts them to h5 then uploads them to google cloud / amazon S3 and/or makes a local network copy
        """

        def processfiles(f, d, bucket=self.bucket, mcspath=self.mcspath, cloudflag=self.cloudflag,
                         localcopyflag=self.localcopyflag):
            """Internal function for multiprocessing conversion and upload
            :param f: file to process
            :param d: generated path by createDirStructure
            :param bucket: bucket to place file
            :param mcspath: path to MCS commandline conversion tool
            :param cloudflag: flag whether to use google ("gcs", "aws" or "local")
            :param localcopyflag: if True, also copies file to a local network direction
            """

            if not os.path.isfile(f.replace('.msrd', '.h5')):
                bashCommand = '%s -t hdf5 "%s"' % (mcspath, f.replace('.msrd', '.msrs'))
                process = subprocess.Popen(bashCommand, stdout=subprocess.PIPE)
                output, error = process.communicate()
                if error is not None:
                    self.logger.error('File failed to convert with error: \n %s' % error)
                else:
                    self.logger.info('Process: Successfully converted file to H5')
                    # Workaround since there is no way to specify output file name with MCS commandline tool
                    os.rename(f.replace('.msrd', '.h5'), os.path.join('//'.join(f.split('//')[:-1]), hashlib.md5(
                        f.split('//')[-1].split('.')[0].encode()).hexdigest() + '.h5'))
                    self.logger.info('Process: File renamed to md5 hash successfully')
                    f = os.path.join('//'.join(f.split('//')[:-1]),
                                     hashlib.md5(f.split('//')[-1].split('.')[0].encode()).hexdigest() + '.h5')

            if (cloudflag == "gcs"):
                self.logger.info('Process: Starting upload of HDF5 to target directory of GS')
                self.uploadFile(file=f, path=d,
                                bucket=bucket)
                # can be improved by composing a composite object containing all the files to upload
                self.logger.info('Process: Successfully uploaded H5 to GS')
            elif (cloudflag == "aws"):
                with open("D:\\code\\user=ND\\ND_AccessKey.csv", 'rb') as csvfile:
                    keys = pd.read_csv(csvfile)
                file_path = f.replace('.msrd', '.h5')
                key_path = d + '/' + f.replace('.msrd', '.h5').split('\\')[-1]
                self.logger.info('Process: %s Uploading' % file_path)
                self.logger.info('Process: File uploading to: %s' % key_path)
                self.multi_part_upload_with_s3(file_path=file_path, key_path=key_path, bucket="meadata")
            else:
                self.logger.info('Process: Not uploading to cloud storage')
            if (localcopyflag == True):
                self.localnetworkcopy(f, d)

        self.logger.info('Start: Conversion from MSDS to H5')
        self.parsefields(getvalchar=['=', '>'], nextfieldchar=[',', '_'])
        self.createDirStructure()
        self.logger.info('Process: Created target directory structure')
        with cf.ProcessPoolExecutor() as executor:
            _ = [executor.submit(processfiles(f=self.files[i], d=self.directory[i])) for i in range(0, len(self.files))]
        self.logger.info('End: Conversion from MSDS to H5')

    def __init__(self, searchpath="D:\\Multi Channel DataManager\\", startfresh=False,
                 suffix='.msrd', gcs_credentials_path="D:\\code\\user=ND\\divine-builder-142611-9884de65797a.json",
                 gcs_project_id='divine-builder-142611', bucketname='meadata',
                 logpath=os.path.join('d:\\', 'code', 'user=ND'),
                 mcspath=os.path.join('d:\\', 'code', 'user=ND', 'McsDataCommandLineConverter',
                                      "McsDataCommandLineConverter.exe"),
                 cloudflag="aws", localcopyflag=True, localsharepath="\\\\132.77.73.171\\MEA_DATA\\",
                 localarchivepath="\\\\data.wexac.weizmann.ac.il\\rivlinlab-arc\\",
                 awsaccesskey="D:\\code\\user=ND\\ND_AccessKey.csv"):
        """Initialize class
        :param searchpath: path to directory with data files
        :param startfresh: flag, if 1 then attempt to upload all files
        :param suffix: suffix of data files
        :param gcs_credentials_path: path to google cloud credential file
        :param gcs_project_id: name of google cloud project
        :param bucketname: name of bucket to upload files
        :param logpath: path to save log file
        :param mcspath: path to mcs data command line converter
        :param cloudflag: flag, if 1 then attempt to upload files to google
        :param localsharepath: path to local network shared directory
        """
        os.environ['GOOGLE_APPLICATION_CREDENTIALS'] = gcs_credentials_path
        self.gcs_client = storage.Client(project=gcs_project_id)
        self.bucket = self.gcs_client.get_bucket(bucketname)

        self.logpath = logpath
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.INFO)

        # create a file handler
        self.handler = logging.FileHandler(os.path.join(self.logpath, 'SpreadSheet_RunLog.log'), mode='a')
        self.handler.setLevel(logging.INFO)

        # create a logging format
        self.formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        self.handler.setFormatter(self.formatter)

        # add the handlers to the logger
        self.logger.addHandler(self.handler)
        self.files = []
        colnames = ['fullpathmsrd', 'exclude', 'folder', 'coating', 'cleaning', 'MEAfiles', 'recFormat', 'recNames',
                    'comments']
        self.field = pd.DataFrame(columns=colnames)
        self.searchpath = searchpath
        self.startfresh = startfresh
        self.suffix = suffix
        self.numfiles = 0
        self.directory = []
        self.mcspath = mcspath
        self.cloudflag = cloudflag
        self.localcopyflag = localcopyflag
        self.localsharepath = localsharepath
        self.localarchivepath = localarchivepath
        self.awsaccesskey = awsaccesskey