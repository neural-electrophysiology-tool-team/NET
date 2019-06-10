import os
import pandas as pd
import xlsxwriter
import logging
import subprocess
from google.cloud import storage
import os
import _pickle as pickle
import numpy as np
import concurrent.futures as cf
import threading
import boto3
import os
import sys
import csv
import re
import hashlib
from boto3.s3.transfer import TransferConfig
from convertuploaddocument import convertuploaddocument


cud = convertuploaddocument(startfresh=True, gcs = 0)
cud.scanforexperiments()
cud.converttoh5()
