This is a local copy of the [AndorSIF repository](https://github.com/RJ3/AndorSIF) placed here for convenience.

# AndorSIF
This is a repository for certain files which allow you to compile a mex file to load the proprietary Andor SIF format into MATLAB

How to compile .mexa64 file on linux and use for sifopen.
Rafael Jaimes
Last updated: April 21st, 2014.

********PART 1: DOWNLOADING AND INSTALLING THE ANDOR LINUX SDK**************  
Obtain the proprietary .so (or .dll for windows) files from the Andor website  
https://www.andor.com/download/  
Software -> Andor SDK  

Andor Linux SDK - 2.96.30004.0  
Last Updated: 20/12/2013 11:30:38  
Size: 13.55 MB  
Filename: andor-2.96.30004.0.tar.gz  
Andor SDK for CCD cameras. Designed to work on generic Linux distributions.   
To install simply extract the archive and follow the instructions in the README contained within.   

tar -zxvf andor-2.96.30004.0.tar.gz  
cd andor  
sudo ./install_andor  
5 (All USB Cameras)  

This will install SDK header files in the following directory:  
/usr/local/include  
And also static libraries in the following directory:  
/usr/local/lib  
For uninstall:  
/etc/andor/andor.install  

********PART 2: COMPILING THE .C FILE TO CREATE .MEXA64*****************  
Create a directory and place in the following files:  
atsifiomex.c  
ATSIFIO.h  
ATSIFTypes.h  
ATSIFProperties.h  
ATSIFErrorCodes.h  
ATPrimitiveTypes.h  
ATLibraryExport.h  

Open matlab, then cd to the directory with the files.  
Run the following command:  
mex atsifiomex.c /usr/local/lib/libatsifio.so  
Which creates:  
atsifiomex.mexa64  

**********PART 3: USING THE .MEXA64 AND SIFOPEN**************  
Make sure you install the Andor Linux SDK on every machine you plan on using  
You don't have to compile the .mexa64 file every time, most of the time the same file will  
work on different machines.  
Just put the atsifio.mexa64 and sifopen99.m in the same folder, and run from there.  
You should be able to open any andor .SIF files directly with matlab.

If you get an error about missing libatsifio.so.2 then you need to add this to your console
before running matlab. This will point matlab to the shared objects.

export LD_LIBRARY_PATH=/usr/local/lib/
