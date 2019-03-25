/* Include for Andor descriptor definition */
/*
#ifdef _WIN32
 #include <windows.h>            // required for all Windows applications
 #include <stdio.h>
 #define SPRF sprintf_s
#else
 /* For UNIX */
 /*
 #include <unistd.h>
 #include <stdlib.h>
 #define SPRF snprintf
#endif
*/

#define SPRF snprintf
#include "ATSIFIO.h"  
#include "ATSIFTypes.h"
#include "mex.h"
#include "matrix.h"
#include <string.h>

#define STRLEN 256

/* Function to check incoming parameters */
void AndorParCheck(char *func, int Anlhs, int Anrhs, char *Snrhs,
			       int nlhs, mxArray *plhs[],int nrhs,const mxArray *prhs[]) {
  int i=0;
  char err[STRLEN];
  char *DS="SDDDDDDDDDDDDDDDDDDDD";
/*Check incoming parameters */
  *err = 0;			/* Set err string empty */
  if(Anrhs > nrhs) {
    SPRF(err,STRLEN,"%s needs %d input parameters",func,Anrhs);
  }
  if(Snrhs == NULL) Snrhs = DS;		/* Default test */
  if(Anrhs > 1)
  while(++i < Anrhs) {
    switch(Snrhs[i]) {
      case 'I':
        if(!mxIsInt32(prhs[i]))
          SPRF(err,STRLEN,"%s- par %d not Integer",func, i+1);
        break;
	  case 'F':
        if(!mxIsSingle(prhs[i]))
		  SPRF(err,STRLEN,"%s- par %d not float",func, i+1);
		break;
      case 'S':
        if(!mxIsChar(prhs[i]))
          SPRF(err,STRLEN,"%s- par %d not String",func, i+1);
        break;
      case 'D':
        if(!mxIsDouble(prhs[i]))
        SPRF(err,STRLEN,"%s- par %d not Double",func, i+1);
		break;
      default:
        SPRF(err,STRLEN,"AndorParCheck [%s] check %d unknown",func,i);
	}
  }
/* Look for O/P parameters */
  if(Anlhs > nlhs)
	SPRF(err,STRLEN,"%s needs %d output parameters",func,Anlhs);
/* If anything written, ABORT */
  if(*err)
	mexErrMsgTxt(err);
}
/* Exit Function When the MEX-file is cleared, must also free permanent memory */ 
void ExitFcn(void)
{
}
/* Matlab Function ================================================================================*/
void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
  char func[STRLEN];
  int  dims1[]={1,0};
  int  ret=999;		/* Return Value */

  int opt = mxIsChar(prhs[0]) ? -1 : (int)*mxGetPr(prhs[0]);	       /* mode */

  *func='\0';
  if(opt < 0) mxGetString(prhs[0],func,STRLEN-1);   /* get conversion string */

/* unsigned int ATSIF_SetFileAccessMode(ATSIF_ReadMode _mode); ===============================================================*/
  if(!strcmp(func,"ATSIF_SetFileAccessMode")) {
    AT_32 i_mode;
    AndorParCheck(func, 1, 2, NULL, nlhs, plhs, nrhs, prhs);
    i_mode = (int)*mxGetPr(prhs[1]);
    i_mode += 0x40000000; 
    ret = ATSIF_SetFileAccessMode((ATSIF_ReadMode)(i_mode));
  }
/* unsigned int ATSIF_ReadFromFile(AT_C * _sz_filename); =====================================================================*/
  else if(!strcmp(func,"ATSIF_ReadFromFile")) {
    char sz_fileName[STRLEN];
    AndorParCheck(func, 1, 2, "SS", nlhs, plhs, nrhs, prhs);
    mxGetString(prhs[1], sz_fileName, STRLEN-1);
    ret = ATSIF_ReadFromFile(sz_fileName);
  }
/* unsigned int ATSIF_CloseFile(); ======================================================================*/
  else if(!strcmp(func,"ATSIF_CloseFile")) {
    ret = ATSIF_CloseFile();
  }
/* unsigned int ATSIF_ReadFromByteArray(AT_U8 * _buffer, AT_U32 _ui_bufferSize); =======================================================================*/
  else if(!strcmp(func,"ATSIF_ReadFromByteArray")) {
    AT_U32 ui_bufferSize;
    AndorParCheck(func, 1, 3, NULL, nlhs, plhs, nrhs, prhs);
    ui_bufferSize = (AT_32)*mxGetPr(prhs[2]);
    ret = ATSIF_ReadFromByteArray((unsigned char *)mxGetPr(prhs[1]),ui_bufferSize);
  }
/* unsigned int ATSIF_IsLoaded(AT_32 * _i_loaded); ======================================================================*/
  else if(!strcmp(func,"ATSIF_IsLoaded")) {
    AT_32 i_loaded;
    AndorParCheck(func, 2, 1, NULL, nlhs, plhs, nrhs, prhs);	  
    ret = ATSIF_IsLoaded(&i_loaded);
    plhs[1] = mxCreateDoubleMatrix(1,1, mxREAL);
    *mxGetPr(plhs[1]) = (double) i_loaded ;
  }
/* unsigned int ATSIF_IsDataSourcePresent(ATSIF_DataSource _source, AT_32 *_i_present); ======================================================*/
  else if(!strcmp(func,"ATSIF_IsDataSourcePresent")) {
    AT_32 i_present;
    AT_32 i_source;
    AndorParCheck(func, 2, 2, NULL, nlhs, plhs, nrhs, prhs);	 
    i_source = (int)*mxGetPr(prhs[1]);
    i_source += 0x40000000;  
    ret = ATSIF_IsDataSourcePresent((ATSIF_DataSource)(i_source),&i_present);
    plhs[1] = mxCreateDoubleMatrix(1,1, mxREAL);
    *mxGetPr(plhs[1]) = (double) i_present ;
  }
/* unsigned int ATSIF_GetStructureVersion(ATSIF_StructureElement _element, AT_U32 * _ui_versionHigh, AT_U32 * _ui_versionLow); =============================================================*/
  else if(!strcmp(func,"ATSIF_GetStructureVersion")) {
    AT_32 i_element;
    AT_U32 ui_versionH, ui_versionL;
    AndorParCheck(func, 3, 2, NULL, nlhs, plhs, nrhs, prhs);	 
    i_element = (int)*mxGetPr(prhs[1]);
    i_element += 0x40000000;  
    ret = ATSIF_GetStructureVersion((ATSIF_StructureElement)(i_element),&ui_versionH, &ui_versionL);
    plhs[1] = mxCreateDoubleMatrix(1,1, mxREAL);
    plhs[2] = mxCreateDoubleMatrix(1,1, mxREAL);
    *mxGetPr(plhs[1]) = (double) ui_versionH ;
    *mxGetPr(plhs[2]) = (double) ui_versionL ;		
  }
/* unsigned int ATSIF_GetFrameSize(ATSIF_DataSource _source, AT_U32 * _ui_size); ======================*/
  else if(!strcmp(func,"ATSIF_GetFrameSize")) {
    AT_U32 ui_size;
    AT_32 i_source;
    AndorParCheck(func, 2, 2, NULL, nlhs, plhs, nrhs, prhs);	 
    i_source = (int)*mxGetPr(prhs[1]);
    i_source += 0x40000000;  
    ret = ATSIF_GetFrameSize((ATSIF_DataSource)(i_source),&ui_size);
    plhs[1] = mxCreateDoubleMatrix(1,1, mxREAL);
    *mxGetPr(plhs[1]) = (double) ui_size ;  
  }
/* unsigned int ATSIF_GetNumberFrames(ATSIF_DataSource _source, AT_U32 * _ui_images); ======================*/
  else if(!strcmp(func,"ATSIF_GetNumberFrames")) {
    AT_U32 ui_noImages;
    AT_32 i_source;
    AndorParCheck(func, 2, 2, NULL, nlhs, plhs, nrhs, prhs);	 
    i_source = (int)*mxGetPr(prhs[1]);
    i_source += 0x40000000;  
    ret = ATSIF_GetNumberFrames((ATSIF_DataSource)(i_source),&ui_noImages);
    plhs[1] = mxCreateDoubleMatrix(1,1, mxREAL);
    *mxGetPr(plhs[1]) = (double) ui_noImages ;  
  }
/* unsigned int ATSIF_GetNumberSubImages(ATSIF_DataSource _source, AT_U32 * _ui_noSubImages); ======================*/
  else if(!strcmp(func,"ATSIF_GetNumberSubImages")) {
    AT_U32 ui_noSubImages;
    AT_32 i_source;
    AndorParCheck(func, 2, 2, NULL, nlhs, plhs, nrhs, prhs);	 
    i_source = (int)*mxGetPr(prhs[1]);
    i_source += 0x40000000;  
    ret = ATSIF_GetNumberSubImages((ATSIF_DataSource)(i_source),&ui_noSubImages);
    plhs[1] = mxCreateDoubleMatrix(1,1, mxREAL);
    *mxGetPr(plhs[1]) = (double) ui_noSubImages ;  
  }
/* unsigned int ATSIF_GetSubImageInfo(ATSIF_DataSource _source, AT_U32 _ui_index,
                                                   AT_U32 * _ui_left, AT_U32 * _ui_bottom,
                                                   AT_U32 * _ui_right, AT_U32 * _ui_top,
                                                   AT_U32 * _ui_hBin, AT_U32 * _ui_vBin); ======================*/
  else if(!strcmp(func,"ATSIF_GetSubImageInfo")) {
    AT_32 i_source;
    AT_U32 ui_index;
    AT_U32 ui_left, ui_right, ui_bottom, ui_top, ui_hBin, ui_vBin;
    AndorParCheck(func, 7, 3, NULL, nlhs, plhs, nrhs, prhs);	 
    i_source = (int)*mxGetPr(prhs[1]);
    i_source += 0x40000000;  
    ui_index = (int)*mxGetPr(prhs[2]);
    ret = ATSIF_GetSubImageInfo((ATSIF_DataSource)(i_source),ui_index, &ui_left, &ui_bottom, &ui_right, &ui_top, &ui_hBin, &ui_vBin);
    plhs[1] = mxCreateDoubleMatrix(1,1, mxREAL);
    plhs[2] = mxCreateDoubleMatrix(1,1, mxREAL);
    plhs[3] = mxCreateDoubleMatrix(1,1, mxREAL);
    plhs[4] = mxCreateDoubleMatrix(1,1, mxREAL);
    plhs[5] = mxCreateDoubleMatrix(1,1, mxREAL);
    plhs[6] = mxCreateDoubleMatrix(1,1, mxREAL);
    *mxGetPr(plhs[1]) = (double) ui_left;  
    *mxGetPr(plhs[2]) = (double) ui_bottom;  
    *mxGetPr(plhs[3]) = (double) ui_right;  
    *mxGetPr(plhs[4]) = (double) ui_top;   
    *mxGetPr(plhs[5]) = (double) ui_hBin;  
    *mxGetPr(plhs[6]) = (double) ui_vBin;                 		
  }
/* unsigned int ATSIF_GetAllFrames(ATSIF_DataSource _source, float * _pf_data, AT_U32 _ui_bufferSize); ========================================*/
  else if(!strcmp(func,"ATSIF_GetAllFrames")) {
    mwSize p[2];
    AT_U32 ui_size;
    AT_32 i_source;
    AndorParCheck(func, 2, 3, NULL, nlhs, plhs, nrhs, prhs);
    i_source = (int)*mxGetPr(prhs[1]);
    i_source += 0x40000000;  
    ui_size = *mxGetPr(prhs[2]);
    p[0] = (mwSize)ui_size;
    p[1] = 0;		/* BBB close the dimension count */
    plhs[1] = mxCreateNumericArray(1, p, mxSINGLE_CLASS, mxREAL);
    ret = ATSIF_GetAllFrames((ATSIF_DataSource)(i_source),mxGetData(plhs[1]), ui_size);
  }
/* unsigned int ATSIF_GetFrame(ATSIF_DataSource _source, AT_U32 _ui_index, float * _pf_data, AT_U32 _ui_bufferSize);========================================*/
  else if(!strcmp(func,"ATSIF_GetFrame")) {
    mwSize p[2];
    AT_U32 ui_size;
    AT_U32 ui_index;
    AT_32 i_source;
    AndorParCheck(func, 2, 4, NULL, nlhs, plhs, nrhs, prhs);
    i_source = (int)*mxGetPr(prhs[1]);
    i_source += 0x40000000;  
    ui_index = *mxGetPr(prhs[2]);
    ui_size = *mxGetPr(prhs[3]);
    p[0] = (mwSize)ui_size;
    p[1] = 0;		/* BBB close the dimension count */
    plhs[1] = mxCreateNumericArray(1, p, mxSINGLE_CLASS, mxREAL);
    ret = ATSIF_GetFrame((ATSIF_DataSource)(i_source),ui_index, mxGetData(plhs[1]), ui_size);
  }
/* unsigned int ATSIF_GetPropertyValue(ATSIF_DataSource _source,const AT_C * _sz_propertyName,AT_C * _sz_propertyValue,AT_U32 _ui_bufferSize); ========================================================*/
  else if(!strcmp(func,"ATSIF_GetPropertyValue")) {
    AT_32 i_source;
    char sz_propertyVal[STRLEN];
    char sz_propertyName[STRLEN];
    AndorParCheck(func, 2, 3, "SDS", nlhs, plhs, nrhs, prhs);
    i_source = (int)*mxGetPr(prhs[1]);
    i_source += 0x40000000;  
    mxGetString(prhs[2], sz_propertyName, STRLEN-1);
    ret = ATSIF_GetPropertyValue((ATSIF_DataSource)(i_source), sz_propertyName, sz_propertyVal, STRLEN-1);
	plhs[1] = mxCreateString(sz_propertyVal);			/* create string */
  }
/* unsigned int ATSIF_GetPropertyType(ATSIF_DataSource _source, const AT_C * _sz_propertyName, ATSIF_PropertyType * _propertyType); ========================================================*/
  else if(!strcmp(func,"ATSIF_GetPropertyType")) {
    AT_32 i_source;
    char sz_propertyName[STRLEN];
    AT_32 i_type;
    AndorParCheck(func, 2, 3, NULL, nlhs, plhs, nrhs, prhs);
    i_source = (int)*mxGetPr(prhs[1]);
    i_source += 0x40000000;  
    mxGetString(prhs[2], sz_propertyName, STRLEN-1);
    ret = ATSIF_GetPropertyType((ATSIF_DataSource)(i_source), sz_propertyName, (ATSIF_PropertyType*)&i_type);
    plhs[1] = mxCreateDoubleMatrix(1,1, mxREAL);
	i_type -= 0x40000000;
    *mxGetPr(plhs[1]) = (double) i_type;  
  }
/* unsigned int ATSIF_GetDataStartBytePosition(ATSIF_DataSource _source, AT_32 * _i_startPosition); ======================*/
  else if(!strcmp(func,"ATSIF_GetDataStartBytePosition")) {
    AT_32 i_startPos, i_source;
    AndorParCheck(func, 2, 2, NULL, nlhs, plhs, nrhs, prhs);	 
    i_source = (int)*mxGetPr(prhs[1]);
    i_source += 0x40000000;  
    ret = ATSIF_GetDataStartBytePosition((ATSIF_DataSource)(i_source),&i_startPos);
    plhs[1] = mxCreateDoubleMatrix(1,1, mxREAL);
    *mxGetPr(plhs[1]) = (double) i_startPos ;  
  }
/* unsigned int ATSIF_GetPixelCalibration(ATSIF_DataSource _source, ATSIF_CalibrationAxis _axis, AT_32 _i_pixel, double * _d_calibValue); ======================*/
  else if(!strcmp(func,"ATSIF_GetPixelCalibration")) {
    AT_32 i_startPos, i_source, i_calibAxis, i_pix;
    double d_calibValue;
    AndorParCheck(func, 2, 4, NULL, nlhs, plhs, nrhs, prhs);	 
    i_source = (int)*mxGetPr(prhs[1]);
    i_source += 0x40000000; 
    i_calibAxis = (int)*mxGetPr(prhs[2]);
    i_calibAxis += 0x40000000;  
    i_pix = (int)*mxGetPr(prhs[3]);
    ret = ATSIF_GetPixelCalibration((ATSIF_DataSource)(i_source), (ATSIF_CalibrationAxis)(i_calibAxis), i_pix, &d_calibValue);
    plhs[1] = mxCreateDoubleMatrix(1,1, mxREAL);
    *mxGetPr(plhs[1]) = (double) d_calibValue ;  
  }
  else {
    mexErrMsgTxt("Command Unknown");
  }
  /* Return, ALWAYS, the status of call */
  plhs[0] = mxCreateDoubleMatrix(1,1, mxREAL);
  *mxGetPr(plhs[0]) = (double) ret ;  
  return;
}
