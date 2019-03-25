// median_filter.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <stdlib.h>

#define MATLAB

#ifdef MATLAB
#include "mex.h"
#include <math.h>
#endif


typedef struct RECORD {
	double val;
	int next;
	int prev;
};

typedef struct _value {
	float val;
	int idx;	
} VALUE;

typedef double elem_type;
int compare ( const void * a, const void * b )
{
	return (( (*( VALUE*)b).val - (*( VALUE*)a).val )>=0) ? -1 : 1 ;
}


void my_median(double in_array[],double out_array[],int array_len,int window_len,int kth_value,int offset)
{
	// define nodes array
	RECORD *nodes = new RECORD[window_len];
	VALUE *local_array  =new VALUE [window_len];
	int last_node_idx = 0;
	int old_median_node;
	int head_node,tail_node;

	int node_idx;
	for (int i=0;i<window_len;i++){
		nodes[i].val = in_array[i];
		local_array[i].val=in_array[i];
		local_array[i].idx = i;
	}

	// sort the first window_len samples
	qsort(local_array, window_len, sizeof(elem_type), compare);

	// initiate node values
	for (int i=0;i<window_len;i++){
		nodes[i].val = in_array[i];
	}
	// forward connection between noodes
	for (int i=0;i<window_len;i++){
		if (i<window_len-1)
			nodes[local_array[i].idx].next = local_array[i+1].idx;
		else{
			nodes[local_array[i].idx].next = -1;
			tail_node = local_array[i].idx;
		}
	}
	// backword connection between noodes
	for (int i=0;i<window_len;i++){
		if (i==kth_value)
			old_median_node = local_array[i].idx;
		if (i>0)
			nodes[local_array[i].idx].prev = local_array[i-1].idx;
		else{
			nodes[local_array[i].idx].prev = -1;
			head_node = local_array[i].idx;
		}
	}



	
#ifndef MATLAB
	int node_ptr ;
	node_ptr = head_node;
	printf("%2.0d %3.0f ",window_len-1,in_array[window_len-1]);
	while (nodes[node_ptr].next!=-1){
		printf("%3.0f ",nodes[node_ptr].val);
		node_ptr=nodes[node_ptr].next;
	} 
	printf("%3.0f MV:%3.0f\n",nodes[node_ptr].val,nodes[old_median_node].val);
#endif
	out_array[window_len+offset-1] = nodes [old_median_node].val;
	for (int i=window_len;i<array_len;i++)	{

		if (nodes[old_median_node].val<=in_array[i]) {
			// if old median is smaller then new value. 
			// this code looks for the first node that is locatef after the median and
			// its value is samller the the new value.
			node_idx = old_median_node;
			while (nodes[node_idx].next!=-1 && nodes[nodes[node_idx].next].val<=in_array[i] )
				node_idx=nodes[node_idx].next;
		}else

			if (nodes[old_median_node].val>in_array[i]) {
				// if old median is smaller then new value. 
				// this code looks for the first node that is locatef after the median and
				// its value is samller the the new value.
				node_idx = old_median_node;
				while (nodes[node_idx].prev!=-1 && nodes[nodes[node_idx].prev].val>=in_array[i] )
					node_idx=nodes[node_idx].prev;
			}
			// if the new number is placed on the sampe place as the current removed number
			// don't change the pointers
			if (node_idx!=last_node_idx){
				// disoconnect old from itt neighborus
				if (nodes[last_node_idx].next!=-1)
					nodes [ nodes[last_node_idx].next ].prev = nodes[last_node_idx].prev;
				else{
					//	nodes [ nodes[last_node_idx].next ].prev  = -1;
					tail_node = nodes[last_node_idx].prev;
				}
				if (nodes[last_node_idx].prev!=-1)
					nodes [ nodes[last_node_idx].prev ].next = nodes[last_node_idx].next; 
				else{
					//	nodes [ nodes[last_node_idx].prev ].next  = -1;
					head_node = nodes[last_node_idx].next;
				}


				// reconnect new node
				if (nodes[old_median_node].val<=in_array[i]) {
					// if old median is smaller then new value. 
					// this code looks for the first node that is locatef after the median and
					// its value is samller the the new value.

					nodes[last_node_idx].prev = node_idx;
					nodes[last_node_idx].next = nodes[node_idx].next;

					if (nodes[node_idx].next!=-1)
						nodes[nodes[node_idx].next].prev = last_node_idx;
					nodes[node_idx].next = last_node_idx;

					
					if (nodes[last_node_idx].next==-1)
						// if new value is add at the end of the vector - set it as tail
						tail_node = last_node_idx;


				}else
					if (nodes[old_median_node].val>in_array[i]) {
						// if old median is smaller then new value. 
						// this code looks for the first node that is locatef after the median and
						// its value is samller the the new value.

						nodes[last_node_idx].next = node_idx;
						nodes[last_node_idx].prev = nodes[node_idx].prev;

						if (nodes[node_idx].prev!=-1) 
							nodes[nodes[node_idx].prev].next = last_node_idx;
						nodes[node_idx].prev = last_node_idx;

						// if new value is add at the end of the vector - set it as tail
						if (nodes[last_node_idx].prev==-1)
							head_node = last_node_idx;

					}

			}





			// set new value to the connected node
			nodes[last_node_idx].val = in_array[i];

			// update median node
			old_median_node = head_node;
			for (int ii=0;ii<kth_value;ii++)
				old_median_node=nodes[old_median_node].next;

			// update output vector 
			out_array[offset+i] = nodes [old_median_node].val;


#ifndef MATLAB
			int node_ptr;
			node_ptr = head_node;
			printf("%2.0d %3.0f ",i,in_array[i]);
			while (nodes[node_ptr].next!=-1){
				printf("%3.0f ",nodes[node_ptr].val);
				node_ptr=nodes[node_ptr].next;
			} 
			printf("%3.0f MV:%3.0f\n",nodes[node_ptr].val,nodes[old_median_node].val);
#endif
			last_node_idx++;
			last_node_idx%=window_len;
	}
	// free mem 
	delete local_array;
	delete nodes;
}


// #ifndef MATLB
// int _tmain(int argc, _TCHAR* argv[])
// {
// 	double in_array[33]={15,2,25,29,8,28,42,21,31,26,17,38,14,49,19,16,35,39,36,46,34,43,7,22,18,3,28,9,47,12};
// 	double median_array[33];
// 	int window_len=5;
// 	int kth_value=2;
// 	int offset = -2;
// 	int array_len = 33;
// 	my_median(in_array, median_array, array_len,window_len, kth_value,offset);
// 	return 0;
// }
// #endif

#ifdef MATLAB
void mexFunction( int nlhs, mxArray *plhs[],
				 int nrhs, const mxArray *prhs[] )
{
	int i,j,k,array_len,window_len,kth_value,offset;
	double *median_array,*x_diff,*in_array;

	/* get input vector contains values of K-S statistics*/
	if (nrhs!=2) {
		mexPrintf("\nThe number of variables entered is: %d\n",nrhs);
		mexErrMsgTxt("The correct format is:\n[median_array]=qumran_median(in_array,window_len)"); 
	}
	/* length of first parametr */
	array_len = mxGetDimensions(prhs[0])[1];  
	if (array_len==1)
		array_len = mxGetDimensions(prhs[0])[0];

	in_array = mxGetPr( prhs[0] );
	window_len = (int)*mxGetPr( prhs[1] );

	kth_value = (int)window_len/2;
	offset = -kth_value;
	// output
	plhs[0] = mxCreateDoubleMatrix(1,array_len, mxREAL);
	median_array = mxGetPr( plhs[0] );

	my_median(in_array, median_array, array_len,window_len, kth_value,offset);            

	//   mxFree(x_diff);

}

#endif

