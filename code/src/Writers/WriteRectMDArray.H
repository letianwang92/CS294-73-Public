#ifndef _WRITERECTMDARRAY_H_
#define _WRITERECTMDARRAY_H_
#include "RectMDArray.H"
#include "Real.H"

#include <cstring>
const char* MDWrite(RectMDArray<float>* a_array);
const char* MDWrite(RectMDArray<double>& a_array);
void MDWrite(const char* filename, RectMDArray<float>* a_array);
void MDWrite(const char* filename, RectMDArray<double>& a_array);
#endif 
