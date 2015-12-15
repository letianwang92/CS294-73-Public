#include "WriteRectMDArray.H"
#include "RectMDArray.H"
#include "VisitWriter.H"
#include "Real.H"
#include <cstdio>

const char* MDWrite(RectMDArray<float>* a_array)
{
  static int fileCount = 0;
  static char nameBuffer[10];
  if(a_array == NULL)
    {
      return nameBuffer;
    }
  sprintf(nameBuffer, "md%d",fileCount);
  MDWrite(nameBuffer, a_array);
  fileCount++;
  return nameBuffer;
}

const char* MDWrite(RectMDArray<double>& a_array)
{
  Box bx = a_array.getBox();
  RectMDArray<float> array(bx);
  for (int k = 0; k < bx.sizeOf();k++)
    {
      array[k] = (float) a_array[k];
    }
  return MDWrite(&array);
  
}

void MDWrite(const char* a_filename, RectMDArray<float>* a_array)
{
  if(a_filename == NULL || a_array == NULL)
    {
      return;
    }
  int dim[3] = {1,1,1};
  int vardims[1] ={1};
  int centering[1]={1};
  float* vars[1];
 
  const char * const varnames[] = { "cellCentered" };
  Point lo, hi;
  const Box& box = a_array->getBox();
  lo = box.getLowCorner();
  hi = box.getHighCorner();
  for(int i=0; i<DIM;i++)
    {
      dim[i] = hi[i]-lo[i]+1;
    }
  float& val = a_array->getPointer()[0];
  vars[0] = &val;
  write_regular_mesh(a_filename, 1, dim, 1, vardims, centering,  varnames, vars);
}

void MDWrite(const char* filename, RectMDArray<double>& a_array)
{
  Box bx = a_array.getBox();
  RectMDArray<float> array(bx);
  for (int k = 0; k < bx.sizeOf();k++)
    {
      array[k] = (float) a_array[k];
    }
  MDWrite(filename,&array);
}



