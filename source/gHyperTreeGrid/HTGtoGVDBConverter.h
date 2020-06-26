//
// Created by antoine on 26/06/2020.
//

#ifndef GVDB_ALL_HTGTOGVDBCONVERTER_H
#define GVDB_ALL_HTGTOGVDBCONVERTER_H

class vtkHyperTreeGrid;
namespace nvdb
{
	class VolumeGVDB;
}
class HTGtoGVDBConverter
{
public:
	static void compute(vtkHyperTreeGrid * htgIn, nvdb::VolumeGVDB & gvdbIn, const char* scalarNameIn);
};


#endif //GVDB_ALL_HTGTOGVDBCONVERTER_H
