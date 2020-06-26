//
// Created by antoine on 26/06/2020.
//

#ifndef GVDB_ALL_HTGTOGVDBCONVERTER_H
#define GVDB_ALL_HTGTOGVDBCONVERTER_H

class openHyperTreeGrid;
namespace nvdb
{
	class VolumeGVDB;
}
class OpenHTGtoGVDBConverter
{
public:
	static void compute(openHyperTreeGrid & htgIn, nvdb::VolumeGVDB & gvdbIn, const char* scalarNameIn);
};


#endif //GVDB_ALL_HTGTOGVDBCONVERTER_H
