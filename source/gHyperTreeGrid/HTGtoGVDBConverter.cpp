//
// Created by antoine on 26/06/2020.
//

#include "vtkHyperTreeGrid.h"

#include "HTGtoGVDBConverter.h"
#include "gvdb.h"

#include <vtkHyperTreeGridNonOrientedGeometryCursor.h>
#include <vtkHyperTree.h>
#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkLookupTable.h>

const uint8_t HTG_ACTIVATE_SPACE = 0;
const uint8_t HTG_FILL_INFO = 1;

typedef struct
{
	float sceneSize;
	float sceneFactor;
	std::vector<nvdb::Vector3DF> htgPositions;
	std::vector<uint32_t> htgColors;
	nvdb::VolumeGVDB *gvdb;
	vtkHyperTreeGrid *htg;
	Vector3DF sceneOffset;
	vtkDataArray *scalars;
	vtkLookupTable *lut;
} data_t;


uint32_t calculColor(data_t &data, uint32_t index)
{
	if (!data.scalars)
	{
		return 0x00808080;
	}
	uint32_t tmp = 0;

	const unsigned char *rgb = data.lut->MapValue(*(data.scalars->GetTuple(index)));
	tmp |= rgb[0];
	tmp |= rgb[1] << 8U;
	tmp |= rgb[2] << 16;
	return tmp;
}


void dfsHT(data_t &data,
           vtkHyperTreeGridNonOrientedGeometryCursor *cursor,
           unsigned int operation)
{

	if (cursor->IsMasked())
	{
		return;
	}
	if (cursor->IsLeaf())
	{
		double bounds[6];
		cursor->GetBounds(bounds);
		Vector3DF minCursor(bounds[0], bounds[2], bounds[4]);
		Vector3DF maxCursor(bounds[1], bounds[3], bounds[5]);
		minCursor *= data.sceneFactor;
		minCursor -= data.sceneOffset;

		maxCursor *= data.sceneFactor;
		maxCursor -= data.sceneOffset;

		uint color;
		if (operation == HTG_FILL_INFO)
		{
			color = calculColor(data, cursor->GetGlobalNodeIndex());
		}

		for (float i = minCursor.x; i < maxCursor.x; ++i)
		{
			for (float j = minCursor.y; j < maxCursor.y; ++j)
			{
				for (float k = minCursor.z; k < maxCursor.z; ++k)
				{
					Vector3DF pos(i, j, k);
					if (operation == HTG_ACTIVATE_SPACE)
					{
						data.gvdb->ActivateSpace(pos);
					}
					else if (operation == HTG_FILL_INFO)
					{
						data.htgPositions.push_back(pos);
						data.htgColors.push_back(color);
					}
				}
			}

		}


	}
	else
	{
		unsigned int nbChild = cursor->GetNumberOfChildren();
		for (unsigned int i = 0; i < nbChild; i++)
		{
			cursor->ToChild(i);
			dfsHT(data, cursor, operation);
			cursor->ToParent();
		}
	}
};

void pushNodesIntoGPU(data_t &data)
{
	DataPtr htgPosPtr;
	DataPtr htgColorPtr;

	unsigned int numberOfVoxels = data.htgPositions.size();

	if (!numberOfVoxels)
	{
		return;
	}

	DataPtr pntpos, pntclr;
	data.gvdb->AllocData(htgPosPtr, numberOfVoxels, sizeof(Vector3DF));
	data.gvdb->AllocData(htgColorPtr, numberOfVoxels, sizeof(uint));
	{
		Vector3DF *tmpPtr = (Vector3DF *) data.gvdb->getDataPtr(0, htgPosPtr);

		for (unsigned int i = 0; i < numberOfVoxels; ++i)
		{
			tmpPtr[i] = data.htgPositions[i];
		}
		data.htgPositions.clear();
	}
	{
		uint *tmpPtr = (uint *) data.gvdb->getDataPtr(0, htgColorPtr);
		for (unsigned int i = 0; i < numberOfVoxels; ++i)
		{
			tmpPtr[i] = data.htgColors[i];
		}
		data.htgColors.clear();
	}

	data.gvdb->CommitData(htgColorPtr);
	data.gvdb->CommitData(htgPosPtr);

	data.gvdb->SetDataGPU(pntpos, numberOfVoxels, htgPosPtr.gpu, 0, sizeof(Vector3DF));
	data.gvdb->SetDataGPU(pntclr, numberOfVoxels, htgColorPtr.gpu, 0, sizeof(uint));

	DataPtr dataPtr;
	data.gvdb->SetPoints(pntpos, dataPtr, pntclr);

	int scPntLen = 0;
	int subcell_size = 8;
	float radius = 1.f;

	data.gvdb->InsertPointsSubcell(subcell_size, numberOfVoxels, radius, Vector3DF(0.f, 0.f, 0.f), scPntLen);
	data.gvdb->GatherDensity(subcell_size, numberOfVoxels, radius, Vector3DF(0, 0, 0), scPntLen, 0, 1,
	                         true); // true = accumulate
}


void HTGtoGVDBConverter::compute(vtkHyperTreeGrid *htgIn, nvdb::VolumeGVDB &gvdbIn, const char *scalarNameIn)
{

	data_t data;

	data.gvdb = &gvdbIn;
	data.htg = htgIn;

	data.scalars = data.htg->GetCellData()->GetScalars(scalarNameIn);

	if (data.scalars)
	{
		std::cout << "\"" << scalarNameIn << "\" scalar exists" << std::endl;
		double *dataRange = data.scalars->GetRange();

		data.lut = vtkLookupTable::New();
		data.lut->SetHueRange(0.66, 0);
		data.lut->SetTableRange(dataRange[0], dataRange[1]);
		data.lut->Build();
	}

	auto cursor = vtkHyperTreeGridNonOrientedGeometryCursor::New();
	auto treeIterator = vtkHyperTreeGrid::vtkHyperTreeGridIterator();
	vtkHyperTree *hyperTree = nullptr;

	treeIterator.Initialize(data.htg);

	Vector3DF minBounds;
	Vector3DF maxBounds;
	unsigned int count = 0;

	while ((hyperTree = treeIterator.GetNextTree()))
	{
		cursor->Initialize(data.htg, hyperTree->GetTreeIndex());
		double bounds[6];
		cursor->GetBounds(bounds);
		if (count++ == 0)
		{
			minBounds.Set(bounds[0], bounds[2], bounds[4]);
			maxBounds.Set(bounds[1], bounds[3], bounds[5]);
		}
		else
		{
			minBounds.Set(std::min(minBounds.x, (float) bounds[0]), std::min(minBounds.y, (float) bounds[2]),
			              std::min(minBounds.z, (float) bounds[4]));
			maxBounds.Set(std::max(maxBounds.x, (float) bounds[1]), std::max(maxBounds.y, (float) bounds[3]),
			              std::max(maxBounds.z, (float) bounds[5]));

		}
	}

	float biggestAxe = std::max(std::max(maxBounds.x - minBounds.x, maxBounds.y - minBounds.y),
	                            maxBounds.z - minBounds.z);

	double htgBounds[6];
	data.htg->GetBounds(htgBounds);

	unsigned int htgDimensions[3];
	data.htg->GetDimensions(htgDimensions);

	Vector3DF nbTrees(htgDimensions[0] - 1, htgDimensions[1] - 1, htgDimensions[2] - 1);

	double totalHTGBiggestAxe = std::max(htgBounds[1] - htgBounds[0],
	                                     std::max(htgBounds[3] - htgBounds[2], htgBounds[5] - htgBounds[4]));

	float diffFactor = biggestAxe / totalHTGBiggestAxe;

	nbTrees = nbTrees * (float) diffFactor;

	unsigned int tmp = 1;

	while (tmp < nbTrees.x || tmp < nbTrees.y || tmp < nbTrees.z)
	{
		tmp *= 2;
	}

	// For the voxel size we considere that the HTG is uniform, for simplicity

	unsigned int htgNbLevel = data.htg->GetNumberOfLevels();

	unsigned int sceneSize = tmp * (int) pow(2, htgNbLevel);

	unsigned int nbLevels = log2(sceneSize);

	data.sceneFactor = sceneSize / (biggestAxe * (tmp / std::max(nbTrees.x, std::max(nbTrees.y, nbTrees.z))));
	data.sceneOffset = minBounds * data.sceneFactor;
	bool bnew;

	unsigned int gridArgs[5] = {0, 0, 0, 0, 3};
	for (unsigned int i = 0; i < nbLevels - 3; i++)
	{
		++gridArgs[i % 4];
	}

	data.gvdb->Configure(gridArgs[0], gridArgs[1], gridArgs[2], gridArgs[3], gridArgs[4]);

	treeIterator.Initialize(data.htg);
	while ((hyperTree = treeIterator.GetNextTree()))
	{
		cursor->Initialize(data.htg, hyperTree->GetTreeIndex());
		dfsHT(data, cursor, HTG_ACTIVATE_SPACE);
	}

	// Finish Topology
	data.gvdb->FinishTopology();

	data.gvdb->DestroyChannels();
	data.gvdb->AddChannel(0, T_FLOAT, 0, F_LINEAR);
	data.gvdb->AddChannel(1, T_UCHAR4, 0, F_LINEAR);
	data.gvdb->UpdateAtlas();

	unsigned int numberOfVoxels;
	count = 0;
	treeIterator.Initialize(data.htg);
	while ((hyperTree = treeIterator.GetNextTree()))
	{
		cursor->Initialize(data.htg, hyperTree->GetTreeIndex());
		dfsHT(data, cursor, HTG_FILL_INFO);

		numberOfVoxels = data.htgPositions.size();
		if (numberOfVoxels < 1000000)
		{
			continue;
		}
		pushNodesIntoGPU(data);
		data.htgPositions.reserve(1000000);
		data.htgColors.reserve(1000000);
	}

	pushNodesIntoGPU(data);

	data.gvdb->UpdateApron();

	data.gvdb->SetColorChannel(1);
}

