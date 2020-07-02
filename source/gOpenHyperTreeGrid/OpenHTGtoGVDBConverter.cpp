//
// Created by antoine on 18/06/2020.
//

#include "OpenHTGtoGVDBConverter.h"
#include "gvdb.h"

//#undef AllValues

#include "openhtg/openHyperTreeGridCursor.h"
#include "openhtg/openHyperTreeGrid.h"
#include "openhtg/openHyperTree.h"

const unsigned char HTG_ACTIVATE_SPACE = 0;
const unsigned char HTG_FILL_INFO = 1;

typedef struct
{
	float sceneFactor;
	std::vector<nvdb::Vector3DF> htgPositions;
	std::vector<uint32_t> htgColors;
	nvdb::VolumeGVDB *gvdb;
	openHyperTreeGrid *openHTG;
	Vector3DF sceneOffset;
	uint64_t numberOfVoxel{0};
} data_t;


void activateSpace(const Vector3DF &pos, unsigned, data_t &data)
{
	data.gvdb->ActivateSpace(pos);
}

void fillInfo(const Vector3DF &pos, uint32_t color, data_t &data)
{
	data.htgPositions.push_back(pos);
	data.htgColors.push_back(color);
}

void processVoxel(const Vector3DF &minCursor, const Vector3DF &maxCursor, uint32_t color, data_t & data,
                  std::function<void(const Vector3DF &pos, unsigned color, data_t &data)> func)
{
	for (float i = minCursor.x; i < maxCursor.x; ++i)
	{
		for (float j = minCursor.y; j < maxCursor.y; ++j)
		{
			for (float k = minCursor.z; k < maxCursor.z; ++k)
			{
				Vector3DF pos(i, j, k);
				func(pos, color, data);
				++data.numberOfVoxel;
			}
		}

	}
}


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
	data.gvdb->AllocData(htgPosPtr, numberOfVoxels, sizeof(nvdb::Vector3DF));
	data.gvdb->AllocData(htgColorPtr, numberOfVoxels, sizeof(uint));
	{
		nvdb::Vector3DF *tmpPtr = (nvdb::Vector3DF *) data.gvdb->getDataPtr(0, htgPosPtr);
		std::copy(data.htgPositions.begin(), data.htgPositions.end(), tmpPtr);
		data.htgPositions.clear();
	}

	{
		uint *tmpPtr = (uint *) data.gvdb->getDataPtr(0, htgColorPtr);
		std::copy(data.htgColors.begin(), data.htgColors.end(), tmpPtr);

		data.htgColors.clear();
	}

	data.gvdb->CommitData(htgColorPtr);
	data.gvdb->CommitData(htgPosPtr);

	data.gvdb->SetDataGPU(pntpos, numberOfVoxels, htgPosPtr.gpu, 0, sizeof(nvdb::Vector3DF));
	data.gvdb->SetDataGPU(pntclr, numberOfVoxels, htgColorPtr.gpu, 0, sizeof(uint));

	DataPtr dataP;
	data.gvdb->SetPoints(pntpos, dataP, pntclr);

	int scPntLen = 0;
	int subcell_size = 8;
	float radius = 1.f;

	data.gvdb->InsertPointsSubcell(subcell_size, numberOfVoxels, radius, nvdb::Vector3DF(0.f, 0.f, 0.f), scPntLen);
	data.gvdb->GatherDensity(subcell_size, numberOfVoxels, radius, nvdb::Vector3DF(0, 0, 0), scPntLen, 0, 1,
	                         true); // true = accumulate
}

void dfsHT(openHyperTreeGridCursorNonOrientedGeometry &cursor,
           unsigned int operation, data_t &data)
{

	if (cursor.isMasked())
	{
		return;
	}
	if (cursor.isLeaf())
	{
		double bounds[6];
		cursor.getBounds(bounds);
		Vector3DF minCursor(bounds[0], bounds[2], bounds[4]);
		Vector3DF maxCursor(bounds[1], bounds[3], bounds[5]);
		minCursor *= data.sceneFactor;
		minCursor -= data.sceneOffset;

		maxCursor *= data.sceneFactor;
		maxCursor -= data.sceneOffset;

		if (operation == HTG_FILL_INFO)
		{
			uint32_t color = 0x00FFFF00;
			processVoxel(minCursor, maxCursor, color, data, fillInfo);
		}
		else if(operation == HTG_ACTIVATE_SPACE)
		{
			processVoxel(minCursor, maxCursor, 0, data, activateSpace);
		}


	}
	else
	{
		const unsigned int nbChild = cursor.getNumberOfChildren();
		for (unsigned int i = 0; i < nbChild; i++)
		{
			cursor.toChild(i);
			dfsHT(cursor, operation, data);
			cursor.toParent();
		}
	}
};


void OpenHTGtoGVDBConverter::compute(openHyperTreeGrid &openHTGIn, nvdb::VolumeGVDB &gvdbIn, const char *scalarName)
{
	data_t data;
	data.openHTG = &openHTGIn;
	data.gvdb = &gvdbIn;
	nvdb::Vector3DF minBounds;
	nvdb::Vector3DF maxBounds;
	unsigned int count = 0;

	auto cursor = openHyperTreeGridCursorNonOrientedGeometry(*data.openHTG);
	for (uint32_t i = 0; i < data.openHTG->getMaxNumberOfTrees(); ++i)
	{
		if (!cursor.toTree(i))
		{
			continue;
		}
		double bounds[6];
		cursor.getBounds(bounds);
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
	data.openHTG->getBounds(htgBounds);

	unsigned int htgDimensions[3];
	data.openHTG->getDimensions(htgDimensions);

	nvdb::Vector3DF nbTrees(htgDimensions[0], htgDimensions[1], htgDimensions[2]);

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

	unsigned int htgNbLevel = data.openHTG->getNumberOfLevels();

	unsigned int sceneSize = tmp * (int) pow(2, htgNbLevel);

	unsigned int nbLevels = log2(sceneSize);

	data.sceneFactor = sceneSize / (biggestAxe * (tmp / std::max(nbTrees.x, std::max(nbTrees.y,
	                                                                                 nbTrees.z))));//std::max(htgBounds[1] - htgBounds[0],
	//        std::max(htgBounds[3] - htgBounds[2], htgBounds[5] - htgBounds[4]));
	data.sceneOffset = minBounds * data.sceneFactor;
	bool bnew;

	unsigned int gridArgs[5] = {0, 0, 0, 0, 3};
	for (unsigned int i = 0; i < std::min(nbLevels - 3, 0u); i++)
	{
		++gridArgs[i % 4];
	}

	data.gvdb->Configure(gridArgs[0], gridArgs[1], gridArgs[2], gridArgs[3], gridArgs[4]);

	for (uint32_t i = 0; i < data.openHTG->getMaxNumberOfTrees(); ++i)
	{
		if (cursor.toTree(i))
		{
			dfsHT(cursor, HTG_ACTIVATE_SPACE, data);
		}
	}
	// Finish Topology
	data.gvdb->FinishTopology();

	data.gvdb->DestroyChannels();
	data.gvdb->AddChannel(0, T_FLOAT, 0, F_LINEAR);
	data.gvdb->AddChannel(1, T_UCHAR4, 0, F_LINEAR);
	data.gvdb->UpdateAtlas();

	unsigned int numberOfVoxels;

	for (uint32_t i = 0; i < data.openHTG->getMaxNumberOfTrees(); ++i)
	{
		if (!cursor.toTree(i))
		{
			continue;
		}
		dfsHT(cursor, HTG_FILL_INFO, data);

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

	std::cout << "Number of voxels in GVDB structure : " << data.numberOfVoxel << std::endl;

	data.gvdb->UpdateApron();

	data.gvdb->SetColorChannel(1);
}


//uint32_t calculColor(unsigned int index)
//{
//	if (!scalars)
//	{
//		return 0x00808080;
//	}
//	uint32_t tmp = 0;
//
//	const unsigned char *rgb = lut->MapValue(*(scalars->GetTuple(index)));
//	tmp |= rgb[0];
//	tmp |= rgb[1] << 8U;
//	tmp |= rgb[2] << 16;
//	return tmp;
//}




