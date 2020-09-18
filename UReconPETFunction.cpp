/********************************************************************

Copyright (C), 2019, All rights reserved

File Name     :    UReconPETFunction.cpp
Description   :
History       :

<author>            <time>            <desc>
Ang Li             2019/6/8           create

********************************************************************/
#include "UReconPETFunction.h"
#include <math.h>



UReconPETFunction::UReconPETFunction()
{
}


UReconPETFunction::~UReconPETFunction()
{
}
/**********************************************************
Description:		Well counter correction of reconstructed image
Arguments:
{
	pImage:					Image to be corrected.
	nImgVoxelNum:			Number of image voxel.
	WellCounterSlope:		Slope of the well counter linear transformation
	WellCounterIntercept:	Intercept of the well counter linear transformation
}
Return:				Status.
**********************************************************/
bool UReconPETFunction::WellCounterCorrection(float * pImage, int nImgVoxelNum, float WellCounterSlope, float WellCounterIntercept)
{
	if (WellCounterSlope <= 0)
		return false;
	for (int i = 0; i < nImgVoxelNum; i++)
	{
		pImage[i] = pImage[i] * WellCounterSlope + WellCounterIntercept;
	}
	return true;
}
/**********************************************************
Description:		Decay correction of reconstructed image. Calibrate activity to the time of injection.
Arguments:
{
	pImage:			Image to be corrected.
	nImgVoxelNum:	Number of image voxel.
	timeBeforeScan:	The time from the injection to the start of the scan
	timeDuringScan:	Scanning time
	halfLife:		The half-life of a nuclide
	brachingRatio:	The braching ratio of a nuclide
}
Return:				Status.
**********************************************************/
bool UReconPETFunction::DecayCorrection(float * pImage, int nImgVoxelNum, float timeBeforeScan, float timeDuringScan, float halfLife, float brachingRatio)
{
	if (timeBeforeScan < 0 || timeDuringScan < 0)
	{
		return false;
	}
	float coff;
	if (timeDuringScan < 0.001f)
		coff = powf(2.0f, timeBeforeScan / halfLife) / brachingRatio;//20190617Updated.之前的写错了 //201908025Updated.加入分支比
	else
		coff = logf(2.0f) / halfLife * timeDuringScan / powf(0.5f, timeBeforeScan / halfLife) / (1 - powf(0.5f, timeDuringScan / halfLife)) / brachingRatio;
	for (int i = 0; i < nImgVoxelNum; i++)
	{
		pImage[i] *= coff;
	}
	return true;
}
/**********************************************************
Description:		Calculate the average activity by dividing the scanning time.
Arguments:
{
	pImage:			Image to be corrected.
	nImgVoxelNum:	Number of image voxel.
	timeDuringScan:	Scanning time
}
Return:				Status.
**********************************************************/
bool UReconPETFunction::CalcAverageActivity(float *pImage, int nImgVoxelNum, float timeDuringScan)
{
	if (timeDuringScan <= 0)
	{
		return false;
	}
	for (int i = 0; i < nImgVoxelNum; i++)
	{
		pImage[i] /= timeDuringScan;
	}
	return true;
}


