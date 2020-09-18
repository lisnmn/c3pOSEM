/********************************************************************

Copyright (C), 2019, All rights reserved

File Name     :    UReconPETFunction.h
Description   :
History       :

<author>            <time>            <desc>
Ang Li             2019/6/8           create

********************************************************************/
#ifndef __RECON_PET_FUNCTION_H
#define __RECON_PET_FUNCTION_H


class UReconPETFunction
{
public:
	UReconPETFunction();
	virtual ~UReconPETFunction();
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
	static bool WellCounterCorrection(float *pImage, int nImgVoxelNum, float WellCounterSlope, float WellCounterIntercept);
	/**********************************************************
	Description:		Decay correction of reconstructed image. Calibrate average activity to the time of injection.
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
	static bool DecayCorrection(float *pImage, int nImgVoxelNum, float timeBeforeScan, float timeDuringScan, float halfLife, float brachingRatio);
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
	static bool CalcAverageActivity(float *pImage, int nImgVoxelNum, float timeDuringScan);
private:
	UReconPETFunction(const UReconPETFunction&){};
	UReconPETFunction& operator=(const UReconPETFunction&){};
};

#endif // !__PET_RECON_FUNCTION_H
