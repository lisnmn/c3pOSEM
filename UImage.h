/********************************************************************

Copyright (C), 2019, All rights reserved

File Name     :    UImage.h
Description   :
History       :

<author>            <time>            <desc>
Ang Li             2019/6/8           create

********************************************************************/
#ifndef __IMAGE_H_
#define __IMAGE_H_

#include "UVolume.h"


class UPETScanner;
class UImage
{


public:
	UImage(void);
	virtual ~UImage(void);
	UImage& operator=(const UImage& img);
	
public:
	//get
	int GetVoxelNumX()const;
	int GetVoxelNumY()const;
	int GetVoxelNumZ()const;
	int GetVoxelNumInAll()const;
	UVolume3D<int> GetVoxelNum()const;
	UVolume3D<double> GetVoxelSize()const;
	UVolume3D<double> GetVoxelCenter()const;
	//(set)configuration function
	void SetCenter(double x, double y, double z);
	bool SetVoxelSize(double x, double y, double z);
	bool SetVoxelNum(int x, int y, int z);
	void SetCylinder(UVolume3D<double> pos, double r, double h, float value = 1);
	//======================================
	// High level getting
	//======================================
	/**********************************************************
	Description:		It will be used when calculate the Solid Angle.
	Arguments:
	{
	index:			The index of the pixel = x + y*XAll + z*XAll*YAll.
	position1:		One edge point of this voxel.
	position2:		Another edge point of this voxel.
	}
	Return:				void.
	**********************************************************/
	void GetVoxelEdgePos(int index, UVolume3D<double> &position1, UVolume3D<double> &position2)const;
	/**********************************************************
	Description:		It will be used when calculate the Solid Angle.
	Arguments:
	{
	index:			The index of the pixel = x + y*XAll + z*XAll*YAll.
	subNumOneLine:	The subvoxel number in one direction, and the total subvoxel number is the cube of subNumOneLine.
	centre:			A array of Positions, and the array size is cube of subNumOneLine.
	}
	Return:				void.
	**********************************************************/
	void GetSubvoxelCentre(int index, int subNumOneLine, UVolume3D<double> *centre)const;//计算分割像素的中心坐标
	//image processing
	/**********************************************************
	Description:		3D Total Variation (TV) smoothing.
	Arguments:
	{
	iter:			The iteration number of TV smoothing.
	dt:				the step lengh of TV.
	epsilon:		A very small number to avoid dividing 0.
	lamdba:			A weight to keep the similarrity between the row image and the smoothed image.
	}
	Return:				Status.
	**********************************************************/
	bool TotalVariationDenoising3D(int iter, float dt, float epsilon = 1, float lamda = 0)const;//全变分滤波
	/**********************************************************
	Description:		3D PDE Anisotropic smoothing. (各向异性滤波)
	Arguments:
	{
	iter:			The iteration number of smoothing.
	K:				Control the smoothing level.
	lamdba:			The step lengh of TV.
	}
	Return:				Status.
	**********************************************************/
	bool PDEAnisotropic3D(int iter, float K, float lamda)const;//3D各项异性扩散滤波
	/**********************************************************
	Description:		3D mean smoothing
	Arguments:
	{
	windowHalf:		The half of the filtering window.
	}
	Return:				Status.
	**********************************************************/
	bool MeanDenoising3D(int windowHalf)const;//均值滤波
	/**********************************************************
	Description:		3D median smoothing
	Arguments:
	{
	windowHalf:		The half of the filtering window.
	}
	Return:				Status.
	**********************************************************/
	bool MedianDenoising3D(int windowHalf)const;//中值滤波
	/**********************************************************
	Description:		3D Guassian smoothing, and each direction is independent.
	Arguments:
	{
	windowHalf:		The half of the filtering window.
	FWHM:			Full width at half maximum of the Guassian kernel.
	}
	Return:				Status.
	**********************************************************/
	bool GuassianDenoising3D(int windowHalf, float FWHM)const;//高斯滤波
	/**********************************************************
	Description:		3D filtering for the Negative pixel.
	Arguments:
	{
	windowHalf:		The half of the filtering window.
	}
	Return:				Status.
	**********************************************************/
	bool NegativeFilter3D(int windowHalf = 1)const;//将负值体素等于领域非负值体素的平均
	/**********************************************************
	Description:		Set the pixel out of the FOV to 0.
	Arguments:
	{

	}
	Return:				Status.
	**********************************************************/
	bool SetOutOfViewZero()const;

public:
	float	*m_pData;//x+y*nX+z*nX*nY
private:
	UImage(const UImage&);
private:
	int		m_nVoxelNumInAll;
	UVolume3D<double>	voxelSize;	//单位mm
	UVolume3D<double>	center;
	UVolume3D<int>		voxelNum;
};

#endif
