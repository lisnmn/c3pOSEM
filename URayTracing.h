/********************************************************************

Copyright (C), 2019, All rights reserved

File Name     :    URayTracing.cpp
Description   :
History       :

<author>            <time>            <desc>
Ang Li             2019/6/8           create

********************************************************************/
/****************************************************************************************************
 * 此代码计算方法的参考文献为:
 * "Fast calculation of the exact radiological path for a three-dimensional CT array"
*****************************************************************************************************/
#ifndef __RAYTRACING_H
#define __RAYTRACING_H

#include "UVolume.h"

class UImage;
class UOneLor;
/*********************************************************************************
 * 提供计算一条光线穿过图像各体素的索引和长度的函数
 * 20190125优化：优化开辟内存耗费的时间
 * 注意：该类需要防止多线程时，子体素写入时的竞争错误。
 * 建议：每个线程分别建立该类的实例
**********************************************************************************/
class URayTracing
{
public:
	URayTracing();
	~URayTracing();
	/**********************************************************
	Description:		Calculate a 3D raytracing for one ray.
	Arguments:
	{
	ray:			Return the lenth of the ray across the each voxel.
	img:			The calculated image object.
	point1:			One end of the ray.
	point2:			Another end of the ray.
	}
	Return:				Status.
	**********************************************************/
	void GetOneLorResponse3D(UOneLor &ray, const UImage &img, UVolume3D<double> & point1, UVolume3D<double> &point2);
private:
	URayTracing(const URayTracing&){};
	URayTracing& operator=(const URayTracing&){};
private:
	//用于存放三个方向的比例系数的数组，为了减少开辟释放的时间，在类初始化和结束时开辟和释放这些数组的内存。
	double *m_ax;
	double *m_ay;
	double *m_az;
	double *m_aAll;
};

#endif
