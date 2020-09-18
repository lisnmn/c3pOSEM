/********************************************************************

Copyright (C), 2019, All rights reserved

File Name     :    URaytracing.cpp
Description   :
History       :

<author>            <time>            <desc>
Ang Li             2019/6/8           create

********************************************************************/
#include "URayTracing.h"
#include "UOneLor.h"
#include "UImage.h"
#include <math.h>
#include <algorithm>

#include <iostream>

using namespace std;
#define ZERO_PLUS 0.0001
#define MAX_VOXEL_NUM 4096
//#define max(x,y)	((x)<(y)?(y):(x)) //using the max,min function from <algorithm>
//#define min(x,y)	((x)<(y)?(x):(y))
//template<class T>
//inline T max(T x, T y){ return ((x) < (y) ? (y) : (x)); }
//template<class T>
//inline T min(T x, T y){ return ((x) < (y) ? (x) : (y)); }



URayTracing::URayTracing()
{
	m_aAll = new double[3 * MAX_VOXEL_NUM];
	m_ax = new double[MAX_VOXEL_NUM];
	m_ay = new double[MAX_VOXEL_NUM];
	m_az = new double[MAX_VOXEL_NUM];
}


URayTracing::~URayTracing()
{
	delete[]m_aAll;
	delete[]m_ax;
	delete[]m_ay;
	delete[]m_az;
}

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
void URayTracing::GetOneLorResponse3D(UOneLor &ray, const UImage & img, UVolume3D<double>& point1, UVolume3D<double>& point2)
{
	double ax0, axn;//光线与体素各方向边界平面交点的a（比例）值
	double ay0, ayn;
	double az0, azn;
	double amax, amin;//a的最大最小值
	UVolume3D<int> voxelNum = img.GetVoxelNum();
	UVolume3D<double> voxCenter = img.GetVoxelCenter();
	UVolume3D<double> voxelSize = img.GetVoxelSize();
	//体素各轴向的边界坐标值
	double xPlane_Nx, xPlane_0;
	double yPlane_Ny, yPlane_0;
	double zPlane_Nz, zPlane_0;
	xPlane_Nx = voxCenter.x + (double)voxelNum.x / 2 * voxelSize.x;
	xPlane_0 = voxCenter.x - (double)voxelNum.x / 2 * voxelSize.x;

	yPlane_0 = voxCenter.y - (double)voxelNum.y / 2 * voxelSize.y;
	yPlane_Ny = voxCenter.y + (double)voxelNum.y / 2 * voxelSize.y;

	zPlane_0 = voxCenter.z - (double)voxelNum.z / 2 * voxelSize.z;
	zPlane_Nz = voxCenter.z + (double)voxelNum.z / 2 * voxelSize.z;
	//判断首位两点是否存在xyz方向上的两点坐标相同
	if (fabs(point1.x - point2.x) < ZERO_PLUS && fabs(point1.z - point2.z) < ZERO_PLUS)//x1==x2 && z1==z2
	{
		if (point1.x > xPlane_Nx || point1.x<xPlane_0 || point1.z>zPlane_Nz || point1.z < zPlane_0)
		{
			ray.SetVoxelNum(0);
			return;
		}
		//简化，懒得算起始点和终止点了
		ray.SetVoxelNum(voxelNum.y);
		int xID, yID, zID;
		xID = static_cast<int>((point1.x - xPlane_0) / voxelSize.x);
		zID = static_cast<int>((point1.z - zPlane_0) / voxelSize.z);
		for (int i = 0; i < ray.GetVoxelNum(); i++)
		{
			ray.pWeight[i] = static_cast<float>(voxelSize.y);
			yID = i;
			ray.pVoxelID[i] = xID + yID * voxelNum.x + zID * voxelNum.x*voxelNum.y;
		}
	}
	else if (fabs(point1.y - point2.y) < ZERO_PLUS && fabs(point1.z - point2.z) < ZERO_PLUS)//y1==y2 && z1==z2
	{
		if (point1.y > yPlane_Ny || point1.y<yPlane_0 || point1.z>zPlane_Nz || point1.z < zPlane_0)
		{
			ray.SetVoxelNum(0);
			return;
		}
		//简化，懒得算起始点和终止点了
		ray.SetVoxelNum(voxelNum.x);
		int xID, yID, zID;
		yID = static_cast<int>((point1.y - yPlane_0) / voxelSize.y);
		zID = static_cast<int>((point1.z - zPlane_0) / voxelSize.z);
		for (int i = 0; i < ray.GetVoxelNum(); i++)
		{
			ray.pWeight[i] = static_cast<float>(voxelSize.x);
			xID = i;
			ray.pVoxelID[i] = xID + yID * voxelNum.x + zID * voxelNum.x*voxelNum.y;
		}
	}
	else if (fabs(point1.x - point2.x) < ZERO_PLUS)//x1==x2
	{
		if (point1.x >= xPlane_Nx || point1.x <= xPlane_0)
		{
			ray.SetVoxelNum(0);
			return;
		}
		//calculate amax amin
		ay0 = (yPlane_0 - point1.y) / (point2.y - point1.y);
		ayn = (yPlane_Ny - point1.y) / (point2.y - point1.y);

		az0 = (zPlane_0 - point1.z) / (point2.z - point1.z);
		azn = (zPlane_Nz - point1.z) / (point2.z - point1.z);

		amax = min(1.0, min(max(ay0, ayn), max(az0, azn)));
		amin = max(0.0, max(min(ay0, ayn), min(az0, azn)));
		//amax<amin时未穿过体素
		if (amax - amin < ZERO_PLUS)
		{
			ray.SetVoxelNum(0);
			return;
		}
		//求各方向体素需要求解的范围
		int ymin, ymax, zmin, zmax;
		double yPlane_i, zPlane_i;
		//计算ax,ay,az的取值范围和范围内的值，需要分X2>X1和X2<X1的情况
		if (point2.y - point1.y > 0)//y
		{
			ymin = static_cast<int>( 1 + voxelNum.y - (yPlane_Ny - amin * (point2.y - point1.y) - point1.y) / voxelSize.y );
			ymax = static_cast<int>( (point1.y + amax * (point2.y - point1.y) - yPlane_0) / voxelSize.y );
			//ay = new double[ymax - ymin + 1];
			for (int i = 0; i < ymax - ymin + 1; i++)
			{//{ay}={ay[ymin],...,ay[ymax]}
				yPlane_i = (ymin + i - (double)voxelNum.y / 2)*voxelSize.y + voxCenter.y;
				m_ay[i] = (yPlane_i - point1.y) / (point2.y - point1.y);
			}
		}
		else
		{
			ymin = static_cast<int>( 1 + voxelNum.y - (yPlane_Ny - amax * (point2.y - point1.y) - point1.y) / voxelSize.y );
			ymax = static_cast<int>( (point1.y + amin * (point2.y - point1.y) - yPlane_0) / voxelSize.y );
			//ay = new double[ymax - ymin + 1];
			for (int i = 0; i < ymax - ymin + 1; i++)
			{//{ay}={ay[ymax],...,ax[ymin]}
				yPlane_i = (ymin + i - (double)voxelNum.y / 2)*voxelSize.y + voxCenter.y;
				m_ay[ymax - ymin - i] = (yPlane_i - point1.y) / (point2.y - point1.y);
			}
		}

		if (point2.z - point1.z > 0)//z
		{
			zmin = static_cast<int>( 1 + voxelNum.z - (zPlane_Nz - amin * (point2.z - point1.z) - point1.z) / voxelSize.z );
			zmax = static_cast<int>( (point1.z + amax * (point2.z - point1.z) - zPlane_0) / voxelSize.z );
			//az = new double[zmax - zmin + 1];
			for (int i = 0; i < zmax - zmin + 1; i++)
			{//{az}={az[zmax],...,az[zmin]}
				zPlane_i = (zmin + i - (double)voxelNum.z / 2)*voxelSize.z + voxCenter.z;
				m_az[i] = (zPlane_i - point1.z) / (point2.z - point1.z);
			}
		}
		else
		{
			zmin = static_cast<int>( 1 + voxelNum.z - (zPlane_Nz - amax * (point2.z - point1.z) - point1.z) / voxelSize.z );
			zmax = static_cast<int>( (point1.z + amin * (point2.z - point1.z) - zPlane_0) / voxelSize.z );
			//az = new double[zmax - zmin + 1];
			for (int i = 0; i < zmax - zmin + 1; i++)
			{//{az}={az[zmax],...,ax[zmin]}
				zPlane_i = (zmin + i - (double)voxelNum.z / 2)*voxelSize.z + voxCenter.z;
				m_az[zmax - zmin - i] = (zPlane_i - point1.z) / (point2.z - point1.z);
			}
		}
		//合并{amin,{ax,ay,az},amax}和排序
		int  yCount, zCount;
		yCount = (ymax - ymin + 1);
		zCount = (zmax - zmin + 1);
		int aCount = 2 + yCount + zCount;
		m_aAll[0] = amin;
		m_aAll[aCount - 1] = amax;
		for (int i = 0; i < yCount; i++)
		{
			m_aAll[1 + i] = m_ay[i];
		}
		for (int i = 0; i < zCount; i++)
		{
			m_aAll[1 + yCount + i] = m_az[i];
		}
		sort(m_aAll, m_aAll + aCount);
		//求线长和体素编号，有可能存在线长为0的体素
		UVolume3D<double> vctL(point2 - point1);//用于求首位两点长度的向量
		double L = sqrt(vctL*vctL);//首位两点的长度
		ray.SetVoxelNum(aCount - 3);//丢掉两边的两个，最大的边有可能出界
		int xID, yID, zID;
		xID = static_cast<int>( (point1.x - xPlane_0) / voxelSize.x );
		for (int i = 0; i < ray.GetVoxelNum(); i++)
		{
			ray.pWeight[i] = static_cast<float>( L * (m_aAll[i + 2] - m_aAll[i + 1]) );
			yID = static_cast<int>( (point1.y + (m_aAll[i + 2] + m_aAll[i + 1]) / 2 * (point2.y - point1.y) - yPlane_0) / voxelSize.y );//与文献上不同，不+1
			zID = static_cast<int>( (point1.z + (m_aAll[i + 2] + m_aAll[i + 1]) / 2 * (point2.z - point1.z) - zPlane_0) / voxelSize.z );
			ray.pVoxelID[i] = xID + yID * voxelNum.x + zID * voxelNum.x*voxelNum.y;
		}
	}
	else if (fabs(point1.y - point2.y) < ZERO_PLUS)//y1==y2
	{
		if (point1.y >= yPlane_Ny || point1.y <= yPlane_0)
		{
			ray.SetVoxelNum(0);
			return;
		}
		//calculate amax amin
		ax0 = (xPlane_0 - point1.x) / (point2.x - point1.x);
		axn = (xPlane_Nx - point1.x) / (point2.x - point1.x);

		az0 = (zPlane_0 - point1.z) / (point2.z - point1.z);
		azn = (zPlane_Nz - point1.z) / (point2.z - point1.z);

		amax = min(min(1.0, max(ax0, axn)), max(az0, azn));
		amin = max(max(0.0, min(ax0, axn)), min(az0, azn));
		//amax<amin时未穿过体素
		if (amax - amin < ZERO_PLUS)
		{
			ray.SetVoxelNum(0);
			return;
		}
		//求各方向体素需要求解的范围
		int xmin, xmax, zmin, zmax;
		double xPlane_i, zPlane_i;
		//计算ax,ay,az的取值范围和范围内的值，需要分X2>X1和X2<X1的情况
		if (point2.x - point1.x > 0)//x
		{
			xmin = static_cast<int>( 1 + voxelNum.x - (xPlane_Nx - amin * (point2.x - point1.x) - point1.x) / voxelSize.x );
			xmax = static_cast<int>( (point1.x + amax * (point2.x - point1.x) - xPlane_0) / voxelSize.x );//+1的位置与文献有一点不同，和从0计数有关
			//ax = new double[xmax - xmin + 1];
			for (int i = 0; i < xmax - xmin + 1; i++)
			{//{ax}={ax[xmin],...,ax[xmax]}
				xPlane_i = (xmin + i - (double)voxelNum.x / 2)*voxelSize.x + voxCenter.x;//要计算体素的中心偏移量
				m_ax[i] = (xPlane_i - point1.x) / (point2.x - point1.x);
			}
		}
		else
		{
			xmin = static_cast<int>( 1 + voxelNum.x - (xPlane_Nx - amax * (point2.x - point1.x) - point1.x) / voxelSize.x );
			xmax = static_cast<int>( (point1.x + amin * (point2.x - point1.x) - xPlane_0) / voxelSize.x );
			//ax = new double[xmax - xmin + 1];
			for (int i = 0; i < xmax - xmin + 1; i++)
			{//{ax}={ax[xmax],...,ax[xmin]}
				xPlane_i = (xmin + i - (double)voxelNum.x / 2)*voxelSize.x + voxCenter.x;
				m_ax[xmax - xmin - i] = (xPlane_i - point1.x) / (point2.x - point1.x);
			}
		}

		if (point2.z - point1.z > 0)//z
		{
			zmin = static_cast<int>( 1 + voxelNum.z - (zPlane_Nz - amin * (point2.z - point1.z) - point1.z) / voxelSize.z );
			zmax = static_cast<int>( (point1.z + amax * (point2.z - point1.z) - zPlane_0) / voxelSize.z );
			//az = new double[zmax - zmin + 1];
			for (int i = 0; i < zmax - zmin + 1; i++)
			{//{az}={az[zmax],...,az[zmin]}
				zPlane_i = (zmin + i - (double)voxelNum.z / 2)*voxelSize.z + voxCenter.z;
				m_az[i] = (zPlane_i - point1.z) / (point2.z - point1.z);
			}
		}
		else
		{
			zmin = static_cast<int>( 1 + voxelNum.z - (zPlane_Nz - amax * (point2.z - point1.z) - point1.z) / voxelSize.z );
			zmax = static_cast<int>( (point1.z + amin * (point2.z - point1.z) - zPlane_0) / voxelSize.z );
			//az = new double[zmax - zmin + 1];
			for (int i = 0; i < zmax - zmin + 1; i++)
			{//{az}={az[zmax],...,ax[zmin]}
				zPlane_i = (zmin + i - (double)voxelNum.z / 2)*voxelSize.z + voxCenter.z;
				m_az[zmax - zmin - i] = (zPlane_i - point1.z) / (point2.z - point1.z);
			}
		}
		//合并{amin,{ax,ay,az},amax}和排序
		int xCount, zCount;
		xCount = (xmax - xmin + 1);
		zCount = (zmax - zmin + 1);
		int aCount = 2 + xCount + zCount;
		m_aAll[0] = amin;
		m_aAll[aCount - 1] = amax;
		for (int i = 0; i < xCount; i++)
		{
			m_aAll[1 + i] = m_ax[i];
		}
		for (int i = 0; i < zCount; i++)
		{
			m_aAll[1 + xCount + i] = m_az[i];
		}
		sort(m_aAll, m_aAll + aCount);
		//求线长和体素编号，有可能存在线长为0的体素
		UVolume3D<double> vctL(point2 - point1);//用于求首位两点长度的向量
		double L = sqrt(vctL*vctL);//首位两点的长度
		ray.SetVoxelNum(aCount - 3);//丢掉两边的两个，最大的边有可能出界
		int xID, yID, zID;
		yID = static_cast<int>( (point1.y - yPlane_0) / voxelSize.y );
		for (int i = 0; i < ray.GetVoxelNum(); i++)
		{
			ray.pWeight[i] = static_cast<float>( L * (m_aAll[i + 2] - m_aAll[i + 1]) );
			xID = static_cast<int>( (point1.x + (m_aAll[i + 2] + m_aAll[i + 1]) / 2 * (point2.x - point1.x) - xPlane_0) / voxelSize.x );//与文献上不同，不+1
			zID = static_cast<int>( (point1.z + (m_aAll[i + 2] + m_aAll[i + 1]) / 2 * (point2.z - point1.z) - zPlane_0) / voxelSize.z );
			ray.pVoxelID[i] = xID + yID * voxelNum.x + zID * voxelNum.x*voxelNum.y;
		}
	}
	else if (fabs(point1.z - point2.z) < ZERO_PLUS)//z1==z2
	{
		if (point1.z > zPlane_Nz || point1.z < zPlane_0)
		{
			ray.SetVoxelNum(0);
			return;
		}
		//calculate amax amin
		ax0 = (xPlane_0 - point1.x) / (point2.x - point1.x);
		axn = (xPlane_Nx - point1.x) / (point2.x - point1.x);

		ay0 = (yPlane_0 - point1.y) / (point2.y - point1.y);
		ayn = (yPlane_Ny - point1.y) / (point2.y - point1.y);

		amax = min(min(1.0, max(ax0, axn)), max(ay0, ayn));
		amin = max(max(0.0, min(ax0, axn)), min(ay0, ayn));
		//amax<amin时未穿过体素
		if (amax - amin < ZERO_PLUS)
		{
			ray.SetVoxelNum(0);
			return;
		}
		//求各方向体素需要求解的范围
		int xmin, xmax, ymin, ymax;
		double xPlane_i, yPlane_i;
		//计算ax,ay的取值范围和范围内的值，需要分X2>X1和X2<X1的情况
		if (point2.x - point1.x > 0)//x
		{
			xmin = static_cast<int>( 1 + voxelNum.x - (xPlane_Nx - amin * (point2.x - point1.x) - point1.x) / voxelSize.x );
			xmax = static_cast<int>( (point1.x + amax * (point2.x - point1.x) - xPlane_0) / voxelSize.x );//+1的位置与文献有一点不同，和从0计数有关
			//ax = new double[xmax - xmin + 1];
			for (int i = 0; i < xmax - xmin + 1; i++)
			{//{ax}={ax[xmin],...,ax[xmax]}
				xPlane_i = (xmin + i - (double)voxelNum.x / 2)*voxelSize.x + voxCenter.x;//要计算体素的中心偏移量
				m_ax[i] = (xPlane_i - point1.x) / (point2.x - point1.x);
			}
		}
		else
		{
			xmin = static_cast<int>( 1 + voxelNum.x - (xPlane_Nx - amax * (point2.x - point1.x) - point1.x) / voxelSize.x );
			xmax = static_cast<int>( (point1.x + amin * (point2.x - point1.x) - xPlane_0) / voxelSize.x );
			//ax = new double[xmax - xmin + 1];
			for (int i = 0; i < xmax - xmin + 1; i++)
			{//{ax}={ax[xmax],...,ax[xmin]}
				xPlane_i = (xmin + i - (double)voxelNum.x / 2)*voxelSize.x + voxCenter.x;
				m_ax[xmax - xmin - i] = (xPlane_i - point1.x) / (point2.x - point1.x);
			}
		}

		if (point2.y - point1.y > 0)//y
		{
			ymin = static_cast<int>( 1 + voxelNum.y - (yPlane_Ny - amin * (point2.y - point1.y) - point1.y) / voxelSize.y );
			ymax = static_cast<int>( (point1.y + amax * (point2.y - point1.y) - yPlane_0) / voxelSize.y );
			//ay = new double[ymax - ymin + 1];
			for (int i = 0; i < ymax - ymin + 1; i++)
			{//{ay}={ay[ymin],...,ay[ymax]}
				yPlane_i = (ymin + i - (double)voxelNum.y / 2)*voxelSize.y + voxCenter.y;
				m_ay[i] = (yPlane_i - point1.y) / (point2.y - point1.y);
			}
		}
		else
		{
			ymin = static_cast<int>( 1 + voxelNum.y - (yPlane_Ny - amax * (point2.y - point1.y) - point1.y) / voxelSize.y );
			ymax = static_cast<int>( (point1.y + amin * (point2.y - point1.y) - yPlane_0) / voxelSize.y );
			//ay = new double[ymax - ymin + 1];
			for (int i = 0; i < ymax - ymin + 1; i++)
			{//{ay}={ay[ymax],...,ax[ymin]}
				yPlane_i = (ymin + i - (double)voxelNum.y / 2)*voxelSize.y + voxCenter.y;
				m_ay[ymax - ymin - i] = (yPlane_i - point1.y) / (point2.y - point1.y);
			}
		}
		//合并{amin,{ax,ay,az},amax}和排序
		int xCount, yCount;// , zCount;
		xCount = (xmax - xmin + 1);
		yCount = (ymax - ymin + 1);
		int aCount = 2 + xCount + yCount;
		m_aAll[0] = amin;
		m_aAll[aCount - 1] = amax;
		for (int i = 0; i < xCount; i++)
		{
			m_aAll[1 + i] = m_ax[i];
		}
		for (int i = 0; i < yCount; i++)
		{
			m_aAll[1 + xCount + i] = m_ay[i];
		}
		sort(m_aAll, m_aAll + aCount);
		//求线长和体素编号，有可能存在线长为0的体素
		UVolume3D<double> vctL(point2 - point1);//用于求首位两点长度的向量
		double L = sqrt(vctL*vctL);//首位两点的长度
		ray.SetVoxelNum(aCount - 3);//丢掉两边的两个，最大的边有可能出界
		int xID, yID, zID;
		zID = static_cast<int>( (point1.z - zPlane_0) / voxelSize.z );
		for (int i = 0; i < ray.GetVoxelNum(); i++)
		{
			ray.pWeight[i] = static_cast<float>( L * (m_aAll[i + 2] - m_aAll[i + 1]) );
			xID = static_cast<int>( (point1.x + (m_aAll[i + 2] + m_aAll[i + 1]) / 2 * (point2.x - point1.x) - xPlane_0) / voxelSize.x );//与文献上不同，不+1
			yID = static_cast<int>( (point1.y + (m_aAll[i + 2] + m_aAll[i + 1]) / 2 * (point2.y - point1.y) - yPlane_0) / voxelSize.y );
			ray.pVoxelID[i] = xID + yID * voxelNum.x + zID * voxelNum.x*voxelNum.y;
		}
	}
	else//normal
	{
		//calculate amax amin
		ax0 = (xPlane_0 - point1.x) / (point2.x - point1.x);
		axn = (xPlane_Nx - point1.x) / (point2.x - point1.x);

		ay0 = (yPlane_0 - point1.y) / (point2.y - point1.y);
		ayn = (yPlane_Ny - point1.y) / (point2.y - point1.y);

		az0 = (zPlane_0 - point1.z) / (point2.z - point1.z);
		azn = (zPlane_Nz - point1.z) / (point2.z - point1.z);

		amax = min(min(1.0, max(ax0, axn)), min(max(ay0, ayn), max(az0, azn)));
		amin = max(max(0.0, min(ax0, axn)), max(min(ay0, ayn), min(az0, azn)));
		//amax<amin时未穿过体素
		if (amax - amin < ZERO_PLUS)
		{
			ray.SetVoxelNum(0);
			return;
		}
		//求各方向体素需要求解的范围
		int xmin, xmax, ymin, ymax, zmin, zmax;
		double xPlane_i, yPlane_i, zPlane_i;
		//计算ax,ay,az的取值范围和范围内的值，需要分X2>X1和X2<X1的情况
		if (point2.x - point1.x > 0)//x
		{
			xmin = static_cast<int>( 1 + voxelNum.x - (xPlane_Nx - amin * (point2.x - point1.x) - point1.x) / voxelSize.x );
			xmax = static_cast<int>( (point1.x + amax * (point2.x - point1.x) - xPlane_0) / voxelSize.x );//+1的位置与文献有一点不同，和从0计数有关
			//ax = new double[xmax - xmin + 1];
			for (int i = 0; i < xmax - xmin + 1; i++)
			{//{ax}={ax[xmin],...,ax[xmax]}
				xPlane_i = (xmin + i - (double)voxelNum.x / 2)*voxelSize.x + voxCenter.x;//要计算体素的中心偏移量
				m_ax[i] = (xPlane_i - point1.x) / (point2.x - point1.x);
			}
		}
		else
		{
			xmin = static_cast<int>( 1 + voxelNum.x - (xPlane_Nx - amax * (point2.x - point1.x) - point1.x) / voxelSize.x );
			xmax = static_cast<int>( (point1.x + amin * (point2.x - point1.x) - xPlane_0) / voxelSize.x );
			//ax = new double[xmax - xmin + 1];
			for (int i = 0; i < xmax - xmin + 1; i++)
			{//{ax}={ax[xmax],...,ax[xmin]}
				xPlane_i = (xmin + i - (double)voxelNum.x / 2)*voxelSize.x + voxCenter.x;
				m_ax[xmax - xmin - i] = (xPlane_i - point1.x) / (point2.x - point1.x);
			}
		}

		if (point2.y - point1.y > 0)//y
		{
			ymin = static_cast<int>( 1 + voxelNum.y - (yPlane_Ny - amin * (point2.y - point1.y) - point1.y) / voxelSize.y );
			ymax = static_cast<int>( (point1.y + amax * (point2.y - point1.y) - yPlane_0) / voxelSize.y );
			//ay = new double[ymax - ymin + 1];
			for (int i = 0; i < ymax - ymin + 1; i++)
			{//{ay}={ay[ymin],...,ay[ymax]}
				yPlane_i = (ymin + i - (double)voxelNum.y / 2)*voxelSize.y + voxCenter.y;
				m_ay[i] = (yPlane_i - point1.y) / (point2.y - point1.y);
			}
		}
		else
		{
			ymin = static_cast<int>( 1 + voxelNum.y - (yPlane_Ny - amax * (point2.y - point1.y) - point1.y) / voxelSize.y );
			ymax = static_cast<int>( (point1.y + amin * (point2.y - point1.y) - yPlane_0) / voxelSize.y );
			//ay = new double[ymax - ymin + 1];
			for (int i = 0; i < ymax - ymin + 1; i++)
			{//{ay}={ay[ymax],...,ax[ymin]}
				yPlane_i = (ymin + i - (double)voxelNum.y / 2)*voxelSize.y + voxCenter.y;
				m_ay[ymax - ymin - i] = (yPlane_i - point1.y) / (point2.y - point1.y);
			}
		}

		if (point2.z - point1.z > 0)//z
		{
			zmin = static_cast<int>( 1 + voxelNum.z - (zPlane_Nz - amin * (point2.z - point1.z) - point1.z) / voxelSize.z );
			zmax = static_cast<int>( (point1.z + amax * (point2.z - point1.z) - zPlane_0) / voxelSize.z );
			//az = new double[zmax - zmin + 1];
			for (int i = 0; i < zmax - zmin + 1; i++)
			{//{az}={az[zmax],...,az[zmin]}
				zPlane_i = (zmin + i - (double)voxelNum.z / 2)*voxelSize.z + voxCenter.z;
				m_az[i] = (zPlane_i - point1.z) / (point2.z - point1.z);
			}
		}
		else
		{
			zmin = static_cast<int>( 1 + voxelNum.z - (zPlane_Nz - amax * (point2.z - point1.z) - point1.z) / voxelSize.z );
			zmax = static_cast<int>( (point1.z + amin * (point2.z - point1.z) - zPlane_0) / voxelSize.z );
			//az = new double[zmax - zmin + 1];
			for (int i = 0; i < zmax - zmin + 1; i++)
			{//{az}={az[zmax],...,ax[zmin]}
				zPlane_i = (zmin + i - (double)voxelNum.z / 2)*voxelSize.z + voxCenter.z;
				m_az[zmax - zmin - i] = (zPlane_i - point1.z) / (point2.z - point1.z);
			}
		}
		//合并{amin,{ax,ay,az},amax}和排序
		int xCount, yCount, zCount;
		xCount = (xmax - xmin + 1);
		yCount = (ymax - ymin + 1);
		zCount = (zmax - zmin + 1);
		int aCount = 2 + xCount + yCount + zCount;
		m_aAll[0] = amin;
		m_aAll[aCount - 1] = amax;
		for (int i = 0; i < xCount; i++)
		{
			m_aAll[1 + i] = m_ax[i];
		}
		for (int i = 0; i < yCount; i++)
		{
			m_aAll[1 + xCount + i] = m_ay[i];
		}
		for (int i = 0; i < zCount; i++)
		{
			m_aAll[1 + xCount + yCount + i] = m_az[i];
		}
		sort(m_aAll, m_aAll + aCount);
		//求线长和体素编号，有可能存在线长为0的体素
		UVolume3D<double> vctL(point2 - point1);//用于求首位两点长度的向量
		double L = sqrt(vctL*vctL);//首位两点的长度
		ray.SetVoxelNum(aCount - 3);//丢掉两边的两个，最大的边有可能出界
		int xID, yID, zID;
		double aMid;
		for (int i = 0; i < ray.GetVoxelNum(); i++)
		{
			aMid = (m_aAll[i + 2] + m_aAll[i + 1]) / 2;
			ray.pWeight[i] = static_cast<float>( L * (m_aAll[i + 2] - m_aAll[i + 1]) );
			xID = static_cast<int>( (point1.x + aMid * (point2.x - point1.x) - xPlane_0) / voxelSize.x );//与文献上不同，不+1
			yID = static_cast<int>( (point1.y + aMid * (point2.y - point1.y) - yPlane_0) / voxelSize.y );
			zID = static_cast<int>( (point1.z + aMid * (point2.z - point1.z) - zPlane_0) / voxelSize.z );
			ray.pVoxelID[i] = xID + yID * voxelNum.x + zID * voxelNum.x*voxelNum.y;
		}
	}
	int nVoxelNum = voxelNum.x*voxelNum.y*voxelNum.z;//为了防止voxelID超出范围//20190330
	for (int i = 0; i < ray.GetVoxelNum(); i++)
	{
		if (ray.pVoxelID[i] < 0 || ray.pVoxelID[i] >= nVoxelNum)
		{
			ray.pWeight[i] = 0;
			ray.pVoxelID[i] = 0;
		}
	}
}
