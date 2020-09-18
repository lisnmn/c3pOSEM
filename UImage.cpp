/********************************************************************

Copyright (C), 2019, All rights reserved

File Name     :    UImage.cpp
Description   :
History       :

<author>            <time>            <desc>
Ang Li             2019/6/8           create

********************************************************************/
#include "UImage.h"

#include <math.h>
#include <time.h>

//#define max(x,y)	((x)<(y)?(y):(x))
//#define min(x,y)	((x)<(y)?(x):(y))
template<class T>
inline T max(T x, T y){ return ((x) < (y) ? (y) : (x)); }
template<class T>
inline T min(T x, T y){ return ((x) < (y) ? (x) : (y)); }


UImage::UImage() :
	m_pData(nullptr), m_nVoxelNumInAll(0)
{
}


UImage::~UImage()
{
	if (m_pData != nullptr)
		delete[]m_pData;
	m_pData = nullptr;
}

UImage& UImage::operator=(const UImage& img)
{
	if (this != &img)
	{
		float *temp = nullptr;
		temp = new float[img.m_nVoxelNumInAll];
		for (int i = 0; i < img.m_nVoxelNumInAll; i++)
		{
			temp[i] = img.m_pData[i];
		}
		this->center = img.center;
		this->voxelNum = img.voxelNum;
		this->voxelSize = img.voxelSize;
		this->m_nVoxelNumInAll = img.m_nVoxelNumInAll;
		if (m_pData != nullptr)
			delete[]m_pData;
		m_pData = temp;
	}
	return *this;
}

int UImage::GetVoxelNumX() const
{
	return voxelNum.x;
}

int UImage::GetVoxelNumY() const
{
	return voxelNum.y;
}

int UImage::GetVoxelNumZ() const
{
	return voxelNum.z;
}

int UImage::GetVoxelNumInAll() const
{
	return m_nVoxelNumInAll;;
}

UVolume3D<int> UImage::GetVoxelNum() const
{
	return voxelNum;
}

UVolume3D<double> UImage::GetVoxelSize() const
{
	return voxelSize;
}

UVolume3D<double> UImage::GetVoxelCenter() const
{
	return center;
}

void UImage::SetCenter(double x, double y, double z)
{
	center.set(x, y, z);
}

bool UImage::SetVoxelSize(double x, double y, double z)
{
	if (x > 0 && y > 0 && z > 0)
	{
		voxelSize.set(x, y, z);
		return true;
	}
	else
		return false;
}
/**********************************************************
Description:		It will initial the memory when you set the voxel number.
Arguments:			Voxel number in 3 dirctions.
Return:				True when it runs succcessfully.
**********************************************************/
bool UImage::SetVoxelNum(int x, int y, int z)
{
	if (x > 0 && y > 0 && z > 0)
	{
		voxelNum.set(x, y, z);
		m_nVoxelNumInAll = voxelNum.x*voxelNum.y*voxelNum.z;
		if (m_pData != nullptr)
		{
			delete[]m_pData;
			m_pData = nullptr;
		}
		m_pData = new float[m_nVoxelNumInAll];
		if (m_pData != nullptr)
			return true;
		else
			return false;
		
	}
	else
		return false;
}
/**********************************************************
Description:		It will initial the 3D image to a cylinder.
Arguments:
{
	pos:			The center position of the cylinder.
	r:				The radius of the cylinder.
	h:				The hight of the cylinder.
	value:			every pixel in the cylinder will be set to 'value'
}
Return:				void.
**********************************************************/
void UImage::SetCylinder(UVolume3D<double> pos, double r, double h, float value)
{
	for (int i = 0; i < m_nVoxelNumInAll; i++)
		m_pData[i] = 0;
	double xx, yy, zz;//体素的中心坐标
	for (int z = 0; z < voxelNum.z; z++)
	{
		for (int x = 0; x < voxelNum.x; x++)
		{
			for (int y = 0; y < voxelNum.y; y++)
			{
				xx = ((double)x - (double)voxelNum.x / 2 + 0.5) * voxelSize.x + center.z;//体素的中心坐标
				yy = ((double)y - (double)voxelNum.y / 2 + 0.5) * voxelSize.y + center.y;
				if ((xx - pos.x)*(xx - pos.x) + (yy - pos.y)*(yy - pos.y) < r*r)
				{
					zz = ((double)z - (double)voxelNum.z / 2 + 0.5) * voxelSize.z + center.z;
					if (fabs(zz - pos.z) < h / 2)
						m_pData[x + y * voxelNum.x + z * voxelNum.x*voxelNum.y] = value;

				}
			}
		}
	}
}
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
void UImage::GetVoxelEdgePos(int index, UVolume3D<double>& position1, UVolume3D<double>& position2) const
{
	int x, y, z;//GATE下坐标系方向像素索引值
	z = index / (voxelNum.x*voxelNum.y);
	y = index % (voxelNum.x*voxelNum.y) / voxelNum.x;
	x = index % voxelNum.x;
	double xx, yy, zz;//以探测器中心为原点坐标//改成了GATE坐标系
	xx = (x - voxelNum.x / 2)*voxelSize.x + center.x;//20190425日改，之前没有算体素大小
	yy = (y - voxelNum.y / 2)*voxelSize.y + center.y;
	zz = (z - voxelNum.z / 2)*voxelSize.x + center.z;

	position1.x = xx;
	position1.y = yy;
	position1.z = zz;

	position2.x = xx + voxelSize.x;
	position2.y = yy + voxelSize.y;
	position2.z = zz + voxelSize.z;
}
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
void UImage::GetSubvoxelCentre(int index, int subNumOneLine, UVolume3D<double>* centre) const
{
	int x, y, z;//GATE下坐标系方向像素索引值
	z = index / (voxelNum.x*voxelNum.y);
	y = index % (voxelNum.x*voxelNum.y) / voxelNum.x;
	x = index % voxelNum.x;
	double xx, yy, zz;//以探测器中心为原点坐标//GATE坐标系
	xx = ((double)x - (double)voxelNum.x / 2)*voxelSize.x + center.x;//20190425日改，之前没有算体素大小
	yy = ((double)y - (double)voxelNum.y / 2)*voxelSize.y + center.y;
	zz = ((double)z - (double)voxelNum.z / 2)*voxelSize.z + center.z;

	int subNumOneLine2 = subNumOneLine * subNumOneLine;//subNumOneLine平方
	double subNumOneLineReciprocal = 1 / double(subNumOneLine);//subNumOneLine倒数
	for (int iz = 0; iz < subNumOneLine; iz++)
		for (int iy = 0; iy < subNumOneLine; iy++)
			for (int ix = 0; ix < subNumOneLine; ix++)
			{
				centre[ix + iy * subNumOneLine + iz * subNumOneLine2].x = xx + (ix + 0.5) * subNumOneLineReciprocal*voxelSize.x;//20190425日改，之前没有算体素大小
				centre[ix + iy * subNumOneLine + iz * subNumOneLine2].y = yy + (iy + 0.5) * subNumOneLineReciprocal*voxelSize.y;
				centre[ix + iy * subNumOneLine + iz * subNumOneLine2].z = zz + (iz + 0.5) * subNumOneLineReciprocal*voxelSize.z;
			}
}
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
bool UImage::TotalVariationDenoising3D(int iter, float dt, float epsilon, float lamdba) const
{
	if (m_pData == nullptr)
	{
		return false;
	}
	if (iter<=0 || dt<=0 || epsilon<=0 || lamdba <=0)
	{
		return false;
	}
	float *pData = m_pData;
	int voxelNumInAll = m_nVoxelNumInAll;

	int yAll = voxelNum.y;
	int xAll = voxelNum.x;
	int zAll = voxelNum.z;

	float *nimg = nullptr, *nimgtemp = nullptr;//updata image & temp image
	nimg = new float[voxelNumInAll];
	nimgtemp = new float[voxelNumInAll];
	if (nimg == nullptr || nimgtemp == nullptr)
		return false;
	for (int i = 0; i < voxelNumInAll; i++)
	{
		nimg[i] = pData[i];
		nimgtemp[i] = pData[i];
	}
	for (int it = 0; it < iter; it++)
	{
#pragma omp parallel for
		for (int z = 1; z < zAll - 1; z++)
		{
			for (int y = 1; y < yAll - 1; y++)
			{
				float fx_p, fy_p, fz_p;
				float fx_n, fx_ny, fx_nz;
				float fy_n, fy_nx, fy_nz;
				float fz_n, fz_ny, fz_nx;
				float TVGradient;
				int id;
				for (int x = 1; x < xAll - 1; x++)
				{
					fx_p = nimg[x + y * xAll + z * xAll*yAll] - nimg[(x - 1) + y * xAll + z * xAll*yAll];
					fy_p = nimg[x + y * xAll + z * xAll*yAll] - nimg[x + (y - 1) * xAll + z * xAll*yAll];
					fz_p = nimg[x + y * xAll + z * xAll*yAll] - nimg[x + y * xAll + (z - 1) * xAll*yAll];

					fx_n = nimg[(x + 1) + y * xAll + z * xAll*yAll] - nimg[x + y * xAll + z * xAll*yAll];
					fx_ny = nimg[(x + 1) + y * xAll + z * xAll*yAll] - nimg[(x + 1) + (y - 1) * xAll + z * xAll*yAll];
					fx_nz = nimg[(x + 1) + y * xAll + z * xAll*yAll] - nimg[(x + 1) + y * xAll + (z - 1) * xAll*yAll];


					fy_n = nimg[x + (y + 1) * xAll + z * xAll*yAll] - nimg[x + y * xAll + z * xAll*yAll];
					fy_nx = nimg[x + (y + 1) * xAll + z * xAll*yAll] - nimg[(x - 1) + (y + 1) * xAll + z * xAll*yAll];
					fy_nz = nimg[x + (y + 1) * xAll + z * xAll*yAll] - nimg[x + (y + 1) * xAll + (z - 1) * xAll*yAll];

					fz_n = nimg[x + y * xAll + (z + 1) * xAll*yAll] - nimg[x + y * xAll + z * xAll*yAll];
					fz_nx = nimg[x + y * xAll + (z + 1) * xAll*yAll] - nimg[(x - 1) + y * xAll + (z + 1) * xAll*yAll];
					fz_ny = nimg[x + y * xAll + (z + 1) * xAll*yAll] - nimg[x + (y - 1) * xAll + (z + 1) * xAll*yAll];

					TVGradient = (fx_p + fy_p + fz_p) / sqrt(fx_p*fx_p + fy_p * fy_p + fz_p * fz_p + epsilon);
					TVGradient -= fx_n / sqrt(fx_n*fx_n + fx_ny * fx_ny + fx_nz * fx_nz + epsilon);
					TVGradient -= fy_n / sqrt(fy_n*fy_n + fy_nx * fy_nx + fy_nz * fy_nz + epsilon);
					TVGradient -= fz_n / sqrt(fz_n*fz_n + fz_nx * fz_nx + fz_ny * fz_ny + epsilon);


					id = x + y * xAll + z * xAll*yAll;
					if (fabs(TVGradient) > pData[id] / dt)//去除奇异点
						TVGradient = pData[id] / dt / 2;
					nimgtemp[id] = nimg[id] - dt * (TVGradient + lamdba * (nimg[id] - pData[id]));
					if (nimgtemp[id] < 0)
						nimgtemp[id] = 0;
				}
			}
		}
		for (int i = 0; i < voxelNumInAll; i++)
			nimg[i] = nimgtemp[i];
	}
	for (int i = 0; i < voxelNumInAll; i++)
		pData[i] = nimg[i];
	delete[]nimg;
	delete[]nimgtemp;
	return true;
}
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
bool UImage::PDEAnisotropic3D(int iter, float K, float lamdba)const
{
	if (m_pData == nullptr)
	{
		return false;
	}
	if (iter<=0 || K<=0 || lamdba<=0)
	{
		return false;
	}
	float *pData = m_pData;
	int voxelNumInAll = m_nVoxelNumInAll;

	int yAll = voxelNum.y;
	int xAll = voxelNum.x;
	int zAll = voxelNum.z;

	float kk = K * K;//K的平方
	//-----------------------------------------
	//开始扩散
	//-----------------------------------------
	for (int it = 0; it < iter; it++)
	{
#pragma omp parallel for
		for (int z = 1; z < zAll - 1; z++)
		{
			float img_xp, c_xp;
			float img_xn, c_xn;
			float img_yp, c_yp;
			float img_yn, c_yn;
			float img_zp, c_zp;
			float img_zn, c_zn;
			int id;

			for (int y = 1; y < yAll - 1; y++)
			{
				for (int x = 1; x < xAll - 1; x++)
				{
					id = x + y * xAll + z * xAll*yAll;

					img_xp = pData[(x + 1) + y * xAll + z * xAll*yAll] - pData[x + y * xAll + z * xAll*yAll];
					img_xn = pData[(x - 1) + y * xAll + z * xAll*yAll] - pData[x + y * xAll + z * xAll*yAll];

					img_yp = pData[x + (y + 1) * xAll + z * xAll*yAll] - pData[x + y * xAll + z * xAll*yAll];
					img_yn = pData[x + (y - 1) * xAll + z * xAll*yAll] - pData[x + y * xAll + z * xAll*yAll];

					img_zp = pData[x + y * xAll + (z + 1) * xAll*yAll] - pData[x + y * xAll + z * xAll*yAll];
					img_zn = pData[x + y * xAll + (z - 1) * xAll*yAll] - pData[x + y * xAll + z * xAll*yAll];

					c_xp = exp(-1 * img_xp*img_xp / kk);
					c_xn = exp(-1 * img_xn*img_xn / kk);

					c_yp = exp(-1 * img_yp*img_yp / kk);
					c_yn = exp(-1 * img_yn*img_yn / kk);

					c_zp = exp(-1 * img_zp*img_zp / kk);
					c_zn = exp(-1 * img_zn*img_zn / kk);

					pData[id] = pData[id] + lamdba * (img_xp * c_xp + img_xn * c_xn + img_yp * c_yp + img_yn * c_yn + img_zp * c_zp + img_zn * c_zn);

				}
			}
		}
	}
	return true;
}
/**********************************************************
Description:		3D mean smoothing
Arguments:
{
	windowHalf:		The half of the filtering window.
}
Return:				Status.
**********************************************************/
bool UImage::MeanDenoising3D(int windowHalf)const
{
	if (m_pData == nullptr)
	{
		return false;
	}
	if (windowHalf < 1)
	{
		return false;
	}
	float *pData = m_pData;
	int voxelNumInAll = m_nVoxelNumInAll;

	int zAll = voxelNum.z;
	int yAll = voxelNum.y;
	int xAll = voxelNum.x;
	int windowVoxelNum = (windowHalf * 2 + 1)*(windowHalf * 2 + 1)*(windowHalf * 2 + 1);

	float *imgtmp = nullptr;//备份图像
	imgtmp = new float[voxelNumInAll];
	if (imgtmp == nullptr)
		return false;
	for (int i = 0; i < voxelNumInAll; i++)
		imgtmp[i] = pData[i];

#pragma omp parallel for
	for (int z = windowHalf; z < zAll - windowHalf; z++)
	{
		for (int y = windowHalf; y < yAll - windowHalf; y++)
		{
			for (int x = windowHalf; x < xAll - windowHalf; x++)
			{
				float temp = 0;
				for (int k = -1 * windowHalf; k <= windowHalf; k++)//求均值
				{
					for (int j = -1 * windowHalf; j <= windowHalf; j++)
					{
						for (int i = -1 * windowHalf; i <= windowHalf; i++)
						{
							temp += imgtmp[(x + i) + (y + j)*xAll + (z + k)*xAll*yAll];
						}
					}
				}
				temp /= windowVoxelNum;
				pData[x + y * xAll + z * xAll*yAll] = temp;
			}
		}
	}
	delete[]imgtmp;
	return true;
}
/**********************************************************
Description:		3D median smoothing
Arguments:
{
	windowHalf:		The half of the filtering window.
}
Return:				Status.
**********************************************************/
bool UImage::MedianDenoising3D(int windowHalf)const
{

	if (m_pData == nullptr)
	{
		return false;
	}
	if (windowHalf < 1)
	{
		return false;
	}
	float *pData = m_pData;
	int voxelNumInAll = m_nVoxelNumInAll;

	int zAll = voxelNum.z;
	int yAll = voxelNum.y;
	int xAll = voxelNum.x;
	int windowVoxelNum = (windowHalf * 2 + 1)*(windowHalf * 2 + 1)*(windowHalf * 2 + 1);

	float *imgtmp = nullptr;//备份图像
	imgtmp = new float[voxelNumInAll];
	if (imgtmp == nullptr)
		return false;
	for (int i = 0; i < voxelNumInAll; i++)
		imgtmp[i] = pData[i];

#pragma omp parallel for
	for (int z = windowHalf; z < zAll - windowHalf; z++)
	{
		float *voxelInWindow = new float[windowVoxelNum];
		for (int y = windowHalf; y < yAll - windowHalf; y++)
		{
			for (int x = windowHalf; x < xAll - windowHalf; x++)
			{
				int id = 0;
				for (int k = -1 * windowHalf; k <= windowHalf; k++)//取出窗口内的值
				{
					for (int j = -1 * windowHalf; j <= windowHalf; j++)
					{
						for (int i = -1 * windowHalf; i <= windowHalf; i++)
						{
							voxelInWindow[id] = imgtmp[(x + i) + (y + j)*xAll + (z + k)*xAll*yAll];
							id++;
						}
					}
				}
				float kk;//exchange buffer
				for (int i = 0; i < windowVoxelNum; i++)//sort
				{
					for (int j = i + 1; j < windowVoxelNum; j++)
					{
						if (voxelInWindow[i] > voxelInWindow[j])
						{
							kk = voxelInWindow[i];
							voxelInWindow[i] = voxelInWindow[j];
							voxelInWindow[j] = kk;
						}
					}
				}
				pData[x + y * xAll + z * xAll*yAll] = voxelInWindow[int((windowVoxelNum + 1) / 2)];
			}
		}
		delete[]voxelInWindow;
	}
	delete[]imgtmp;
	return true;
}
/**********************************************************
Description:		3D Guassian smoothing, and each direction is independent.
Arguments:
{
	windowHalf:		The half of the filtering window.
	FWHM:			Full width at half maximum of the Guassian kernel.
}
Return:				Status.
**********************************************************/
bool UImage::GuassianDenoising3D(int windowHalf, float FWHM)const
{
	if (m_pData == nullptr)
	{
		return false;
	}
	if (windowHalf < 1 || FWHM <= 0)
	{
		return false;
	}
	float sita = FWHM/2.355f;
	float *pData = m_pData;
	int voxelNumInAll = m_nVoxelNumInAll;

	int zAll = voxelNum.z;
	int yAll = voxelNum.y;
	int xAll = voxelNum.x;
	int windowVoxelNum = (windowHalf * 2 + 1)*(windowHalf * 2 + 1)*(windowHalf * 2 + 1);
	float *kernel = new float[windowVoxelNum];
	float kernelAll = 0;//高斯核分母
	float pi = 3.1415926535f;
	for (int k = 0; k < windowHalf * 2 + 1; k++)//计算高斯核
	{
		for (int j = 0; j < windowHalf * 2 + 1; j++)
		{
			for (int i = 0; i < windowHalf * 2 + 1; i++)
			{
				int id = i + j * (windowHalf * 2 + 1) + k * (windowHalf * 2 + 1)*(windowHalf * 2 + 1);
				float coff = static_cast<float>((i - windowHalf)*(i - windowHalf)*voxelSize.x*voxelSize.x + 
					(j - windowHalf)*(j - windowHalf)*voxelSize.y*voxelSize.y + (k - windowHalf)*(k - windowHalf)*voxelSize.z*voxelSize.z);
				coff = coff / (2 * sita*sita);
				coff = exp(-1 * coff) / (2 * pi * sita*sita);
				kernel[id] = coff;
				kernelAll += kernel[id];
			}
		}
	}
	for (int i = 0; i < windowVoxelNum; i++)//将高斯核归一化
		kernel[i] = kernel[i] / kernelAll;

	float *imgtmp = nullptr;//备份图像
	imgtmp = new float[voxelNumInAll];
	if (imgtmp == nullptr)
		return false;
	for (int i = 0; i < voxelNumInAll; i++)
		imgtmp[i] = pData[i];

#pragma omp parallel for
	for (int z = windowHalf; z < zAll - windowHalf; z++)
	{
		for (int y = windowHalf; y < yAll - windowHalf; y++)
		{
			for (int x = windowHalf; x < xAll - windowHalf; x++)
			{
				float filted = 0;//滤波后结果
				for (int k = -1 * windowHalf; k <= windowHalf; k++)//取出窗口内的值乘上高斯核
				{
					for (int j = -1 * windowHalf; j <= windowHalf; j++)
					{
						for (int i = -1 * windowHalf; i <= windowHalf; i++)
						{
							filted += imgtmp[(x + i) + (y + j)*xAll + (z + k)*xAll*yAll] *
								kernel[(i + windowHalf) + (j + windowHalf)*(2 * windowHalf + 1) + (k + windowHalf)*(2 * windowHalf + 1)*(2 * windowHalf + 1)];
						}
					}
				}
				pData[x + y * xAll + z * xAll*yAll] = filted;
			}
		}
	}
	delete[] imgtmp;
	delete[] kernel;
	return true;
}
/**********************************************************
Description:		3D filtering for the Negative pixel.
Arguments:
{
	windowHalf:		The half of the filtering window.
}
Return:				Status.
**********************************************************/
bool UImage::NegativeFilter3D(int windowHalf)const
{
	if (m_pData == nullptr)
	{
		return false;
	}
	if (windowHalf < 1)
	{
		return false;
	}
	float *pData = m_pData;
	int voxelNumInAll = m_nVoxelNumInAll;

	int zAll = voxelNum.z;
	int yAll = voxelNum.y;
	int xAll = voxelNum.x;
	int windowVoxelNum = (windowHalf * 2 + 1)*(windowHalf * 2 + 1)*(windowHalf * 2 + 1);

	float *imgtmp = nullptr;//备份图像
	imgtmp = new float[voxelNumInAll];
	if (imgtmp == nullptr)
		return false;
	for (int i = 0; i < voxelNumInAll; i++)
		imgtmp[i] = pData[i];

#pragma omp parallel for
	for (int z = windowHalf; z < zAll - windowHalf; z++)
	{
		int *kernel = new int[windowVoxelNum];//记录窗口内大于0的体素id
		for (int y = windowHalf; y < yAll - windowHalf; y++)
		{
			for (int x = windowHalf; x < xAll - windowHalf; x++)
			{
				int positiveNum = 0;//窗口内大于0的体素个数
				if (imgtmp[x + y * xAll + z * xAll*yAll] >= 0)//大于等于0的不管
					continue;
				for (int k = -1 * windowHalf; k <= windowHalf; k++)//取出窗口内大于等于0的体素id
				{
					for (int j = -1 * windowHalf; j <= windowHalf; j++)
					{
						for (int i = -1 * windowHalf; i <= windowHalf; i++)
						{
							if (imgtmp[(x + i) + (y + j) * xAll + (z + k) * xAll*yAll] > 0)
							{
								kernel[positiveNum] = (x + i) + (y + j) * xAll + (z + k) * xAll*yAll;
								positiveNum++;
							}
						}
					}
				}
				pData[x + y * xAll + z * xAll*yAll] = 0;
				for (int i = 0; i < positiveNum; i++)
					pData[x + y * xAll + z * xAll*yAll] += imgtmp[kernel[i]] / positiveNum;
			}
		}
		delete[] kernel;
	}

	for (int i = 0; i < voxelNumInAll; i++)
		if (pData[i] < 0)
			pData[i] = 0;
	delete[] imgtmp;
	return true;
}
/**********************************************************
Description:		Set the pixel out of the FOV to 0.
Arguments:
{

}
Return:				Status.
**********************************************************/
bool UImage::SetOutOfViewZero()const
{
	if (m_pData == nullptr)
	{
		return false;
	}
	int rr = (max(voxelNum.x, voxelNum.y) - 1) * (max(voxelNum.x, voxelNum.y) - 1) / 4;
	for (int zz = 0; zz < voxelNum.z; zz++)
	{
		for (int yy = 0; yy < voxelNum.y; yy++)
		{
			for (int xx = 0; xx < voxelNum.x; xx++)
			{
				if ((xx - voxelNum.x / 2)*(xx - voxelNum.x / 2) + ((yy - voxelNum.y / 2)*(yy - voxelNum.y / 2)) > rr)
				{
					m_pData[xx + yy * voxelNum.x + zz * voxelNum.y*voxelNum.x] = 0;
				}
			}
		}
	}
	return true;
}
