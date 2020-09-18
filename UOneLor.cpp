/********************************************************************

Copyright (C), 2019, All rights reserved

File Name     :    UOneLor.cpp
Description   :
History       :

<author>            <time>            <desc>
Ang Li             2019/6/8           create

********************************************************************/
#include "UOneLor.h"
#include <string.h>

UOneLor::UOneLor(int nVoxelNumMax) :
m_nVoxelNum(0), m_nVoxelNumMax(nVoxelNumMax)
{
	pVoxelID = new int[m_nVoxelNumMax];
	pWeight = new float[m_nVoxelNumMax];
	memset(pVoxelID, 0, m_nVoxelNumMax * sizeof(int));
	memset(pWeight, 0, m_nVoxelNumMax * sizeof(float));
}

UOneLor::~UOneLor()
{
	if (pVoxelID != nullptr)delete[] pVoxelID;
	pVoxelID = nullptr;
	if (pWeight != nullptr)delete[] pWeight;
	pWeight = nullptr;
}

UOneLor::UOneLor(const UOneLor& lor)
{
	m_nVoxelNumMax = lor.m_nVoxelNumMax;
	m_nVoxelNum = lor.m_nVoxelNum;
	pVoxelID = new int[m_nVoxelNumMax];
	pWeight = new float[m_nVoxelNumMax];
	memset(pVoxelID, 0, m_nVoxelNumMax * sizeof(int));
	memset(pWeight, 0, m_nVoxelNumMax * sizeof(float));
	for (int i = 0; i < m_nVoxelNum; i++)
	{
		pVoxelID[i] = lor.pVoxelID[i];
		pWeight[i] = lor.pWeight[i];
	}
}
/**********************************************************
Description:		Set the number of pixel on a LOR. It will automatically expand the memory when memory are not enough.
Arguments:
{
	nVoxelNum:		The number of voxel in LOR will be set.
}
Return:				Status.
**********************************************************/
bool UOneLor::SetVoxelNum(int nVoxelNum)
{
	if (nVoxelNum <= m_nVoxelNumMax && nVoxelNum >= 0)
	{
		m_nVoxelNum = nVoxelNum;
		return true;
	}
	else if (nVoxelNum > m_nVoxelNumMax)//数组最大长度扩展成设置大小的两倍
	{
		m_nVoxelNumMax = 2 * nVoxelNum;
		int *pVoxelIDTemp = nullptr;
		float *pWeightTemp = nullptr;
		pVoxelIDTemp = new int[m_nVoxelNumMax];
		pWeightTemp = new float[m_nVoxelNumMax];
		//如果new失败
		if (pVoxelIDTemp == nullptr || pWeightTemp == nullptr)
			return false;
		memset(pVoxelID, 0, m_nVoxelNumMax * sizeof(int));
		memset(pWeight, 0, m_nVoxelNumMax * sizeof(float));
		//复制原始的数组内容到新的数组
		memcpy(pVoxelIDTemp, pVoxelID, m_nVoxelNum * sizeof(int));
		memcpy(pWeightTemp, pWeight, m_nVoxelNum * sizeof(float));
		//释放原始数组
		delete[]pVoxelID;
		delete[]pWeight;
		//将指针指向新的数
		pVoxelID = pVoxelIDTemp;
		pWeight = pWeightTemp;
		//数组扩展完成
		m_nVoxelNum = nVoxelNum;
		return true;
	}
	else//nVoxelNum < 0
	{
		m_nVoxelNum = 0;
		return true;
	}
}
/**********************************************************
Description:		Add a voxel and its weight in this LOR.
Arguments:
{
	nVoxelID:		Voxel index
	fWeight:		System response weight of the voxel in this LOR
}
Return:				Status.
**********************************************************/
bool UOneLor::AddVoxel(int nVoxelID, float fWeight)
{
	if ((m_nVoxelNum+1) < m_nVoxelNumMax)
	{
		pVoxelID[m_nVoxelNum] = nVoxelID;
		pWeight[m_nVoxelNum] = fWeight;
		m_nVoxelNum++;//已存的像素数+1
		return true;
	}
	else//(m_nVoxelNum+1) >= m_nVoxelNumMax;如果内存已满，则扩大2倍内存容量
	{
		m_nVoxelNumMax *= 2;
		int *pVoxelIDTemp = nullptr;
		float *pWeightTemp = nullptr;
		pVoxelIDTemp = new int[m_nVoxelNumMax];
		pWeightTemp = new float[m_nVoxelNumMax];
		//如果new失败
		if (pVoxelIDTemp == nullptr || pWeightTemp == nullptr)
			return false;
		memset(pVoxelID, 0, m_nVoxelNumMax * sizeof(int));
		memset(pWeight, 0, m_nVoxelNumMax * sizeof(float));
		//复制原始的数组内容到新的数组
		memcpy(pVoxelIDTemp, pVoxelID, m_nVoxelNum * sizeof(int));
		memcpy(pWeightTemp, pWeight, m_nVoxelNum * sizeof(float));
		//释放原始数组
		delete[]pVoxelID;
		delete[]pWeight;
		//将指针指向新的数
		pVoxelID = pVoxelIDTemp;
		pWeight = pWeightTemp;
		//数组扩展完成，执行add操作
		pVoxelID[m_nVoxelNum] = nVoxelID;
		pWeight[m_nVoxelNum] = fWeight;
		m_nVoxelNum++;//已存的像素数+1
		return true;
	}
}
/**********************************************************
Description:		Clear all the voxels and weights in this LOR.
Arguments:
{

}
Return:				void.
**********************************************************/
void UOneLor::Clear()
{
	memset(pVoxelID, 0, m_nVoxelNum*sizeof(int));
	memset(pWeight, 0, m_nVoxelNum*sizeof(float));
	m_nVoxelNum = 0;
}
/**********************************************************
Description:		Get the number of voxels in this LOR.
Arguments:
{

}
Return:				The number of voxels in this LOR.
**********************************************************/
int UOneLor::GetVoxelNum() const
{
	return m_nVoxelNum;
}