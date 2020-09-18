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
	else if (nVoxelNum > m_nVoxelNumMax)//������󳤶���չ�����ô�С������
	{
		m_nVoxelNumMax = 2 * nVoxelNum;
		int *pVoxelIDTemp = nullptr;
		float *pWeightTemp = nullptr;
		pVoxelIDTemp = new int[m_nVoxelNumMax];
		pWeightTemp = new float[m_nVoxelNumMax];
		//���newʧ��
		if (pVoxelIDTemp == nullptr || pWeightTemp == nullptr)
			return false;
		memset(pVoxelID, 0, m_nVoxelNumMax * sizeof(int));
		memset(pWeight, 0, m_nVoxelNumMax * sizeof(float));
		//����ԭʼ���������ݵ��µ�����
		memcpy(pVoxelIDTemp, pVoxelID, m_nVoxelNum * sizeof(int));
		memcpy(pWeightTemp, pWeight, m_nVoxelNum * sizeof(float));
		//�ͷ�ԭʼ����
		delete[]pVoxelID;
		delete[]pWeight;
		//��ָ��ָ���µ���
		pVoxelID = pVoxelIDTemp;
		pWeight = pWeightTemp;
		//������չ���
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
		m_nVoxelNum++;//�Ѵ��������+1
		return true;
	}
	else//(m_nVoxelNum+1) >= m_nVoxelNumMax;����ڴ�������������2���ڴ�����
	{
		m_nVoxelNumMax *= 2;
		int *pVoxelIDTemp = nullptr;
		float *pWeightTemp = nullptr;
		pVoxelIDTemp = new int[m_nVoxelNumMax];
		pWeightTemp = new float[m_nVoxelNumMax];
		//���newʧ��
		if (pVoxelIDTemp == nullptr || pWeightTemp == nullptr)
			return false;
		memset(pVoxelID, 0, m_nVoxelNumMax * sizeof(int));
		memset(pWeight, 0, m_nVoxelNumMax * sizeof(float));
		//����ԭʼ���������ݵ��µ�����
		memcpy(pVoxelIDTemp, pVoxelID, m_nVoxelNum * sizeof(int));
		memcpy(pWeightTemp, pWeight, m_nVoxelNum * sizeof(float));
		//�ͷ�ԭʼ����
		delete[]pVoxelID;
		delete[]pWeight;
		//��ָ��ָ���µ���
		pVoxelID = pVoxelIDTemp;
		pWeight = pWeightTemp;
		//������չ��ɣ�ִ��add����
		pVoxelID[m_nVoxelNum] = nVoxelID;
		pWeight[m_nVoxelNum] = fWeight;
		m_nVoxelNum++;//�Ѵ��������+1
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