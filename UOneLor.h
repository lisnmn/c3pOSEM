/********************************************************************

Copyright (C), 2019, All rights reserved

File Name     :    UOneLor.h
Description   :
History       :

<author>            <time>            <desc>
Ang Li             2019/6/8           create

********************************************************************/
#ifndef __ONE_LOR_H
#define __ONE_LOR_H

/*********************************************************************************
* 记录一条光线穿过图像各体素的索引和长度
* 20190125优化：实例化时初始化数组最大长度，使用时若最大长度不够，自动扩展
**********************************************************************************/
class UOneLor
{
public:
	UOneLor(int nVoxelNumMax);		//初始化赋予一个初始最大数组长度
	~UOneLor();
	UOneLor(const UOneLor& lor);		//拷贝构造函数

public:
	/**********************************************************
	Description:		Set the number of pixel on a LOR. It will automatically expand the memory when memory are not enough.
	Arguments:
	{
	nVoxelNum:		The number of voxel in LOR will be set.
	}
	Return:				Status.
	**********************************************************/
	bool SetVoxelNum(int nVoxelNum);
	/**********************************************************
	Description:		Add a voxel and its weight in this LOR.
	Arguments:
	{
	nVoxelID:		Voxel index
	fWeight:		System response weight of the voxel in this LOR
	}
	Return:				Status.
	**********************************************************/
	bool AddVoxel(int nVoxelID, float weight);
	/**********************************************************
	Description:		Clear all the voxels and weights in this LOR.
	Arguments:
	{

	}
	Return:				void.
	**********************************************************/
	void Clear();
	/**********************************************************
	Description:		Get the number of voxels in this LOR.
	Arguments:
	{

	}
	Return:				The number of voxels in this LOR.
	**********************************************************/
	int GetVoxelNum()const;

public:
	//为了方便使用public变量不加[m_]前缀
	float *pWeight;					//表示系统响应矩阵的概率值，在raytracing中表示光线穿过体素的长度，单位设为mm减少计算量
	int *pVoxelID;
private:
	UOneLor& operator=(const UOneLor& lor){};
private:
	int m_nVoxelNum;
	int m_nVoxelNumMax;
};

#endif
