/********************************************************************

Copyright (C), 2019, All rights reserved

File Name     :    UPETScanner.h
Description   :    
History       :    

<author>            <time>            <desc>
Ang Li             2019/6/8           create
405                2020/3/27          modify

********************************************************************/

/************************************************************************************
* PETScanner：定义PET探测器结构的类
* 使用方法：1.调用预置初始化函数；2.继承该类并重写Init()函数初始化探测器结构信息
* 预留配置函数：InitD80Scanner()InitD180Scanner()InitE80Scanner()
* 注：UVolume3D<double>是模板类，需要设置xyz坐标的信息。方向同GATE下的坐标系。
**************************************************************************************/

#ifndef __PET_SCANNER_H
#define __PET_SCANNER_H

#include "UVolume.h"
#include <iostream>
#include <stdio.h>
#include <string.h>

class UPETScanner
{
public:
	UPETScanner(void) {};//inline
	~UPETScanner(void) {};
	//set&get
	int GetCrystalNumOneRing()const;
	int GetRingNum()const;
	int GetBinNum()const;//20190612Updated
	int GetViewNum()const;
	int GetSliceNum()const;
	int GetLORNum()const;
	int GetCrystalNum()const;//20190612Updated
	int GetCrystalNumZInPanel()const;
	int GetCrystalNumYInPanel()const;
	int GetCrystalNumZInModule()const;
	int GetCrystalNumYInModule()const;
	int GetCrystalNumZInBlock()const;
	int GetCrystalNumYInBlock()const;
	int GetPanelNum()const;
	double GetLengthZ()const;
	double GetPanelSizeZ()const;
	double GetPanelSizeY()const;
	double GetModulePitchZ()const;
	double GetModulePitchY()const;
	double GetModuleSizeZ()const;//20190616Updated
	double GetModuleSizeY()const;
	double GetBlockPitchZ()const;
	double GetBlockPitchY()const;
	double GetBlockSizeZ()const;
	double GetBlockSizeY()const;
	double GetCrystalPitchY()const;
	double GetCrystalPitchZ()const;
	double GetCrystalSizeY()const;
	double GetCrystalSizeZ()const;
	double GetRadius()const;
	int GetCrystalClockwiseOffset()const;
	//根据晶体所在的环和晶体在环上的编号0-GetCrystalNumOneRing()-1求晶体坐标，返回四个点的顺序为0,1,2,3
	/**********************************************************
	Description:		Get the 4 vertices positions of one crystal.
	Arguments:
	{
	crystalInRing:	Crystal index in ring. At the range of [0,crystalNumOneRing-1]
	ring:			Ring index of this crystal. At the range of [0,ringNum-1]
	crystalPos:		get the crystal position which has 4 element: UVolume3D<double> crystalPos[4];
	}
	Return:				Status.
	**********************************************************/
	void GetCrystalPosition(int crystalInRing, int ring, UVolume3D<double> *crystal)const;
	//获得环的最小z轴坐标
	/**********************************************************
	Description:		Get the minimum z coordinate of one ring.
	Arguments:
	{
	ring:			Ring index of this crystal. At the range of [0,ringNum-1]
	}
	Return:				Status.
	**********************************************************/
	double GetRingMinCoordinateZ(int ring)const;
	//根据LOR编号得到LOR上两个晶体的全局晶体编号
	/**********************************************************
	Description:		LORID to global crystal ID.
	Arguments:
	{
	LORID:			LOR index. At the range of [0,LORNum-1]
	crystalID1:		Get the global crystal1 index. At the range of [0,crystalNum-1]
	crystalID2:		Get the global crystal2 index. At the range of [0,crystalNum-1]
	}
	Return:				void.
	**********************************************************/
	void GetCrystalIDFromLORID(int LORID, int &crystalID1, int &crystalID2)const;
	//根据LOR上两个晶体的全局编号得到LOR编号
	/**********************************************************
	Description:		global crystal ID to LOR ID.
	Arguments:
	{
	crystalID1:		Get the global crystal1 index. At the range of [0,crystalNum-1]
	crystalID2:		Get the global crystal2 index. At the range of [0,crystalNum-1]
	LORID:			Get the LOR index. At the range of [0,LORNum-1]
	}
	Return:				void.
	**********************************************************/
	void GetLORIDFromCrystalID(int crystalID1, int crystalID2, int &LORID)const;
	//根据slice得到ring1和ring2
	/**********************************************************
	Description:		LORID to crystal ID. Slice to ring1 and ring2
	Arguments:
	{
	slice:			slice index of this LOR.
	ring1:			Get the ring1 index. At the range of [0,ringNum-1]
	ring2:			Get the ring2 index. At the range of [0,ringNum-1]
	}
	Return:				Status.
	**********************************************************/
	void GetRing1Ring2FromSlice(int slice, int & ring1, int & ring2)const;
	//根据view和bin获得在ring1上的晶体cry1，ring2上的晶体cry2
	/**********************************************************
	Description:		LORID to crystal ID. view,bin to crystal1,crystal2
	Arguments:
	{
	view:			view index of this LOR.
	bin:			bin index of this LOR
	cry1:			Get the index of crystal1 on the ring. At the range of [0,crystalOneRing-1]
	cry2:			Get the index of crystal2 on the ring. At the range of [0,crystalOneRing-1]
	}
	Return:				Status.
	**********************************************************/
	void GetCrystalIDInRingFromViewBin(int view, int bin, int & cry1, int & cry2)const;
	//根据两晶体所在的环ring和晶体在一个环上的晶体编cry得到LORID
	/**********************************************************
	Description:		Crystal ID to LOR ID.
	Arguments:
	{
	ring1:			The ring1 index. At the range of [0,ringNum-1]
	cry1:			The index of crystal1 on the ring. At the range of [0,crystalOneRing-1]
	ring2:			The ring2 index. At the range of [0,ringNum-1]
	cry2:			The index of crystal2 on the ring. At the range of [0,crystalOneRing-1]
	LORID:			Get the LOR index.
	}
	Return:				Status.
	**********************************************************/
	void GetLORIDFromRingAndCrystalInRing(int ring1, int cry1, int ring2, int cry2, int & LORID)const;
	//判断是否为优质符合20180723
	/**********************************************************
	Description:		Determine if the 2 panel has a large distance.
	Arguments:
	{
	panel1:			The panel1 index.
	panel2:			The panel2 index.
	minSectorDifference:	The minimum difference between the 2 panel.
	}
	Return:				1 when get true. 0 when get false.
	**********************************************************/
	int IsGoodPair(int panel1, int panel2, int minSectorDifference)const;

public://提供的探测器结构设置接口
	void InitD80Scanner();
	void InitD180Scanner();
	void InitE180Scanner();
	void InitTest();
	virtual void Init(){};

	UVolume3D<double>	crystalSize; //晶体条真实尺寸(mm)
	UVolume3D<double>	crystalPitch;//晶体条包装尺寸(mm)
	UVolume3D<int>   	crystalNum;  //一个block中所含有的晶体条个数

	UVolume3D<double>	blockSize;   //block真实尺寸(mm)，可以从晶体条相关信息计算出来
	UVolume3D<double>	blockPitch;  //block包装尺寸(mm)
	UVolume3D<int>		blockNum;    //一个module中所含有的block个数

	UVolume3D<double>	moduleSize;  //module真实尺寸(mm)，可以从block相关信息计算出来
	UVolume3D<double>	modulePitch; //module包装尺寸(mm)
	UVolume3D<int>   	moduleNum;   //一个panel中所含有的module个数

	UVolume3D<double>	panelSize;   //panel真实尺寸(mm)，可以从module相关信息计算出来
	UVolume3D<double>	panelPitch;  //panel包装尺寸(mm)
	int	   		panelNum;    //探测环中所含有的panel个数
	int			crystalClockwiseOffset;		//晶体顺时针方向偏移量
	double 		scannerRadius; //必须要大于等于最小内切圆半径

private:
	UPETScanner(const UPETScanner&){};
	UPETScanner& operator=(const UPETScanner&){};
};

#endif
