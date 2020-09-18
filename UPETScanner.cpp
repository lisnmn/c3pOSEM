/********************************************************************

Copyright (C), 2019, All rights reserved

File Name     :    UPETScanner.cpp
Description   :
History       :

<author>            <time>            <desc>
Ang Li             2019/6/8           create
405                2020/3/27          modify

********************************************************************/
#include "UPETScanner.h"
#include "UVolume.h"
#include <stdlib.h>
#include <math.h>



#define	PI			3.141592654


void UPETScanner::InitD80Scanner()
{
	crystalNum.set(1, 13, 13);//这里同Gate坐标系
	crystalSize.set(13, 1.89, 1.89);
	crystalPitch.set(13, 2.0, 2.0);

	blockNum.set(1, 1, 4);
	blockSize.set(13.3, crystalPitch.y * crystalNum.y, crystalPitch.z * crystalNum.z);
	blockPitch.set(13.3, 26.5, 26.5);//不忽略block edge，单位mm

	moduleNum.set(1, 1, 1);
	moduleSize.set(20, blockPitch.y * blockNum.y, blockPitch.z * blockNum.z);
	modulePitch.set(20, moduleSize.y, moduleSize.z + 2);

	panelNum = 12;
	panelSize.set(20, modulePitch.y * moduleNum.y, modulePitch.z * moduleNum.z);
	panelPitch = panelSize;

	crystalClockwiseOffset = 1;
	scannerRadius = 105.8 / 2;
}

void UPETScanner::InitD180Scanner()
{
	crystalNum.x = 1;
	crystalNum.y = 13;//这里同Gate坐标系
	crystalNum.z = 13;
	crystalSize.x = 13;
	crystalSize.y = 1.89;
	crystalSize.z = 1.89;
	crystalPitch.x = 13;
	crystalPitch.y = 2.0;
	crystalPitch.z = 2.0;

	blockNum.x = 1;
	blockNum.y = 1;
	blockNum.z = 4;
	blockSize.x = 13.3;
	blockSize.y = crystalPitch.y * crystalNum.y;
	blockSize.z = crystalPitch.z * crystalNum.z;
	blockPitch.x = 13.3;
	blockPitch.y = 26.5;
	blockPitch.z = 26.5;//不忽略block edge，单位mm

	moduleNum.x = 1;
	moduleNum.y = 1;
	moduleNum.z = 1;
	moduleSize.x = 20;
	moduleSize.y = blockPitch.y * blockNum.y;
	moduleSize.z = blockPitch.z * blockNum.z;
	modulePitch.x = 20;
	modulePitch.y = moduleSize.y;
	modulePitch.z = moduleSize.z + 2;

	panelNum = 24;
	panelSize.x = 20;
	panelSize.y = modulePitch.y * moduleNum.y;
	panelSize.z = modulePitch.z * moduleNum.z;
	panelPitch.x = 20;
	panelPitch.y = panelSize.y;
	panelPitch.z = panelSize.z;

	crystalClockwiseOffset = 1;
	scannerRadius = 106.5;//real is 106.5

}

void UPETScanner::InitE180Scanner()
{
	crystalNum.x = 1;
	crystalNum.y = 13;//这里同Gate坐标系
	crystalNum.z = 13;
	crystalSize.x = 13;
	crystalSize.y = 1.89;
	crystalSize.z = 1.89;
	crystalPitch.x = 13;
	crystalPitch.y = 2.0;
	crystalPitch.z = 2.0;

	blockNum.x = 1;
	blockNum.y = 1;
	blockNum.z = 4;
	blockSize.x = 13.3;
	blockSize.y = crystalPitch.y * crystalNum.y;
	blockSize.z = crystalPitch.z * crystalNum.z;
	blockPitch.x = 13.3;
	blockPitch.y = 26.5;
	blockPitch.z = 26.5;//不忽略block edge，单位mm

	moduleNum.x = 1;
	moduleNum.y = 1;
	moduleNum.z = 2;
	moduleSize.x = 20;
	moduleSize.y = blockPitch.y * blockNum.y;
	moduleSize.z = blockPitch.z * blockNum.z;
	modulePitch.x = 20;
	modulePitch.y = moduleSize.y;
	modulePitch.z = moduleSize.z + 2;

	panelNum = 24;
	panelSize.x = 20;
	panelSize.y = modulePitch.y * moduleNum.y;
	panelSize.z = modulePitch.z * moduleNum.z;
	panelPitch.x = 20;
	panelPitch.y = panelSize.y;
	panelPitch.z = panelSize.z;

	crystalClockwiseOffset = 1;
	scannerRadius = 106.5;//real is 106.5
}

void UPETScanner::InitTest()
{
	crystalNum.set(1, 13, 13);//这里同Gate坐标系
	crystalSize.set(13, 1.89, 1.89);
	crystalPitch.set(13, 2.0, 2.0);

	blockNum.set(1, 1, 2);
	blockSize.set(13.3, crystalPitch.y * crystalNum.y, crystalPitch.z * crystalNum.z);
	blockPitch.set(13.3, 26.5, 26.5);//不忽略block edge，单位mm

	moduleNum.set(1, 1, 1);
	moduleSize.set(20, blockPitch.y * blockNum.y, blockPitch.z * blockNum.z);
	modulePitch.set(20, moduleSize.y, moduleSize.z + 2);

	panelNum = 12;
	panelSize.set(20, modulePitch.y * moduleNum.y, modulePitch.z * moduleNum.z);
	panelPitch = panelSize;

	crystalClockwiseOffset = 1;
	scannerRadius = 105.8 / 2;
}

int UPETScanner::GetCrystalNumOneRing() const
{
	return panelNum * moduleNum.y * blockNum.y * crystalNum.y;
}

int UPETScanner::GetRingNum() const
{
	return moduleNum.z * blockNum.z * crystalNum.z;
}

int UPETScanner::GetBinNum()const
{
	return GetCrystalNumOneRing() - 1;
}

int UPETScanner::GetViewNum()const
{
	return GetCrystalNumOneRing() / 2;
}

int UPETScanner::GetSliceNum()const
{
	return GetRingNum()*GetRingNum();
}

int UPETScanner::GetLORNum()const
{
	return GetBinNum()*GetViewNum()*GetSliceNum();
}

int UPETScanner::GetCrystalNum()const
{
	return GetCrystalNumOneRing()*GetRingNum();
}

int UPETScanner::GetCrystalNumZInPanel() const
{
	return GetRingNum();
}

int UPETScanner::GetCrystalNumYInPanel() const
{
	return crystalNum.y*blockNum.y*moduleNum.y;
}

int UPETScanner::GetCrystalNumZInModule() const
{
	return crystalNum.z*blockNum.z;
}

int UPETScanner::GetCrystalNumYInModule() const
{
	return crystalNum.y*blockNum.y;
}

int UPETScanner::GetCrystalNumZInBlock() const
{
	return crystalNum.z;
}

int UPETScanner::GetCrystalNumYInBlock() const
{
	return crystalNum.y;
}

int UPETScanner::GetPanelNum() const
{
	return panelNum;
}

double UPETScanner::GetLengthZ() const
{
	return modulePitch.z*moduleNum.z;
}

double UPETScanner::GetPanelSizeZ() const
{
	return panelSize.z;
}

double UPETScanner::GetPanelSizeY() const
{
	return panelSize.y;
}

double UPETScanner::GetModulePitchZ() const
{
	return modulePitch.z;
}

double UPETScanner::GetModulePitchY() const
{
	return modulePitch.y;
}

double UPETScanner::GetModuleSizeZ()const
{
	return moduleSize.z;
}

double UPETScanner::GetModuleSizeY() const
{
	return moduleSize.y;
}

double UPETScanner::GetBlockPitchZ() const
{
	return blockPitch.z;
}

double UPETScanner::GetBlockPitchY() const
{
	return blockPitch.y;
}

double UPETScanner::GetBlockSizeZ() const
{
	return blockSize.z;
}

double UPETScanner::GetBlockSizeY() const
{
	return blockSize.y;
}

double UPETScanner::GetCrystalPitchY() const
{
	return crystalPitch.y;
}

double UPETScanner::GetCrystalPitchZ() const
{
	return crystalPitch.z;
}

double UPETScanner::GetCrystalSizeY() const
{
	return crystalSize.y;
}

double UPETScanner::GetCrystalSizeZ() const
{
	return crystalSize.z;
}

double UPETScanner::GetRadius() const
{
	return scannerRadius;
}

int UPETScanner::GetCrystalClockwiseOffset() const
{
	return crystalClockwiseOffset;
}
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
void UPETScanner::GetCrystalPosition(int crystalInRing, int ring, UVolume3D<double> *crystalPos)const
{
	//加入晶体编号的偏置
	crystalInRing = (crystalInRing - crystalClockwiseOffset + GetCrystalNumOneRing()) % GetCrystalNumOneRing();
	//先计算的z坐标
	double z = GetRingMinCoordinateZ(ring);

	//计算晶体四个点的z坐标值
	crystalPos[0].z = z;//将坐标原点移至中心
	crystalPos[1].z = crystalPos[0].z;
	crystalPos[2].z = z + crystalSize.z;
	crystalPos[3].z = crystalPos[2].z;

	//下面计算x和y的坐标值
	int panel1, module1, block1, cry1;
	panel1 = crystalInRing / GetCrystalNumYInPanel();//计算晶体在环的哪一个panel上
	module1 = crystalInRing % GetCrystalNumYInPanel() / GetCrystalNumYInModule();//计算晶体在(panel)的哪一个module上
	block1 = crystalInRing % GetCrystalNumYInModule() / GetCrystalNumYInBlock();//计算晶体在module的哪一个个block上
	cry1 = crystalInRing % GetCrystalNumYInBlock();//计算晶体是block中的第几个晶体
												   //在panel上的偏移值
	double sita1;//panel中心所对的角度
	sita1 = double(panel1) * 2 * PI / panelNum;

	double panelOffset1;//环内晶体的一边偏离所在block中心的距离
	panelOffset1 = module1 * modulePitch.y + (modulePitch.y - moduleSize.y) / 2 +
		block1 * blockPitch.y + (blockPitch.y - blockSize.y) / 2 +
		cry1 * crystalPitch.y + (crystalPitch.y - crystalSize.y) / 2 -
		GetPanelSizeY() / 2;
	//得到xy的坐标
	crystalPos[0].x = scannerRadius * cos(sita1) + panelOffset1 * cos(sita1 + PI / 2);
	crystalPos[0].y = scannerRadius * sin(sita1) + panelOffset1 * sin(sita1 + PI / 2);
	crystalPos[3].x = crystalPos[0].x;
	crystalPos[3].y = crystalPos[0].y;

	crystalPos[1].x = scannerRadius * cos(sita1) + (panelOffset1 + crystalSize.y)*cos(sita1 + PI / 2);
	crystalPos[1].y = scannerRadius * sin(sita1) + (panelOffset1 + crystalSize.y)*sin(sita1 + PI / 2);
	crystalPos[2].x = crystalPos[1].x;
	crystalPos[2].y = crystalPos[1].y;
}
/**********************************************************
Description:		Get the minimum z coordinate of one ring.
Arguments:
{
	ring:			Ring index of this crystal. At the range of [0,ringNum-1]
}
Return:				Status.
**********************************************************/
double UPETScanner::GetRingMinCoordinateZ(int ring) const
{
	double z = int(ring / GetCrystalNumZInModule()) * modulePitch.z + (modulePitch.z - moduleSize.z) / 2 +					//module从中心开始，两边各一半的pitch-size
		int(ring%GetCrystalNumZInModule() / GetCrystalNumZInBlock())*blockPitch.z + (blockPitch.z - blockSize.z) / 2 +		//block和Crystal从中心开始，两边各一半的pitch-size
		ring % GetCrystalNumZInBlock() * crystalPitch.z + (crystalPitch.z - crystalSize.z) / 2 -
		GetLengthZ() / 2;					//将坐标移至中心
	return z;
}
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
void UPETScanner::GetCrystalIDFromLORID(int LORID, int &crystalID1, int &crystalID2)const
{
	int crystalNumOneRing = GetCrystalNumOneRing();
	int binNum = GetBinNum();
	int viewNum = GetViewNum();
	int bin = LORID%binNum;
	int view = LORID / binNum%viewNum;
	int slice = LORID / (binNum*viewNum);
	int cry1, cry2, ring1, ring2;
	GetCrystalIDInRingFromViewBin(view, bin, cry1, cry2);
	GetRing1Ring2FromSlice(slice, ring1, ring2);
	crystalID1 = cry1 + ring1*crystalNumOneRing;
	crystalID2 = cry2 + ring2*crystalNumOneRing;
}
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
void UPETScanner::GetLORIDFromCrystalID(int crystalID1, int crystalID2, int &LORID)const
{
	int crystalNumOneRing = GetCrystalNumOneRing();
	int cry1 = crystalID1%crystalNumOneRing;
	int cry2 = crystalID2%crystalNumOneRing;
	int ring1 = crystalID1 / crystalNumOneRing;
	int ring2 = crystalID2 / crystalNumOneRing;
	GetLORIDFromRingAndCrystalInRing(ring1, cry1, ring2, cry2, LORID);
}
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
void UPETScanner::GetRing1Ring2FromSlice(int slice, int & ring1, int & ring2)const
{
	ring1 = slice / GetRingNum();
	ring2 = slice % GetRingNum();
}
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
void UPETScanner::GetCrystalIDInRingFromViewBin(int view, int bin, int & cry1, int & cry2)const
{
	int crystalOneRing = GetCrystalNumOneRing();
	cry2 = bin / 2 + 1;
	cry1 = crystalOneRing + (1 - bin % 2) - cry2;

	cry2 = (cry2 + view) % crystalOneRing;
	cry1 = (cry1 + view) % crystalOneRing;
}
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
void UPETScanner::GetLORIDFromRingAndCrystalInRing(int ring1, int cry1, int ring2, int cry2, int & LORID)const
{
	//此处输入的ring1和ring2不能确定哪个是构成LOR的ring1和ring2，需要先根据cry1和cry2来判断。
	//因此需要先计算view和bin

	int cry1t = cry1;
	int cry2t = cry2;
	int crystalOneRing = GetCrystalNumOneRing();
	int view = (cry1 + cry2) % crystalOneRing / 2;
	//将cry1和cry2都还原到view为0的情况
	cry1 -= view;
	cry2 -= view;
	if (cry1 <= 0)//等于0也需要
		cry1 += crystalOneRing;			//20180524日改crystalOneRing / 2;
	if (cry2 <= 0)
		cry2 += crystalOneRing;			//20180524日改crystalOneRing / 2;
	int ring1real, ring2real, cry1real, cry2real;//对应LOR的ring1和ring2
	if (cry1 > cry2)
	{
		cry1real = cry1;
		cry2real = cry2;
		ring1real = ring1;
		ring2real = ring2;
	}
	else
	{
		cry1real = cry2;
		cry2real = cry1;
		ring1real = ring2;
		ring2real = ring1;
	}
	int bin = (crystalOneRing - 1) - (cry1real - cry2real);
	int slice = ring1real * GetRingNum() + ring2real;
	LORID = slice * (crystalOneRing - 1)*(crystalOneRing / 2) + view * (crystalOneRing - 1) + bin;
}
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
int UPETScanner::IsGoodPair(int panel1, int panel2, int minSectorDifference)const
{
	return (abs(panel1 - panel2) > minSectorDifference && panelNum - abs(panel1 - panel2) > minSectorDifference) ? 1 : 0;
}
