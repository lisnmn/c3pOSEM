/********************************************************************

Copyright (C), 2019, All rights reserved

File Name     :    UReconPETOSEMAlgor.cpp
Description   :
History       :

<author>            <time>            <desc>
Ang Li             2019/6/8           create

********************************************************************/
#include "UReconPETOSEMAlgor.h"
#include "UProgressNotify.h"
#include "UFileReadSave.h"
#include "UReconPETFunction.h"
#include "UPETScanner.h"
#include "UImage.h"
#include "URayTracing.h"
#include "UOneLor.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

//#define max(x,y)	((x)<(y)?(y):(x))
//#define min(x,y)	((x)<(y)?(x):(y))
template<class T>
inline T max(T x, T y){ return ((x) < (y) ? (y) : (x)); }
template<class T>
inline T min(T x, T y){ return ((x) < (y) ? (x) : (y)); }

UReconPETOSEMPara::UReconPETOSEMPara()
{
	m_pScan = nullptr;
	m_nImgWidth = 160;//image pixel num
	m_nImgHeight = 160;
	m_nImgDepth = 100;
	m_fVoxelSizeXY = 1;//image pixel size
	m_fVoxelSizeZ = 1;

	m_nIterNum = 2;
	m_nSubsetNum = 12;
	m_nRingDifference = 100000;//a large number means slice with any ring difference will be used

	m_fTimeBeforeScan=1;
	m_fTimeDuringScan=1;
	m_fSourceHalfLife = 108*60.0f;
	m_fSourceBrachingRatio = 0.967f;
	m_fWellCounterSlope = 1;// well counter factor
	m_fWellCounterIntercept = 0;

	m_strMichPath = "./mich.dat";//mich file path
	m_strSenMapPath = "./senMap.dat";//sensitivity map path

	m_strNormFactorPath = "./norm.dat";//normalization file path
	m_strAttnFactorPath = "./attn.dat";//attenuation file path
	m_strRandEstimatePath = "./rand.dat";//random estimate file path
	m_strSavePath = "./img.dat";

	m_bCalcSenMap = true;
	m_bRandomCorr = false;
	m_bScatterCorr = false;
	m_bNormCorr = false;
	m_bAttnCorr = false;
	m_bWellCounterCorr = false;
	m_bDecayCorr = false;
	m_DeadTimeCorr = false;
}

UReconPETOSEMPara::~UReconPETOSEMPara()
{
}

//=========================================
//			UReconPETOSEMAlgor
//=========================================
UReconPETOSEMRayTracingAlgor::UReconPETOSEMRayTracingAlgor() :
	m_pPara(nullptr), m_pNotify(nullptr)
{
}

UReconPETOSEMRayTracingAlgor::~UReconPETOSEMRayTracingAlgor()
{
}
/**********************************************************
Description:		Initialize the PET reconstruction parameter.
Arguments:
{
	pPara:			Reconstraction parameter.
	pNotify:		Object to record the progress.
}
Return:				Status.
**********************************************************/
bool UReconPETOSEMRayTracingAlgor::Initial(UReconPETOSEMPara * pPara, UProgressNotify * pNotify)
{
	if (pPara == nullptr || pNotify == nullptr)
		return false;
	if (pPara->m_pScan == nullptr)
		return false;
	m_pPara = pPara;
	m_pNotify = pNotify;
	return true;
}
/**********************************************************
Description:		Reconstruct and save the PET image.
Arguments:
{

}
Return:				Status.
**********************************************************/
bool UReconPETOSEMRayTracingAlgor::DoRecon()
{
	//==================================
	//get information
	//==================================
	const UPETScanner *pScan = m_pPara->m_pScan;
	const int iterNum = m_pPara->m_nIterNum;
	const int subsetNum = m_pPara->m_nSubsetNum;
	const int crystalNumOneRing = pScan->GetCrystalNumOneRing();
	const int ringNum = pScan->GetRingNum();
	const int binNum = pScan->GetCrystalNumOneRing() - 1;
	const int viewNum = pScan->GetCrystalNumOneRing() / 2;
	const int sliceNum = pScan->GetRingNum()*pScan->GetRingNum();
	const int LORNum = binNum * viewNum * sliceNum;
	const int nImgVoxelNum = m_pPara->m_nImgWidth * m_pPara->m_nImgHeight * m_pPara->m_nImgDepth;
	const int threadNum = max(omp_get_num_procs() - 2, 1);
	UImage img;//image
	img.SetCenter(0, 0, 0);
	img.SetVoxelSize(m_pPara->m_fVoxelSizeXY, m_pPara->m_fVoxelSizeXY, m_pPara->m_fVoxelSizeZ);
	img.SetVoxelNum(m_pPara->m_nImgWidth, m_pPara->m_nImgHeight, m_pPara->m_nImgDepth);
	for (int i = 0; i < nImgVoxelNum; i++)
	{
		img.m_pData[i] = 1;
	}
	//==================================
	//alloc
	//==================================
	float *nimageMultiThreads = nullptr, *senMap = nullptr, *mich = nullptr, *michNorm = nullptr, *michAttn = nullptr, *michRand = nullptr;//update image and mich
	nimageMultiThreads = new float[nImgVoxelNum*threadNum];
	senMap = new float[nImgVoxelNum*subsetNum];
	mich = new float[LORNum];
	michNorm = new float[LORNum];
	michAttn = new float[LORNum];
	michRand = new float[LORNum];
	if (nimageMultiThreads == nullptr || senMap == nullptr || mich == nullptr || michNorm == nullptr || michAttn == nullptr || michRand == nullptr)
		return false;
	//==================================
	//initial array
	//==================================
	memset(nimageMultiThreads, 0, nImgVoxelNum*threadNum * sizeof(float));
	memset(mich, 0, LORNum * sizeof(float));
	memset(michRand, 0, LORNum * sizeof(float));
	for (int i = 0; i < LORNum; i++)
	{
		michAttn[i] = 1;
		michNorm[i] = 1;
	}
	//==================================
	//read mich
	//==================================
	if (!UFileReadSave::ReadFile((char*)mich, m_pPara->m_strMichPath, LORNum * sizeof(float)))
		return false;
	if (m_pPara->m_bAttnCorr)
	{
		if (!UFileReadSave::ReadFile((char*)michAttn, m_pPara->m_strAttnFactorPath, LORNum * sizeof(float)))
			return false;
	}
	if (m_pPara->m_bNormCorr)
	{
		if (!UFileReadSave::ReadFile((char*)michNorm, m_pPara->m_strNormFactorPath, LORNum * sizeof(float)))
			return false;
	}
	if (m_pPara->m_bRandomCorr)
	{
		if (!UFileReadSave::ReadFile((char*)michRand, m_pPara->m_strRandEstimatePath, LORNum * sizeof(float)))
			return false;
	}
	//==================================
	//dead time correction
	//==================================
	if (m_pPara->m_DeadTimeCorr)
	{
		printf("Using Dead Time Correction\n");
	}
	//==================================
	//calculate and read sensitivity map
	//==================================
	if (m_pPara->m_bCalcSenMap)
	{
		CalcSenMap(subsetNum, michNorm, michAttn);
		UFileReadSave::ReadFile((char*)senMap, m_pPara->m_strSenMapPath, nImgVoxelNum*subsetNum * sizeof(float));
	}
	else
	{
		if (!UFileReadSave::ReadFile((char*)senMap, m_pPara->m_strSenMapPath, nImgVoxelNum*subsetNum * sizeof(float)))
			return false;//read failure
	}
	//==================================
	//calculate crystal center position//节省计算量
	//==================================
	UVolume3D<double> *cryPosCenterArray = new UVolume3D<double>[crystalNumOneRing*ringNum];
	UVolume3D<double> cryPos[4];
	for (int ring = 0; ring < ringNum; ring++)
	{
		for (int cry = 0; cry < crystalNumOneRing; cry++)
		{
			pScan->GetCrystalPosition(cry, ring, cryPos);
			cryPosCenterArray[cry + ring * crystalNumOneRing] = (cryPos[0] + cryPos[1] + cryPos[2] + cryPos[3])*0.25;
		}
	}
	//==================================
	//recon
	//==================================
	for (int it = 0; it < iterNum; it++)
	{
		for (int sub = 0; sub < subsetNum; sub++)
		{
			//initial updating image
			memset(nimageMultiThreads, 0, nImgVoxelNum*threadNum * sizeof(float));
			#pragma omp parallel for num_threads(threadNum)
			for (int sl = 0; sl < sliceNum; sl++)
			{
				int currentThread = omp_get_thread_num();
				int LORId, ring1, ring2, cry1, cry2;
				pScan->GetRing1Ring2FromSlice(sl, ring1, ring2);
				//set ring difference//20190330
				if (abs(ring1 - ring2) > m_pPara->m_nRingDifference)
					continue;
				UOneLor ray(m_pPara->m_nImgWidth * 6);
				URayTracing rayTracing;
				float projection, SystemResponseCoff;
				for (int vi = sub; vi < viewNum; vi += subsetNum)
				{
					for (int bi = 0; bi < binNum; bi++)
					{
						//calculate crystal position
						LORId = bi + vi * binNum + sl * binNum*viewNum;
						pScan->GetRing1Ring2FromSlice(sl, ring1, ring2);
						pScan->GetCrystalIDInRingFromViewBin(vi, bi, cry1, cry2);
						rayTracing.GetOneLorResponse3D(ray, img, cryPosCenterArray[cry1+ring1*crystalNumOneRing], cryPosCenterArray[cry2+ring2*crystalNumOneRing]);//calculate SRM
						projection = 0;
						SystemResponseCoff = michAttn[LORId] * michNorm[LORId];//减少计算量
						//forward projection
						for (int i = 0; i < ray.GetVoxelNum(); i++)
						{
							projection += ray.pWeight[i] * img.m_pData[ray.pVoxelID[i]];
						}
						projection = (projection + michRand[LORId] / SystemResponseCoff);//正投影反投影，分子分母同时除以衰减因子和归一化因子
						if (projection > 0.0000001f)//正投影为0，则比值为1
							projection = mich[LORId] / projection;
						else
							projection = 1.0f;
						//back projection
						for (int i = 0; i < ray.GetVoxelNum(); i++)
						{
							nimageMultiThreads[ray.pVoxelID[i] + currentThread * nImgVoxelNum] += ray.pWeight[i] * projection;
						}
					}
				}
			}
			//更新图像在各线程上反投影的值合并到第一个线程，构成完整的更新图像
			for (int i = 1; i < threadNum; i++)
			{
				for (int j = 0; j < nImgVoxelNum; j++)
				{
					nimageMultiThreads[j] += nimageMultiThreads[j + i * nImgVoxelNum];
				}
			}
			//update
			for (int i = 0; i < nImgVoxelNum; i++)
			{
				if (senMap[i + sub * nImgVoxelNum] > 0.0000001f)
				{
					img.m_pData[i] *= nimageMultiThreads[i] / senMap[i + sub * nImgVoxelNum];
				}
			}
			//notify
			double step = double(sub+1 + it * subsetNum) / (iterNum*subsetNum);
			m_pNotify->step(step);
		}
	}
	img.SetOutOfViewZero();
	UReconPETFunction::CalcAverageActivity(img.m_pData, nImgVoxelNum, m_pPara->m_fTimeDuringScan);
	if (m_pPara->m_bDecayCorr)
		UReconPETFunction::DecayCorrection(img.m_pData, nImgVoxelNum, m_pPara->m_fTimeBeforeScan, m_pPara->m_fTimeDuringScan,m_pPara->m_fSourceHalfLife, m_pPara->m_fSourceBrachingRatio);
	if (m_pPara->m_bWellCounterCorr)
		UReconPETFunction::WellCounterCorrection(img.m_pData, nImgVoxelNum, m_pPara->m_fWellCounterIntercept, m_pPara->m_fWellCounterIntercept);
	UFileReadSave::SaveFile((char*)img.m_pData, m_pPara->m_strSavePath, nImgVoxelNum * sizeof(float));
	delete[]nimageMultiThreads;
	delete[]senMap;
	delete[]mich;
	delete[]michNorm;
	delete[]michAttn;
	delete[]michRand;
	delete[]cryPosCenterArray;

    m_pNotify->finish();

	return true;
}
/**********************************************************
Description:		Calculate the sensitivity map, and save the result.
Arguments:
{
	subsetNum:		Number of subset in iteration
	michNorm:		Normalization coefficient
	michAttn:		Attenuation coefficient
}
Return:				Status.
**********************************************************/
bool UReconPETOSEMRayTracingAlgor::CalcSenMap(int subsetNum, float *michNorm, float *michAttn)
{
	//get information
	const UPETScanner *pScan = m_pPara->m_pScan;
	const int crystalNumOneRing = pScan->GetCrystalNumOneRing();
	const int ringNum = pScan->GetRingNum();
	const int binNum = pScan->GetCrystalNumOneRing() - 1;
	const int viewNum = pScan->GetCrystalNumOneRing() / 2;
	const int sliceNum = pScan->GetRingNum()*pScan->GetRingNum();
	const int LORNum = binNum * viewNum * sliceNum;
	const int nImgVoxelNum = m_pPara->m_nImgWidth * m_pPara->m_nImgHeight * m_pPara->m_nImgDepth;
	const int threadNum = max(omp_get_num_procs() - 2,1);
	UImage img;//image
	img.SetCenter(0, 0, 0);
	img.SetVoxelSize(m_pPara->m_fVoxelSizeXY, m_pPara->m_fVoxelSizeXY, m_pPara->m_fVoxelSizeZ);
	img.SetVoxelNum(m_pPara->m_nImgWidth, m_pPara->m_nImgHeight, m_pPara->m_nImgDepth);
	//alloc
	float *senMapAllSubset = nullptr;//最终所有子集的灵敏度图像
	float *senMapMultiThreads = nullptr;//为了防止竞争，反投影时，每个线程分配一个图像缓冲区，最后将每个线程反投影的结果合并
	senMapAllSubset = new float[nImgVoxelNum*subsetNum];
	senMapMultiThreads = new float[nImgVoxelNum*threadNum];
	if (senMapMultiThreads == nullptr || senMapAllSubset == nullptr)
		return false;
	memset(senMapAllSubset, 0, nImgVoxelNum*subsetNum * sizeof(float));
	//calculate crystal center position
	UVolume3D<double> *cryPosCenterArray = new UVolume3D<double>[crystalNumOneRing*ringNum];
	UVolume3D<double> cryPos[4];
	for (int ring = 0; ring < ringNum; ring++)
	{
		for (int cry = 0; cry < crystalNumOneRing; cry++)
		{
			pScan->GetCrystalPosition(cry, ring, cryPos);
			cryPosCenterArray[cry + ring * crystalNumOneRing] = (cryPos[0] + cryPos[1] + cryPos[2] + cryPos[3])*0.25;
		}
	}
	//calculate
	for (int sub = 0; sub < subsetNum; sub++)
	{
		memset(senMapMultiThreads, 0, nImgVoxelNum*threadNum * sizeof(float));
		#pragma omp parallel for num_threads(threadNum)
		for (int sl = 0; sl < sliceNum; sl++)
		{
			int currentThread = omp_get_thread_num();
			int LORId, ring1, ring2, cry1, cry2;
			pScan->GetRing1Ring2FromSlice(sl, ring1, ring2);
			if (abs(ring1 - ring2) > m_pPara->m_nRingDifference)//set ring difference//20190330
				continue;
			UOneLor ray(m_pPara->m_nImgWidth * 6);
			URayTracing rayTracing;
			float SystemResponseCoff;
			for (int vi = sub; vi < viewNum; vi += subsetNum)
			{
				for (int bi = 0; bi < binNum; bi++)
				{
					//calculate crystal position
					LORId = bi + vi * binNum + sl * binNum*viewNum;
					pScan->GetRing1Ring2FromSlice(sl, ring1, ring2);
					pScan->GetCrystalIDInRingFromViewBin(vi, bi, cry1, cry2);
					rayTracing.GetOneLorResponse3D(ray, img, cryPosCenterArray[cry1 + ring1 * crystalNumOneRing], cryPosCenterArray[cry2 + ring2 * crystalNumOneRing]);//calculate SRM
					SystemResponseCoff = michAttn[LORId] * michNorm[LORId];//减少计算量
					//back projection
					for (int i = 0; i < ray.GetVoxelNum(); i++)
					{
						senMapMultiThreads[ray.pVoxelID[i] + currentThread*nImgVoxelNum] += ray.pWeight[i] * SystemResponseCoff;
					}
				}
			}
		}
		//合并该子集下，各线程的结果
		for (int i = 0; i < threadNum; i++)
		{
			for (int j = 0; j < nImgVoxelNum; j++)
			{
				senMapAllSubset[j + sub*nImgVoxelNum] += senMapMultiThreads[j + i * nImgVoxelNum];
			}
		}
		double step = double(sub + 1) / subsetNum*0.99;
		m_pNotify->step(step);
	}
	UFileReadSave::SaveFile((char*)senMapAllSubset, m_pPara->m_strSenMapPath, nImgVoxelNum*subsetNum * sizeof(float));
	delete[]senMapAllSubset;
	delete[]senMapMultiThreads;
	delete[]cryPosCenterArray;
	return true;
}
/**********************************************************
Description:		Calculate the forward projection.
Arguments:
{
	pScan:			A scanner object pointer
	img:			The forward projected image 
	path:			The path to save the result
	pNotify:		Object to record the progress
}
Return:				Status.
**********************************************************/
bool UReconPETOSEMRayTracingAlgor::ForwardProjection(const UPETScanner *pScan, UImage &img, std::string savePath, UProgressNotify *pNotify)
{
	//get information
	const int crystalNumOneRing = pScan->GetCrystalNumOneRing();
	const int ringNum = pScan->GetRingNum();
	const int binNum = pScan->GetCrystalNumOneRing() - 1;
	const int viewNum = pScan->GetCrystalNumOneRing() / 2;
	const int sliceNum = pScan->GetRingNum()*pScan->GetRingNum();
	const int LORNum = binNum * viewNum * sliceNum;
	const int nImgVoxelNum = img.GetVoxelNumInAll();
	const int threads = max(omp_get_num_procs() - 2, 1);
	//alloc
	float *image = img.m_pData;
	float *mich = new float[LORNum];
	memset(mich, 0, LORNum * sizeof(float));
	//calculate crystal center position
	UVolume3D<double> *cryPosCenterArray = new UVolume3D<double>[crystalNumOneRing*ringNum];
	UVolume3D<double> cryPos[4];
	for (int ring = 0; ring < ringNum; ring++)
	{
		for (int cry = 0; cry < crystalNumOneRing; cry++)
		{
			pScan->GetCrystalPosition(cry, ring, cryPos);
			cryPosCenterArray[cry + ring * crystalNumOneRing] = (cryPos[0] + cryPos[1] + cryPos[2] + cryPos[3])*0.25;
		}
	}
	//calculate in parallel
#pragma omp parallel for num_threads(threads)
	for (int vi = 0; vi < viewNum; vi++)
	{
		int LORId, ring1, ring2, cry1, cry2;
		int currentThread = omp_get_thread_num();
		UOneLor ray(img.GetVoxelNumX() * 5);
		URayTracing rayTracing;
		for (int sl = 0; sl < sliceNum; sl++)
		{
			for (int bi = 0; bi < binNum; bi++)
			{
				//calculate crystal position
				LORId = bi + vi * binNum + sl * binNum*viewNum;
				pScan->GetRing1Ring2FromSlice(sl, ring1, ring2);
				pScan->GetCrystalIDInRingFromViewBin(vi, bi, cry1, cry2);
				rayTracing.GetOneLorResponse3D(ray, img, cryPosCenterArray[cry1 + ring1 * crystalNumOneRing], cryPosCenterArray[cry2 + ring2 * crystalNumOneRing]);//calculate SRM
				//back projection
				for (int i = 0; i < ray.GetVoxelNum(); i++)
				{
					mich[LORId] += ray.pWeight[i] * image[ray.pVoxelID[i]];
				}
			}
		}
		if (pNotify != nullptr || pNotify != 0)
		{
			if (currentThread == 0)
			{
				double step = double(vi*threads) / viewNum;
				pNotify->step(step);
			}
		}
		
	}
	UFileReadSave::SaveFile((char*)mich, savePath, LORNum * sizeof(float));
	pNotify->finish();
	delete[]mich;
	delete[]cryPosCenterArray;
	return true;
}


