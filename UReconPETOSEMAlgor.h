/********************************************************************

Copyright (C), 2019, All rights reserved

File Name     :    UPETReconOSEMAlgor.h
Description   :
History       :

<author>            <time>            <desc>
Ang Li             2019/6/8           create

********************************************************************/
#ifndef __RECON_PET_OSEM_ALGOR_H
#define __RECON_PET_OSEM_ALGOR_H

#include <string>

class UPETScanner;
class UImage;
class UProgressNotify;

class UReconPETOSEMPara
{
public:
	UReconPETOSEMPara();
	~UReconPETOSEMPara();


	UPETScanner					*m_pScan;

	int                         m_nImgWidth;//image pixel num
	int                         m_nImgHeight;
	int                         m_nImgDepth;
	float                       m_fVoxelSizeXY;//image pixel size
	float                       m_fVoxelSizeZ;

	int                         m_nIterNum;
	int                         m_nSubsetNum;
	int							m_nRingDifference;//max ring difference slice in mich will be used

	float						m_fTimeBeforeScan;
	float						m_fTimeDuringScan;
	float						m_fSourceHalfLife;
	float						m_fSourceBrachingRatio;
	float						m_fWellCounterSlope;// well counter factor
	float						m_fWellCounterIntercept;
    
	std::string                 m_strMichPath;//mich file path
	std::string                 m_strSenMapPath;//sensitivity map path

	std::string                 m_strNormFactorPath;//normalization file path
	std::string					m_strAttnFactorPath;//attenuation file path
	std::string					m_strRandEstimatePath;//random estimate file path
	std::string                 m_strSavePath;

	bool						m_bCalcSenMap;
	bool                        m_bRandomCorr;
	bool                        m_bAttnCorr;
	bool                        m_bScatterCorr;
	bool                        m_bNormCorr;
	bool						m_bWellCounterCorr;
	bool						m_bDecayCorr;
	bool						m_DeadTimeCorr;
    
    //bool                        m_bGaussionFWHM;
    //int                         m_nGaussion;
private:
	UReconPETOSEMPara(const UReconPETOSEMPara&){};
	UReconPETOSEMPara& operator=(const UReconPETOSEMPara&){};
};

class UReconPETOSEMRayTracingAlgor
{
public:
	UReconPETOSEMRayTracingAlgor();
	virtual ~UReconPETOSEMRayTracingAlgor();
	//interface
	/**********************************************************
	Description:		Initialize the PET reconstruction parameter.
	Arguments:
	{
	pPara:			Reconstraction parameter.
	pNotify:		Object to record the progress.
	}
	Return:				Status.notify
	**********************************************************/
	bool Initial(UReconPETOSEMPara * pPara, UProgressNotify * pNotify);
	/**********************************************************
	Description:		Reconstruct and save the PET image.
	Arguments:
	{

	}
	Return:				Status.
	**********************************************************/
	bool DoRecon();
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
	bool CalcSenMap(int subsetNum, float *michNorm, float *michAttn);//calculate sensitivity map
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
	static bool ForwardProjection(const UPETScanner *pScan, UImage &img, std::string savePath, UProgressNotify *pNotify=nullptr);
private:
	UReconPETOSEMRayTracingAlgor(const UReconPETOSEMRayTracingAlgor&){};
	UReconPETOSEMRayTracingAlgor& operator=(const UReconPETOSEMRayTracingAlgor&){};
private:
	UReconPETOSEMPara * m_pPara;
	UProgressNotify * m_pNotify;
};

#endif // !__PET_RECON_OSEM_ALOGR_H
