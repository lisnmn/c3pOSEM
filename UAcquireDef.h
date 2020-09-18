/********************************************************************

  Copyright (C), 2019, All rights reserved

  File Name     :    UAcquireDef.h
  Description   :
  History       :

  <author>            <time>            <desc>
  Lingyiqing          2019/1/1          create

********************************************************************/

#ifndef __UACQUIREDEF_H
#define __UACQUIREDEF_H

#include "UTypeDef.h"
#include <string>
using namespace std;

#pragma pack (push, 1)

/* Frame received from BDM */
typedef struct _DataFrameV2
{
    unsigned char   nHeadAndDU;
    unsigned char   nBDM;
    unsigned char   nTime[8];
    unsigned char   X;
    unsigned char   Y;
    unsigned char   Energy[2];
    char            nTemperatureInt;
    unsigned char   nTemperatureAndTail;
}DataFrameV2;

/* Used in the old version */
typedef struct _ResolveFrame
{
    ULONGLONG       nTime;
    BYTE            X;
    BYTE            Y;
    WORD            nEnergy;
    BYTE            PMT;
    BYTE            BDM;
    float           fTemperature;
}ResolveFrame;

/* Used in clinical PET */
typedef struct _SamplesStruct
{
    unsigned short globalBDMIndex;
    unsigned short localCrystalIndex;
    double         timevalue[8];
}SamplesStruct;

/* Temp structure. Local -> Global */
typedef struct _TempSinglesStruct
{
    unsigned short globalBDMIndex;
    unsigned short localDUIndex;
    unsigned short localCrystalIndex;
    float          energy;
    double         timevalue;
}TempSinglesStruct;

/* Used in the new version */
typedef struct _SinglesStruct
{
    unsigned int   globalCrystalIndex;
    float          energy;
    double         timevalue;
}SinglesStruct;

/* Used in listmode */
typedef struct _CoinStruct
{
    SinglesStruct nCoinStruct[2];
}CoinStruct;

#pragma pack(pop)

/* The size of crystal in one DU */
const int CRYSTAL_SIZE = 13;
/* The size of position in one DU */
const int POSITION_SIZE = 256;
/* The num of DU */
const int DU_NUM = 2;
/* The num of BDM */
const int BDM_NUM = 12;
/* The length of time window */
const int TIME_WINDOW = 600;
/*************************************/
/* To be normalized */
/* The minimum energy, 0.4f*/
const float MIN_ENERGY = 2000;
/* The maximum energy, 0.8f*/
const float MAX_ENERGY = 3000;
/*************************************/
/* The num of a common data buffer, 150MB for DataFrameV2 */
const int FRAME_NUM_ONE_BUFFER = 64 * 1024 * 150 ;
/* The num of total packages in one channel, 150MB */
const int PACKAGE_NUM = 1024 * 150;
/* The size of one package, 1KB */
const int PACKAGE_SIZE = 1024;
/* The size of one frame */
const int FRAME_SIZE = sizeof(DataFrameV2);
/* The num of frames in one package */
const int FRAME_NUM_ONE_PKG = PACKAGE_SIZE / FRAME_SIZE;
/* The version of the udp frame */
const int UDP_FRAME_VERSION = 2;

/* To save raw data received from BDM Udp pkg */
static const string strRawDataPath = "../PET/RawData/";
/* To save coin data */
static const string strCoinDataPath = "../PET/CoinData/Coin.dat";
/* To save mich data */
static const string strMichDataPath = "../PET/MichData/Mich.dat";
/* To save sort data */
static const string strSortDataPath = "../PET/SortData/Sort.dat";
/* To save CT data received from Tcp pkg */
static const string strCTDataPath = "../CT/";
/* To save preset position table */
static const string strPositionTablePath = "../Preset/Position.dat";
#endif

