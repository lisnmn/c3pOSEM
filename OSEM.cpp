//
// Created by svf on 2020/9/14.
//

#include "OSEM.h"
#include <iostream>


bool OSEM::start(const char* path) {
    // 1.Read in paras.
    int error = mIniFile.Load(path);
    if(error){
        std::cout << "Error while reading " << path << std::endl;
        return false;
    }

    initScanner();
    this->mOSEMPara.m_pScan = &mScanner;
    std::string sectionName = "parameters";
    mIniFile.GetIntValue(sectionName, "int_imageWidth", &mOSEMPara.m_nImgWidth);
    mIniFile.GetIntValue(sectionName, "int_imageHeight", &mOSEMPara.m_nImgHeight);
    mIniFile.GetIntValue(sectionName, "int_imageDepth", &mOSEMPara.m_nImgDepth);

    double doubleTemp;
    mIniFile.GetDoubleValue(sectionName, "double_voxelSizeXY", &doubleTemp);
    mOSEMPara.m_fVoxelSizeXY = (float)doubleTemp;
    mIniFile.GetDoubleValue(sectionName, "double_voxelSizeZ", &doubleTemp);
    mOSEMPara.m_fVoxelSizeZ = (float)doubleTemp;
    mIniFile.GetIntValue(sectionName, "int_iterNum", &mOSEMPara.m_nIterNum);
    mIniFile.GetIntValue(sectionName, "int_subsetNum", &mOSEMPara.m_nSubsetNum);
    mIniFile.GetBoolValue(sectionName, "bool_calSenMap", &mOSEMPara.m_bCalcSenMap);

    mIniFile.GetIntValue("reconstructSetting", "ring_difference", &mOSEMPara.m_nRingDifference);

    sectionName = "reconstructInfo";
    mIniFile.GetStringValue(sectionName, "mich_path", &mOSEMPara.m_strMichPath);
    mIniFile.GetStringValue(sectionName, "save_path", &mOSEMPara.m_strSavePath);
    if(mOSEMPara.m_bCalcSenMap) {
        std::string coinPath;
        mIniFile.GetStringValue(sectionName, "coin_path", &coinPath);
        mOSEMPara.m_strSenMapPath = coinPath + "/bed0-senMap.dat";
    }

    // 2.Run algorithm IN THREAD.
    m_pThreadOSEM = std::make_shared<std::thread>(std::thread([this]{
        this->mOSEM.Initial(&this->mOSEMPara, &this->mNotify);
        this->mOSEM.DoRecon();
    }));
    m_pThreadOSEM->detach();
    return true;
}

bool OSEM::stop() {
    m_pThreadOSEM.~shared_ptr();
    return true;
}

double OSEM::progress() {
    return mNotify.getStep();
}

void OSEM::initScanner() {
    //mScanner.InitD180Scanner();
    std::string sectionName = "scannerPETSetting";
    mIniFile.GetIntValue(sectionName, "crystal_num_X", &mScanner.crystalNum.x);
    mIniFile.GetIntValue(sectionName, "crystal_num_Y", &mScanner.crystalNum.y);
    mIniFile.GetIntValue(sectionName, "crystal_num_Z", &mScanner.crystalNum.z);
    mIniFile.GetDoubleValue(sectionName, "crystal_size_X", &mScanner.crystalSize.x);
    mIniFile.GetDoubleValue(sectionName, "crystal_size_Y", &mScanner.crystalSize.y);
    mIniFile.GetDoubleValue(sectionName, "crystal_size_Z", &mScanner.crystalSize.z);
    mIniFile.GetDoubleValue(sectionName, "crystal_pitch_X", &mScanner.crystalPitch.x);
    mIniFile.GetDoubleValue(sectionName, "crystal_pitch_Y", &mScanner.crystalPitch.y);
    mIniFile.GetDoubleValue(sectionName, "crystal_pitch_Z", &mScanner.crystalPitch.z);

    mIniFile.GetIntValue(sectionName, "block_num_X", &mScanner.blockNum.x);
    mIniFile.GetIntValue(sectionName, "block_num_Y", &mScanner.blockNum.y);
    mIniFile.GetIntValue(sectionName, "block_num_Z", &mScanner.blockNum.z);
    mIniFile.GetDoubleValue(sectionName, "block_size_X", &mScanner.blockSize.x);
    mIniFile.GetDoubleValue(sectionName, "block_size_Y", &mScanner.blockSize.y);
    mIniFile.GetDoubleValue(sectionName, "block_size_Z", &mScanner.blockSize.z);
    mIniFile.GetDoubleValue(sectionName, "block_pitch_X", &mScanner.blockPitch.x);
    mIniFile.GetDoubleValue(sectionName, "block_pitch_Y", &mScanner.blockPitch.y);
    mIniFile.GetDoubleValue(sectionName, "block_pitch_Z", &mScanner.blockPitch.z);

    mIniFile.GetIntValue(sectionName, "module_num_X", &mScanner.moduleNum.x);
    mIniFile.GetIntValue(sectionName, "module_num_Y", &mScanner.moduleNum.y);
    mIniFile.GetIntValue(sectionName, "module_num_Z", &mScanner.moduleNum.z);
    mIniFile.GetDoubleValue(sectionName, "module_size_X", &mScanner.moduleSize.x);
    mIniFile.GetDoubleValue(sectionName, "module_size_Y", &mScanner.moduleSize.y);
    mIniFile.GetDoubleValue(sectionName, "module_size_Z", &mScanner.moduleSize.z);
    mIniFile.GetDoubleValue(sectionName, "module_pitch_X", &mScanner.modulePitch.x);
    mIniFile.GetDoubleValue(sectionName, "module_pitch_Y", &mScanner.modulePitch.y);
    mIniFile.GetDoubleValue(sectionName, "module_pitch_Z", &mScanner.modulePitch.z);

    mIniFile.GetIntValue(sectionName, "panel_num", &mScanner.panelNum);
    mIniFile.GetDoubleValue(sectionName, "panel_size_X", &mScanner.panelSize.x);
    mIniFile.GetDoubleValue(sectionName, "panel_size_Y", &mScanner.panelSize.y);
    mIniFile.GetDoubleValue(sectionName, "panel_size_Z", &mScanner.panelSize.z);
    mIniFile.GetDoubleValue(sectionName, "panel_pitch_X", &mScanner.panelPitch.x);
    mIniFile.GetDoubleValue(sectionName, "panel_pitch_Y", &mScanner.panelPitch.y);
    mIniFile.GetDoubleValue(sectionName, "panel_pitch_Z", &mScanner.panelPitch.z);

    mIniFile.GetIntValue(sectionName, "crystal_offset", &mScanner.crystalClockwiseOffset);
    mIniFile.GetDoubleValue(sectionName, "scanner_radius", &mScanner.scannerRadius);
}