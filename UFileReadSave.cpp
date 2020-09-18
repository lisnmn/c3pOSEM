/********************************************************************

Copyright (C), 2019, All rights reserved

File Name     :    UFileReadSave.cpp
Description   :
History       :

<author>            <time>            <desc>
Ang Li             2019/6/8           create

********************************************************************/
#include "UFileReadSave.h"
#include <fstream>
#include <iostream>
#include <string>

using namespace std;


/**********************************************************
Description:		Read the binary file to an memory.
Arguments:
{
	pDst:			The destination memory.
	strFilePath:	File path.
	fileByteNum:	The Bytes will be read.
}
Return:				Status.
**********************************************************/
bool UFileReadSave::ReadFile(char * pDst, std::string & strFilePath, size_t fileByteNum)
{
	ifstream projFile(strFilePath.c_str(), ios::binary);
	if (!projFile.is_open())
	{
		cout << strFilePath << " is losted!" << endl;//file losted
		return false;
	}
	projFile.seekg(0, ios::end);
	size_t fileLen = projFile.tellg();
	if (fileLen < fileByteNum)//file too small
	{
		cout << strFilePath << " only has " << fileLen << "bytes which is less than required " << fileByteNum << " bytes." << endl;
	}
	projFile.seekg(0, ios::beg);
	projFile.read(pDst, fileByteNum);
	projFile.close();
	return true;
}
/**********************************************************
Description:		Write the binary file to the disk.
Arguments:
{
	pSrc:			The source data.
	strFilePath:	File path.
	fileByteNum:	The Bytes will be saved.
}
Return:				Status.
**********************************************************/
bool UFileReadSave::SaveFile(char * pSrc, std::string & strFilePath, size_t fileByteNum)
{
	ofstream projFile(strFilePath.c_str(), ios::binary);
	if (!projFile.is_open())
	{
		cout << strFilePath << " is losted!" << endl;//file losted
		return false;
	}
	projFile.write(pSrc, fileByteNum);
	projFile.close();
	return true;
}