/********************************************************************

Copyright (C), 2019, All rights reserved

File Name     :    UFileReadSave.h
Description   :
History       :

<author>            <time>            <desc>
Ang Li             2019/6/8           create

********************************************************************/
#ifndef __FILE_READ_SAVE_H
#define __FILE_READ_SAVE_H

#include <string>

class UFileReadSave
{
public:
	UFileReadSave(void) {};//inline
	~UFileReadSave(void) {};
public:
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
	static bool ReadFile(char *pSrc, std::string &strFilePath, size_t fileByteNum);
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
	static bool SaveFile(char *pSrc, std::string &strFilePath, size_t fileByteNum);
private:
	UFileReadSave(const UFileReadSave&){};
	UFileReadSave& operator=(const UFileReadSave&){};
};


#endif
