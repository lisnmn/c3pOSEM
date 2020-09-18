/********************************************************************

Copyright (C), 2019, All rights reserved

File Name     :    UProgressNotify.cpp
Description   :
History       :

<author>            <time>            <desc>
Ang Li             2019/6/8           create

********************************************************************/
#include "UProgressNotify.h"

UProgressNotify::UProgressNotify() :
m_fProgress(0)
{
}


UProgressNotify::~UProgressNotify()
{
}

void UProgressNotify::start()
{
	m_fProgress = 0.0;
}

void UProgressNotify::step(const double fStep)
{
	if (fStep >= 0 && fStep <= 1)
		m_fProgress = fStep;
	else if (fStep < 0)
		m_fProgress = 0;
	else//fStep>1
		m_fProgress = 1;
}

double UProgressNotify::getStep() const
{
	return m_fProgress;
}

void UProgressNotify::finish()
{
	m_fProgress = 2.0;
}
