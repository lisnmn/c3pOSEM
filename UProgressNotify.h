/********************************************************************

Copyright (C), 2019, All rights reserved

File Name     :    UProgressNotify.h
Description   :
History       :

<author>            <time>            <desc>
Ang Li             2019/6/8           create

********************************************************************/
#ifndef __PROGRESS_NOTIFY_H
#define __PROGRESS_NOTIFY_H

class UProgressNotify
{
public:
	UProgressNotify();
	virtual ~UProgressNotify();
	void start();
	void step(const double fStep);
	double getStep()const;
	void finish();

private:
	UProgressNotify(const UProgressNotify&){};
	UProgressNotify& operator=(const UProgressNotify&){};
private:
	double m_fProgress;
};

#endif
