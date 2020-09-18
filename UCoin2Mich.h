#ifndef __COIN_TO_MICH_H
#define __COIN_TO_MICH_H

#include <string>

//#include "UTypeDef.h"
#include "UAcquireDef.h"

class UPETScanner;
namespace UCoin2Mich
{
	int UCoin2Mich(UPETScanner &scanner, CoinStruct *coin, size_t coinNum, std::string michPath);
       
}


#endif
