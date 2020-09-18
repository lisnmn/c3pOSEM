#include "UCoin2Mich.h"
#include "UPETScanner.h"
#include "UFileReadSave.h"

namespace UCoin2Mich
{
	
    int UCoin2Mich(UPETScanner &scanner, CoinStruct *coin, size_t coinNum, std::string michPath)
	{
		size_t lorNum = scanner.GetLORNum();
		int ringNum = scanner.GetRingNum();
		int cryNumOneRing = scanner.GetCrystalNumOneRing();
		float* mich = nullptr;
		mich = new float[lorNum];
		if (mich == nullptr)
			return 1;
		for (int i = 0; i < lorNum; i++)
			mich[i] = 0.0f;
		for (int i = 0; i < coinNum; i++)
		{
			int crystalId1 = coin[i].nCoinStruct[0].globalCrystalIndex;
			int crystalId2 = coin[i].nCoinStruct[1].globalCrystalIndex;
			int lorId = 0;
			scanner.GetLORIDFromCrystalID(crystalId1, crystalId2, lorId);
			mich[lorId]++;
		}
		bool status = UFileReadSave::SaveFile((char*)mich, michPath, lorNum*sizeof(float));
		delete[] mich;
		if (status)
			return 0;
		else
			return 2;
	}
       

} // namespace UCoin2Mich
