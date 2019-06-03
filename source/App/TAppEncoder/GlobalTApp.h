#pragma
#include "TAppEncTop.h"

class GlobalTApp {
public:
	static TAppEncTop *getTApp();
	static int frameCnt;
	static int offset;
	virtual ~GlobalTApp();
	GlobalTApp();
	static TAppEncTop *m_cTAppEncTop;
};
