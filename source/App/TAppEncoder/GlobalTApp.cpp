#include "GlobalTApp.h"

GlobalTApp::GlobalTApp()
{
	GlobalTApp::m_cTAppEncTop = NULL;
}

TAppEncTop *GlobalTApp::getTApp()
{
	if (m_cTAppEncTop == NULL) {
		m_cTAppEncTop = new TAppEncTop();
	}
	return m_cTAppEncTop;
}

GlobalTApp::~GlobalTApp()
{
}
TAppEncTop* GlobalTApp::m_cTAppEncTop = nullptr;
int GlobalTApp::frameCnt = -1;
int GlobalTApp::offset = 0;