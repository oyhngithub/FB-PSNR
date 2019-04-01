#ifndef __TFBPSNRCALC__
#define __TFBPSNRCALC__
#include "TLib360/TGeometry.h"

class TFBPSNRMetric
{
private:
	Bool      m_bSPSNREnabled;
	Double    m_dFBPSNR[3];

	CPos2D*   m_pCart2D;
	IPos2D*   m_fpTable;
	Int       m_iSphNumPoints;

	Int       m_outputBitDepth[MAX_NUM_CHANNEL_TYPE];         ///< bit-depth of output file
	Int       m_referenceBitDepth[MAX_NUM_CHANNEL_TYPE];      ///< bit-depth of reference file
public:
	TFBPSNRMetric();
	virtual ~TFBPSNRMetric();
	Bool    getSPSNREnabled() { return m_bSPSNREnabled; }
	Void    setSPSNREnabledFlag(Bool bEnabledFlag) { m_bSPSNREnabled = bEnabledFlag; }
	Void    setOutputBitDepth(Int iOutputBitDepth[MAX_NUM_CHANNEL_TYPE]);
	Void    setReferenceBitDepth(Int iReferenceBitDepth[MAX_NUM_CHANNEL_TYPE]);
	Double* getFBPSNR() { return m_dFBPSNR; }
	Void    sphSampoints(const std::string &cSphDataFile);
	Void    sphToCart(CPos2D*, CPos3D*);
	Void    createTable(TGeometry *pcCodingGeomtry);
	Void    xCalculateFBPSNR(TComPicYuv* pcOrgPicYuv, TComPicYuv* pcPicD);

#if !SVIDEO_ROUND_FIX
	inline Int round(POSType t) { return (Int)(t + (t >= 0 ? 0.5 : -0.5)); };
#endif

};

#endif