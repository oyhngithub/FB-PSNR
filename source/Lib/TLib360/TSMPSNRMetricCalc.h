#pragma once
#include "TLib360/TGeometry.h"

class TSMPSNRMetric
{
public:
	Bool      m_bSMPSNREnabled;
	Double    m_dSMPSNR[3];

	// m_response stores response wseight
	CPos3D*   m_response;
	// m_pCart3D stores radius weights
	CPos3D*   m_pCart3D;
	IPos2D*   m_fpTable;
	Int       m_iSphNumPoints;
	Int       m_iFeaturePoints;

	Int       m_outputBitDepth[MAX_NUM_CHANNEL_TYPE];         ///< bit-depth of output file
	Int       m_referenceBitDepth[MAX_NUM_CHANNEL_TYPE];      ///< bit-depth of reference file


public:
	TSMPSNRMetric();
	~TSMPSNRMetric();
	double min, max;
	Bool    getSMPSNREnabled() { return m_bSMPSNREnabled; }
	Void    setSMPSNREnabledFlag(Bool bEnabledFlag) { m_bSMPSNREnabled = bEnabledFlag; }
	Void    setOutputBitDepth(Int iOutputBitDepth[MAX_NUM_CHANNEL_TYPE]);
	Void    setReferenceBitDepth(Int iReferenceBitDepth[MAX_NUM_CHANNEL_TYPE]);
	Double* getSMPSNR() { return m_dSMPSNR; }
	//Void    sphSampoints(const std::string &cSphDataFile);
	Void    sphSampoints(const std::string &cSphDataFile, const std::string &);
	Void    sphToCart(CPos2D*, CPos3D*);
	Void    createTable(TGeometry *pcCodingGeomtry);
	Void    xCalculateSMPSNR(TComPicYuv* pcOrgPicYuv, TComPicYuv* pcPicD);

#if !SVIDEO_ROUND_FIX
	inline Int round(POSType t) { return (Int)(t + (t >= 0 ? 0.5 : -0.5)); };
#endif
};

