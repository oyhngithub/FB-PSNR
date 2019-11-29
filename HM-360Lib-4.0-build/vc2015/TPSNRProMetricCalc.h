#pragma once
#include "TLib360/TGeometry.h"
class TPSNRProMetricCalc
{
public:
	Bool      m_bPSNRProEnabled;
	Double    m_dPSNRPro[3];

	// longtitude, latitude, weight
	CPos3D* m_frameInfo3D;
	// u, v
	IPos2D* m_frameInfo2D;

	Int       m_FeaturePoints;

	Int       m_outputBitDepth[MAX_NUM_CHANNEL_TYPE];         ///< bit-depth of output file
	Int       m_referenceBitDepth[MAX_NUM_CHANNEL_TYPE];      ///< bit-depth of reference file

	// Weight-to-Spherical
	Double* m_fErpWeight_Y;
	Double* m_fErpWeight_C;
	Double* m_fCubeWeight_Y;
	Double* m_fCubeWeight_C;
	Double* m_fEapWeight_Y;
	Double* m_fEapWeight_C;
	Double* m_fOctaWeight_Y;
	Double* m_fOctaWeight_C;
	Double* m_fIcoWeight_Y;
	Double* m_fIcoWeight_C;
#if SVIDEO_WSPSNR_SSP
	Double* m_fSspWeight_Y;
	Double* m_fSspWeight_C;
#endif
#if SVIDEO_ROTATED_SPHERE
	Double* m_fRspWeight_Y;
	Double* m_fRspWeight_C;
#endif

	Int     m_codingGeoType;
	Int     m_iCodingFaceWidth;
	Int     m_iCodingFaceHeight;
	Int     m_iChromaSampleLocType;
#if SVIDEO_ERP_PADDING
	Bool    m_bPERP;
#endif
#if SVIDEO_WSPSNR_E2E
	//for E2E WS-PSNR calculation;
#if !SVIDEO_E2E_METRICS
	TVideoIOYuv* m_pcTVideoIOYuvInputFile;  //note: reference;
	TGeometry* m_pRefGeometry;
	TGeometry* m_pRecGeometry;
	TComPicYuv* m_pcOrgPicYuv;
	TComPicYuv* m_pcRecPicYuv;             //in original geometry domain;
#endif
#if !SVIDEO_E2E_METRICS
	Int         m_iLastFrmPOC;
	UInt        m_temporalSubsampleRatio;
	Int         m_iInputWidth;
	Int         m_iInputHeight;
	ChromaFormat m_inputChromaFomat;
#endif
#endif


public:
	TPSNRProMetricCalc();
	~TPSNRProMetricCalc();
	Bool    getPSNRProEnabled() { return m_bPSNRProEnabled; }
	Void    setPSNRProEnabledFlag(Bool bEnabledFlag) { m_bPSNRProEnabled = bEnabledFlag; }
	Void    setOutputBitDepth(Int iOutputBitDepth[MAX_NUM_CHANNEL_TYPE]);
	Void    setReferenceBitDepth(Int iReferenceBitDepth[MAX_NUM_CHANNEL_TYPE]);
	Double* getPSNRPro() { return m_dPSNRPro; }
	//Void    sphSampoints(const std::string &cSphDataFile);
	//Void    sphToCart(CPos2D*, CPos3D*);
	//Void    createTable(TGeometry *pcCodingGeomtry);
	Void    xCalculatePSNRPro(TComPicYuv* pcOrgPicYuv, TComPicYuv* pcPicD);

	Void    setCodingGeoInfo2(SVideoInfo& sRefVideoInfo, SVideoInfo& sRecVideoInfo, InputGeoParam* pInGeoParam);
	Void    xCalculateE2EPSNRPro(TComPicYuv* pcRecPicYuv, TComPicYuv* pcOrigPicYuv);

#if !SVIDEO_ROUND_FIX
	inline Int round(POSType t) { return (Int)(t + (t >= 0 ? 0.5 : -0.5)); };
#endif
};

