#pragma once
#include "TLib360/TGeometry.h"

class TSMPSNRVMetric
{
public:
	Bool      m_bSMPSNRVEnabled;
	Double    m_dSMPSNRV[3];

	// change Cart2D to Cart3D for storing response
	CPos3D*   m_pCart3D;
	IPos2D*   m_fpTable;
	Int       m_iSphNumPoints;

	Int       m_outputBitDepth[MAX_NUM_CHANNEL_TYPE];         ///< bit-depth of output file
	Int       m_referenceBitDepth[MAX_NUM_CHANNEL_TYPE];      ///< bit-depth of reference file

	// WS-PSNR member
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
	TVideoIOYuv *m_pcTVideoIOYuvInputFile;  //note: reference;
	TGeometry   *m_pRefGeometry;
	TGeometry   *m_pRecGeometry;
	TComPicYuv  *m_pcOrgPicYuv;
	TComPicYuv  *m_pcRecPicYuv;             //in original geometry domain;
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
	TSMPSNRVMetric();
	~TSMPSNRVMetric();
	Bool    getSMPSNRVEnabled() { return m_bSMPSNRVEnabled; }
	Void    setSMPSNRVEnabledFlag(Bool bEnabledFlag) { m_bSMPSNRVEnabled = bEnabledFlag; }
	Void    setOutputBitDepth(Int iOutputBitDepth[MAX_NUM_CHANNEL_TYPE]);
	Void    setReferenceBitDepth(Int iReferenceBitDepth[MAX_NUM_CHANNEL_TYPE]);
	Double* getSMPSNRV() { return m_dSMPSNRV; }
	//Void    sphSampoints(const std::string &cSphDataFile);
	//Void    sphToCart(CPos2D*, CPos3D*);
	//Void    createTable(TGeometry *pcCodingGeomtry);
	Void    xCalculateSMPSNRV(TComPicYuv* pcOrgPicYuv, TComPicYuv* pcPicD);

#if !SVIDEO_ROUND_FIX
	inline Int round(POSType t) { return (Int)(t + (t >= 0 ? 0.5 : -0.5)); };
#endif
};
