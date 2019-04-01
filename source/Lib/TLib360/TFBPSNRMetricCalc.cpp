#include "TFBPSNRMetricCalc.h"
//#include "TLibEncoder/TEncRateCtrl.cpp"
#include "TAppEncHelper360/TExt360AppEncCfg.h"

// ====================================================================================================================
// DIY
#include <sstream>
#include <fstream>
#include <iostream>


TFBPSNRMetric::TFBPSNRMetric()
	: m_bSPSNREnabled(false)
	, m_pCart2D(NULL)
	, m_fpTable(NULL)
{
	m_dFBPSNR[0] = m_dFBPSNR[1] = m_dFBPSNR[2] = 0;
}

TFBPSNRMetric::~TFBPSNRMetric()
{
	if (m_pCart2D)
	{
		free(m_pCart2D); m_pCart2D = NULL;
	}
	if (m_fpTable)
	{
		free(m_fpTable); m_fpTable = NULL;
	}
}

Void TFBPSNRMetric::setOutputBitDepth(Int iOutputBitDepth[MAX_NUM_CHANNEL_TYPE])
{
	for (Int i = 0; i < MAX_NUM_CHANNEL_TYPE; i++)
	{
		m_outputBitDepth[i] = iOutputBitDepth[i];
	}
}

Void TFBPSNRMetric::setReferenceBitDepth(Int iReferenceBitDepth[MAX_NUM_CHANNEL_TYPE])
{
	for (Int i = 0; i < MAX_NUM_CHANNEL_TYPE; i++)
	{
		m_referenceBitDepth[i] = iReferenceBitDepth[i];
	}
}

Void TFBPSNRMetric::sphSampoints(const std::string &cSphDataFile)
{
	FILE* fp = fopen(TExt360AppEncCfg::m_featureFileName.c_str(), "r");
	int x, y;
	x = y = 0;
	int edgeNumbers, featureNumbers;
	// read longtitude, latitude
	if (fscanf(fp, "%d %d", &edgeNumbers, &featureNumbers) != 2)
	{
		printf("SphData file does not exist.\n");
		exit(EXIT_FAILURE);
		fclose(fp);
	}
	m_iSphNumPoints = edgeNumbers + featureNumbers;

	m_pCart2D = (CPos2D*)malloc(sizeof(CPos2D)*(m_iSphNumPoints));
	memset(m_pCart2D, 0, (sizeof(CPos2D) * m_iSphNumPoints));
	for (Int z = 0; z < m_iSphNumPoints; z++)
	{
		if (fscanf(fp, "%lf %lf", &m_pCart2D[z].y, &m_pCart2D[z].x) != 2)
		{
			printf("Format error SphData in sphSampoints().\n");
			exit(EXIT_FAILURE);
		}
	}
	fclose(fp);
}

void TFBPSNRMetric::sphToCart(CPos2D* sph, CPos3D* out)
{
	POSType fLat = (POSType)(sph->x*S_PI / 180.0);
	POSType fLon = (POSType)(sph->y*S_PI / 180.0);

	out->x = ssin(fLon) * scos(fLat);
	out->y = ssin(fLat);
	out->z = -scos(fLon) * scos(fLat);
}

void TFBPSNRMetric::createTable(TGeometry *pcCodingGeomtry)
{
	Int iNumPoints = m_iSphNumPoints;
	CPos2D In2d;
	CPos3D Out3d;
	SPos posIn, posOut;
	m_fpTable = (IPos2D*)malloc(iNumPoints * sizeof(IPos2D));

	for (Int np = 0; np < iNumPoints; np++)
	{
		In2d.x = m_pCart2D[np].x;
		In2d.y = m_pCart2D[np].y;

		//get cartesian coordinates
		sphToCart(&In2d, &Out3d);
		posIn.x = Out3d.x; posIn.y = Out3d.y; posIn.z = Out3d.z;
		assert(posIn.x < 1 && posIn.x > -1 && posIn.y > -1 && posIn.y < 1 && posIn.z < 1 && posIn.z > -1);
		pcCodingGeomtry->map3DTo2D(&posIn, &posOut);

		posOut.x = (POSType)TGeometry::round(posOut.x);
		posOut.y = (POSType)TGeometry::round(posOut.y);
		IPos tmpPos;
		tmpPos.faceIdx = posOut.faceIdx;
		tmpPos.u = (Int)(posOut.x);
		tmpPos.v = (Int)(posOut.y);
		pcCodingGeomtry->clamp(&tmpPos);
		pcCodingGeomtry->geoToFramePack(&tmpPos, &m_fpTable[np]);
	}
	// ====================================================================================================================
	// DIY
	//std::ofstream outfile("CountData.txt", std::ios::out | std::ios::binary);
	//if (!outfile.is_open())
	//{
	   // std::cout << " the file open fail" << std::endl;
	   // exit(1);
	//}
	//point2DCount countVector = pcCodingGeomtry->getPoint2DCount();
	//for (int i = 0; i < countVector.size(); ++i)
	//{
	   // for (int j = 0; j < countVector[i].size(); ++j)
	   // {
		  //  outfile << countVector[i][j] << " ";
	   // }
	   // outfile << "\r\n";
	//}
	//outfile.close();

	// ====================================================================================================================
}

// pass the pics before and after the encode
// m_fpTable stores all the data that already transform into different projections
Void TFBPSNRMetric::xCalculateFBPSNR(TComPicYuv* pcOrgPicYuv, TComPicYuv* pcPicD)
{
	Int iNumPoints = m_iSphNumPoints;
	Int iBitDepthForPSNRCalc[MAX_NUM_CHANNEL_TYPE];
	Int iReferenceBitShift[MAX_NUM_CHANNEL_TYPE];
	Int iOutputBitShift[MAX_NUM_CHANNEL_TYPE];
	iBitDepthForPSNRCalc[CHANNEL_TYPE_LUMA] = std::max(m_outputBitDepth[CHANNEL_TYPE_LUMA], m_referenceBitDepth[CHANNEL_TYPE_LUMA]);
	iBitDepthForPSNRCalc[CHANNEL_TYPE_CHROMA] = std::max(m_outputBitDepth[CHANNEL_TYPE_CHROMA], m_referenceBitDepth[CHANNEL_TYPE_CHROMA]);
	iReferenceBitShift[CHANNEL_TYPE_LUMA] = iBitDepthForPSNRCalc[CHANNEL_TYPE_LUMA] - m_referenceBitDepth[CHANNEL_TYPE_LUMA];
	iReferenceBitShift[CHANNEL_TYPE_CHROMA] = iBitDepthForPSNRCalc[CHANNEL_TYPE_CHROMA] - m_referenceBitDepth[CHANNEL_TYPE_CHROMA];
	iOutputBitShift[CHANNEL_TYPE_LUMA] = iBitDepthForPSNRCalc[CHANNEL_TYPE_LUMA] - m_outputBitDepth[CHANNEL_TYPE_LUMA];
	iOutputBitShift[CHANNEL_TYPE_CHROMA] = iBitDepthForPSNRCalc[CHANNEL_TYPE_CHROMA] - m_outputBitDepth[CHANNEL_TYPE_CHROMA];

	memset(m_dFBPSNR, 0, sizeof(Double) * 3);
	TComPicYuv &picd = *pcPicD;
	Double SSDspsnr[3] = { 0, 0 ,0 };

	for (Int chan = 0; chan < pcPicD->getNumberValidComponents(); chan++)
	{
		const ComponentID ch = ComponentID(chan);
		const Pel*  pOrg = pcOrgPicYuv->getAddr(ch);
		const Int   iOrgStride = pcOrgPicYuv->getStride(ch);
		const Pel*  pRec = picd.getAddr(ch);
		const Int   iRecStride = picd.getStride(ch);

		for (Int np = 0; np < iNumPoints; np++)
		{
			if (!chan)
			{
				Int x_loc = (Int)(m_fpTable[np].x);
				Int y_loc = (Int)(m_fpTable[np].y);
				Intermediate_Int iDifflp = (pOrg[x_loc + (y_loc*iOrgStride)] << iReferenceBitShift[toChannelType(ch)]) - (pRec[x_loc + (y_loc*iRecStride)] << iOutputBitShift[toChannelType(ch)]);
				SSDspsnr[chan] += iDifflp * iDifflp;
			}
			else
			{
				Int x_loc = Int(m_fpTable[np].x >> pcPicD->getComponentScaleX(COMPONENT_Cb));
				Int y_loc = Int(m_fpTable[np].y >> pcPicD->getComponentScaleY(COMPONENT_Cb));
				Intermediate_Int iDifflp = (pOrg[x_loc + (y_loc*iOrgStride)] << iReferenceBitShift[toChannelType(ch)]) - (pRec[x_loc + (y_loc*iRecStride)] << iOutputBitShift[toChannelType(ch)]);
				SSDspsnr[chan] += iDifflp * iDifflp;
			}
		}
	}

	for (Int ch_indx = 0; ch_indx < pcPicD->getNumberValidComponents(); ch_indx++)
	{
		const ComponentID ch = ComponentID(ch_indx);
		const Int maxval = 255 << (iBitDepthForPSNRCalc[toChannelType(ch)] - 8);

		Double fReflpsnr = Double(iNumPoints)*maxval*maxval;
		m_dFBPSNR[ch_indx] = (SSDspsnr[ch_indx] ? 10.0 * log10(fReflpsnr / (Double)SSDspsnr[ch_indx]) : 999.99);
	}
}

// create 2D table 
//void TEncGOP::createTable(TGeometry *pcCodingGeomtry)
//{
//	//-- Added by Ma Guilong
//
//
//	// 3D - cartesian - faceID, x, y
//
//	CPos2D In2d;
//	CPos3D Out3d;
//	SPos posIn, posOut;
//	m_fpTable = (IPos2D*)malloc(iNumPoints * sizeof(IPos2D));
//
//	for (Int np = 0; np < iNumPoints; np++)
//	{
//		In2d.x = m_pCart2D[np].x;
//		In2d.y = m_pCart2D[np].y;
//
//		//get cartesian coordinates
//		sphToCart(&In2d, &Out3d);
//		posIn.x = Out3d.x; posIn.y = Out3d.y; posIn.z = Out3d.z;
//		TGeometry *pcCodingGeomtry;
//
//		pcCodingGeomtry->map3DTo2D(&posIn, &posOut);
//
//		posOut.x = (POSType)TGeometry::round(posOut.x);
//		posOut.y = (POSType)TGeometry::round(posOut.y);
//		IPos tmpPos;
//		tmpPos.faceIdx = posOut.faceIdx;
//		tmpPos.u = (Int)(posOut.x);
//		tmpPos.v = (Int)(posOut.y);
//		pcCodingGeomtry->clamp(&tmpPos);
//		pcCodingGeomtry->geoToFramePack(&tmpPos, &m_fpTable[np]);
//	}
//}
//
//Void TEncGOP::xCalculateFBPSNR(TComPicYuv* pcOrgPicYuv, TComPicYuv* pcPicD)
//{
//	Int iBitDepthForPSNRCalc[MAX_NUM_CHANNEL_TYPE];
//	Int iReferenceBitShift[MAX_NUM_CHANNEL_TYPE];
//	Int iOutputBitShift[MAX_NUM_CHANNEL_TYPE];
//	iBitDepthForPSNRCalc[CHANNEL_TYPE_LUMA] = std::max(m_outputBitDepth[CHANNEL_TYPE_LUMA], m_referenceBitDepth[CHANNEL_TYPE_LUMA]);
//	iBitDepthForPSNRCalc[CHANNEL_TYPE_CHROMA] = std::max(m_outputBitDepth[CHANNEL_TYPE_CHROMA], m_referenceBitDepth[CHANNEL_TYPE_CHROMA]);
//	iReferenceBitShift[CHANNEL_TYPE_LUMA] = iBitDepthForPSNRCalc[CHANNEL_TYPE_LUMA] - m_referenceBitDepth[CHANNEL_TYPE_LUMA];
//	iReferenceBitShift[CHANNEL_TYPE_CHROMA] = iBitDepthForPSNRCalc[CHANNEL_TYPE_CHROMA] - m_referenceBitDepth[CHANNEL_TYPE_CHROMA];
//	iOutputBitShift[CHANNEL_TYPE_LUMA] = iBitDepthForPSNRCalc[CHANNEL_TYPE_LUMA] - m_outputBitDepth[CHANNEL_TYPE_LUMA];
//	iOutputBitShift[CHANNEL_TYPE_CHROMA] = iBitDepthForPSNRCalc[CHANNEL_TYPE_CHROMA] - m_outputBitDepth[CHANNEL_TYPE_CHROMA];
//
//	memset(m_dSPSNR, 0, sizeof(Double) * 3);
//	TComPicYuv &picd = *pcPicD;
//	Double SSDspsnr[3] = { 0, 0 ,0 };
//
//	for (Int chan = 0; chan < pcPicD->getNumberValidComponents(); chan++)
//	{
//		const ComponentID ch = ComponentID(chan);
//		const Pel*  pOrg = pcOrgPicYuv->getAddr(ch);
//		const Int   iOrgStride = pcOrgPicYuv->getStride(ch);
//		const Pel*  pRec = picd.getAddr(ch);
//		const Int   iRecStride = picd.getStride(ch);
//
//		for (Int np = 0; np < iNumPoints; np++)
//		{
//			if (!chan)
//			{
//				Int x_loc = (Int)(m_fpTable[np].x);
//				Int y_loc = (Int)(m_fpTable[np].y);
//				Intermediate_Int iDifflp = (pOrg[x_loc + (y_loc*iOrgStride)] << iReferenceBitShift[toChannelType(ch)]) - (pRec[x_loc + (y_loc*iRecStride)] << iOutputBitShift[toChannelType(ch)]);
//				SSDspsnr[chan] += iDifflp * iDifflp;
//			}
//			else
//			{
//				Int x_loc = Int(m_fpTable[np].x >> pcPicD->getComponentScaleX(COMPONENT_Cb));
//				Int y_loc = Int(m_fpTable[np].y >> pcPicD->getComponentScaleY(COMPONENT_Cb));
//				Intermediate_Int iDifflp = (pOrg[x_loc + (y_loc*iOrgStride)] << iReferenceBitShift[toChannelType(ch)]) - (pRec[x_loc + (y_loc*iRecStride)] << iOutputBitShift[toChannelType(ch)]);
//				SSDspsnr[chan] += iDifflp * iDifflp;
//			}
//		}
//	}
//
//	for (Int ch_indx = 0; ch_indx < pcPicD->getNumberValidComponents(); ch_indx++)
//	{
//		const ComponentID ch = ComponentID(ch_indx);
//		const Int maxval = 255 << (iBitDepthForPSNRCalc[toChannelType(ch)] - 8);
//
//		Double fReflpsnr = Double(iNumPoints)*maxval*maxval;
//		m_dSPSNR[ch_indx] = (SSDspsnr[ch_indx] ? 10.0 * log10(fReflpsnr / (Double)SSDspsnr[ch_indx]) : 999.99);
//	}
//}
