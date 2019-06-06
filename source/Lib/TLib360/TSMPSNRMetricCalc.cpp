/**
 * @file 
 * @brief uniformly-RoI
 * @author Gl.M
 * Doc:
 */

#include "TSMPSNRMetricCalc.h"
#include "TAppEncHelper360/TExt360AppEncCfg.h"
#include <sstream>
#include <fstream>
#include <iostream>


TSMPSNRMetric::TSMPSNRMetric()
	: m_bSMPSNREnabled(false)
	, m_pCart3D(NULL)
	, m_fpTable(NULL)
	, m_response(NULL)
{
	m_dSMPSNR[0] = m_dSMPSNR[1] = m_dSMPSNR[2] = 0;
}


TSMPSNRMetric::~TSMPSNRMetric()
{
	if (m_pCart3D)
	{
		free(m_pCart3D); m_pCart3D = NULL;
	}
	if (m_fpTable)
	{
		free(m_fpTable); m_fpTable = NULL;
	}
	if (m_response)
	{
		free(m_response); m_response = NULL;
	}
}




Void TSMPSNRMetric::setOutputBitDepth(Int iOutputBitDepth[MAX_NUM_CHANNEL_TYPE])
{
	for (Int i = 0; i < MAX_NUM_CHANNEL_TYPE; i++)
	{
		m_outputBitDepth[i] = iOutputBitDepth[i];
	}
}

Void TSMPSNRMetric::setReferenceBitDepth(Int iReferenceBitDepth[MAX_NUM_CHANNEL_TYPE])
{
	for (Int i = 0; i < MAX_NUM_CHANNEL_TYPE; i++)
	{
		m_referenceBitDepth[i] = iReferenceBitDepth[i];
	}
}

// get sphere points    TODO
//Void TSMPSNRMetric::sphSampoints(const std::string &cSphDataFile)
//{
//	FILE* fp = fopen(cSphDataFile.c_str(), "r");
//	int x, y;
//	x = y = 0;
//	// read longtitude, latitude
//	if (fscanf(fp, "%d", &m_iSphNumPoints) != 1)
//	{
//		printf("SphData file does not exist.\n");
//		exit(EXIT_FAILURE);
//		fclose(fp);
//	}
//
//	m_pCart3D = (CPos3D*)malloc(sizeof(CPos3D)*(m_iSphNumPoints));
//	memset(m_pCart3D, 0, (sizeof(CPos3D) * m_iSphNumPoints));
//	for (Int z = 0; z < m_iSphNumPoints; z++)
//	{
//		// Reading from latitude,longtitude
//		if (fscanf(fp, "%lf %lf %lf", &m_pCart3D[z].y, &m_pCart3D[z].x, &m_pCart3D[z].z) != 3)
//		{
//			printf("Format error SphData in sphSampoints().\n");
//			exit(EXIT_FAILURE);
//		}
//	}
//	fclose(fp);
//}

Void TSMPSNRMetric::sphSampoints(const std::string &cSphDataFile, const std::string &responseFile)
{
	FILE* fp = fopen(cSphDataFile.c_str(), "r");
	// read longtitude, latitude
	if (fscanf(fp, "%d", &m_iSphNumPoints) != 1)
	{
		printf("SphData file does not exist.\n");
		exit(EXIT_FAILURE);
		fclose(fp);
	}

	m_pCart3D = (CPos3D*)malloc(sizeof(CPos3D)*(m_iSphNumPoints));
	memset(m_pCart3D, 0, (sizeof(CPos3D) * m_iSphNumPoints));
	for (Int z = 0; z < m_iSphNumPoints; z++)
	{
		// Reading from latitude,longtitude
		if (fscanf(fp, "%lf %lf %lf", &m_pCart3D[z].y, &m_pCart3D[z].x, &m_pCart3D[z].z) != 3)
		{
			printf("Format error SphData in sphSampoints().\n");
			exit(EXIT_FAILURE);
		}
	}
	fclose(fp);

	fp = fopen(responseFile.c_str(), "r");
	// read longtitude, latitude
	if (fscanf(fp, "%d", &m_iFeaturePoints) != 1)
	{
		printf("FeatureData file does not exist.\n");
		exit(EXIT_FAILURE);
		fclose(fp);
	}

	m_response = (CPos3D*)malloc(sizeof(CPos3D)*(m_iFeaturePoints));
	memset(m_response, 0, (sizeof(CPos3D) * m_iFeaturePoints));
	for (Int z = 0; z < m_iFeaturePoints; z++)
	{
		// Reading from latitude,longtitude
		if (fscanf(fp, "%lf %lf %lf", &m_response[z].y, &m_response[z].x, &m_response[z].z) != 3)
		{
			printf("Format error SphData in sphSampoints().\n");
			exit(EXIT_FAILURE);
		}
	}
	fclose(fp);

}

void TSMPSNRMetric::sphToCart(CPos2D* sph, CPos3D* out)
{
	POSType fLat = (POSType)(sph->x*S_PI / 180.0);
	POSType fLon = (POSType)(sph->y*S_PI / 180.0);

	out->x = ssin(fLon) * scos(fLat);
	out->y = ssin(fLat);
	out->z = -scos(fLon) * scos(fLat);
}

// create saliencymap table   TODO
void TSMPSNRMetric::createTable(TGeometry *pcCodingGeomtry)
{
	Int iNumPoints = m_iSphNumPoints;
	CPos2D In2d;
	CPos3D Out3d;
	SPos posIn, posOut;
	m_fpTable = (IPos2D*)malloc(iNumPoints * sizeof(IPos2D));
	m_ffTable = (IPos2D*)malloc(m_iFeaturePoints * sizeof(IPos2D));

	for (Int np = 0; np < iNumPoints; np++)
	{
		In2d.x = m_pCart3D[np].x;
		In2d.y = m_pCart3D[np].y;

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
	for (Int np = 0; np < m_iFeaturePoints; np++)
	{
		In2d.x = m_response[np].x;
		In2d.y = m_response[np].y;

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
		pcCodingGeomtry->geoToFramePack(&tmpPos, &m_ffTable[np]);
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
Void TSMPSNRMetric::xCalculateSMPSNR(TComPicYuv* pcOrgPicYuv, TComPicYuv* pcPicD)
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

	memset(m_dSMPSNR, 0, sizeof(Double) * 3);
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
		m_dSMPSNR[ch_indx] = (SSDspsnr[ch_indx] ? 10.0 * log10(fReflpsnr / (Double)SSDspsnr[ch_indx]) : 999.99);
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
