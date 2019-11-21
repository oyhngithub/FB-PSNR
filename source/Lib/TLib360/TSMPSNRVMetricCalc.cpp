/**
 * @file 
 * @brief uniformly-RoI
 * @author Gl.M
 * Doc:
 */

#include "TSMPSNRVMetricCalc.h"
#include "TAppEncHelper360/TExt360AppEncCfg.h"
#include <sstream>
#include <fstream>
#include <iostream>


TSMPSNRVMetric::TSMPSNRVMetric()
	: m_bSMPSNRVEnabled(false)
	, m_pCart3D(NULL)
	, m_fpTable(NULL)
	, m_response(NULL)
{
	m_dSMPSNRV[0] = m_dSMPSNRV[1] = m_dSMPSNRV[2] = 0;
}


TSMPSNRVMetric::~TSMPSNRVMetric()
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




Void TSMPSNRVMetric::setOutputBitDepth(Int iOutputBitDepth[MAX_NUM_CHANNEL_TYPE])
{
	for (Int i = 0; i < MAX_NUM_CHANNEL_TYPE; i++)
	{
		m_outputBitDepth[i] = iOutputBitDepth[i];
	}
}

Void TSMPSNRVMetric::setReferenceBitDepth(Int iReferenceBitDepth[MAX_NUM_CHANNEL_TYPE])
{
	for (Int i = 0; i < MAX_NUM_CHANNEL_TYPE; i++)
	{
		m_referenceBitDepth[i] = iReferenceBitDepth[i];
	}
}

// get sphere points    TODO
/*Void TSMPSNRVMetric::sphSampoints(const std::string &cSphDataFile)
{
	FILE* fp = fopen(cSphDataFile.c_str(), "r");
	int x, y;
	x = y = 0;
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
}

void TSMPSNRVMetric::sphToCart(CPos2D* sph, CPos3D* out)
{
	POSType fLat = (POSType)(sph->x*S_PI / 180.0);
	POSType fLon = (POSType)(sph->y*S_PI / 180.0);

	out->x = ssin(fLon) * scos(fLat);
	out->y = ssin(fLat);
	out->z = -scos(fLon) * scos(fLat);
}

// create saliencymap table   TODO
void TSMPSNRVMetric::createTable(TGeometry *pcCodingGeomtry)
{
	Int iNumPoints = m_iSphNumPoints;
	CPos2D In2d;
	CPos3D Out3d;
	SPos posIn, posOut;
	m_fpTable = (IPos2D*)malloc(iNumPoints * sizeof(IPos2D));

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
	}*/
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
//}

// pass the pics before and after the encode
// m_fpTable stores all the data that already transform into different projections
Void TSMPSNRVMetric::xCalculateSMPSNRV(TComPicYuv* pcOrgPicYuv, TComPicYuv* pcPicD)
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

	memset(m_dSMPSNRV, 0, sizeof(Double) * 3);
	TComPicYuv &picd = *pcPicD;

	for (Int chan = 0; chan < pcPicD->getNumberValidComponents(); chan++)
	{
		Double trueWeight = 1;
		const ComponentID ch = ComponentID(chan);
		const Pel*  pOrg = pcOrgPicYuv->getAddr(ch);
		const Int   iOrgStride = pcOrgPicYuv->getStride(ch);
		const Pel*  pRec = picd.getAddr(ch);
		const Int   iRecStride = picd.getStride(ch);

		const Int   iWidth = pcPicD->getWidth(ch);
		const Int   iHeight = pcPicD->getHeight(ch);
		Double fWeightSum = 0;
		Double fWeight = 1;

		//Int   iSize   = iWidth*iHeight;

		//Double SSDwpsnr = 0;

		Double SSDwpsnr = 0;
		for (Int np = 0; np < iNumPoints; np++)
		{
			Int x_loc;
			Int y_loc;
			if (!chan)
			{
				x_loc = (Int)(m_fpTable[np].x);
				y_loc = (Int)(m_fpTable[np].y);
			}
			else
			{
				x_loc = Int(m_fpTable[np].x >> pcPicD->getComponentScaleX(COMPONENT_Cb));
				y_loc = Int(m_fpTable[np].y >> pcPicD->getComponentScaleY(COMPONENT_Cb));
			}

			int x, y;
			x = x_loc;
			y = y_loc;

			Double fWeightResponse = m_pCart3D[np].z;

			// read W1
			if (fWeightResponse > 0)
				trueWeight = 0.5 * fWeightResponse;

				if (m_codingGeoType == SVIDEO_EQUIRECT || m_codingGeoType == SVIDEO_NEWUNIFORMMAP)
				{
					if (!chan)
					{
						fWeight = m_fErpWeight_Y[y];
					}
					else
					{
						fWeight = m_fErpWeight_C[y];
					}
				}
				//if (((pOrg + y * iOrgStride)[x]) > 255 || ((pOrg + y * iOrgStride)[x]) < 0 || ((pRec + y * iRecStride)[x]) > 255 || ((pRec + y * iRecStride)[x]) < 0)
				//	printf("x:%d ,y:%d\norg:%d, rec:%d\n", x, y, (pOrg + y * iOrgStride)[x], (pRec + y * iRecStride)[x]);
				//int test = iReferenceBitShift[toChannelType(ch)];
				//assert(y < iHeight && x < iWidth);

				int org = ((pOrg + y * iOrgStride)[x]);
				int rec = ((pRec + y * iRecStride)[x]);
				//int refShift = iReferenceBitShift[toChannelType(ch)];
				//int outputShift = iOutputBitShift[toChannelType(ch)];
				Intermediate_Int iDiff = (Intermediate_Int)((org << iReferenceBitShift[toChannelType(ch)]) - (rec << iOutputBitShift[toChannelType(ch)]));
				//Intermediate_Int iDiff = (Intermediate_Int)((pOrg[x] << iReferenceBitShift[toChannelType(ch)]) - (pRec[x] << iOutputBitShift[toChannelType(ch)]));

				if ((m_codingGeoType == SVIDEO_CUBEMAP)
#if SVIDEO_ADJUSTED_CUBEMAP
					|| (m_codingGeoType == SVIDEO_ADJUSTEDCUBEMAP)
#endif
#if SVIDEO_EQUATORIAL_CYLINDRICAL
					|| (m_codingGeoType == SVIDEO_EQUATORIALCYLINDRICAL)
#endif
#if SVIDEO_EQUIANGULAR_CUBEMAP
					|| (m_codingGeoType == SVIDEO_EQUIANGULARCUBEMAP)
#endif
					)
				{
					if (!chan)
					{
						if (iWidth / 4 == iHeight / 3 && x >= iWidth / 4 && (y < iHeight / 3 || y >= 2 * iHeight / 3))
						{
							fWeight = 0;
						}
						else
						{
							fWeight = m_fCubeWeight_Y[(m_iCodingFaceWidth)*(y % (m_iCodingFaceHeight)) + (x % (m_iCodingFaceWidth))];
						}

					}
					else
					{
						if (iWidth / 4 == iHeight / 3 && x >= iWidth / 4 && (y < iHeight / 3 || y >= 2 * iHeight / 3))
						{
							fWeight = 0;
						}
						else
						{
							fWeight = m_fCubeWeight_C[(m_iCodingFaceWidth >> (pcPicD->getComponentScaleX(COMPONENT_Cb)))*(y % (m_iCodingFaceHeight >> (pcPicD->getComponentScaleY(COMPONENT_Cb)))) + (x % (m_iCodingFaceWidth >> (pcPicD->getComponentScaleX(COMPONENT_Cb))))];
						}
					}
				}
#if SVIDEO_ADJUSTED_EQUALAREA
				else if (m_codingGeoType == SVIDEO_ADJUSTEDEQUALAREA)
#else
				else if (m_codingGeoType == SVIDEO_EQUALAREA)
#endif
				{
					if (!chan)
					{
						fWeight = m_fEapWeight_Y[y*iWidth + x];
					}
					else
					{
						fWeight = m_fEapWeight_C[y*iWidth + x];
					}
				}
				else if (m_codingGeoType == SVIDEO_OCTAHEDRON)
				{
					if (!chan)
					{
						fWeight = m_fOctaWeight_Y[iWidth*y + x];
					}
					else
					{
						fWeight = m_fOctaWeight_C[y*iWidth + x];
					}
				}
				else if (m_codingGeoType == SVIDEO_ICOSAHEDRON)
				{
					if (!chan)
					{
						fWeight = m_fIcoWeight_Y[iWidth*y + x];
					}
					else
					{
						fWeight = m_fIcoWeight_C[y*iWidth + x];
					}
				}
#if SVIDEO_WSPSNR_SSP
				else if (m_codingGeoType == SVIDEO_SEGMENTEDSPHERE)
				{
					if (!chan)
					{
						fWeight = m_fSspWeight_Y[iWidth*y + x];
					}
					else
					{
						fWeight = m_fSspWeight_C[y*iWidth + x];
					}
				}
#endif
#if SVIDEO_ROTATED_SPHERE
				else if (m_codingGeoType == SVIDEO_ROTATEDSPHERE)
				{
					if (!chan)
					{
						fWeight = m_fRspWeight_Y[iWidth*y + x];
					}
					else
					{
						fWeight = m_fRspWeight_C[y*iWidth + x];
					}
				}
#endif
#if SVIDEO_ERP_PADDING
				else if (m_codingGeoType == SVIDEO_EQUIRECT && m_bPERP)
				{
#if 1
					ChromaFormat fmt = pcPicD->getChromaFormat();
					if ((x < (SVIDEO_ERP_PAD_L >> getComponentScaleX(ch, fmt))) || (x >= (iWidth - (SVIDEO_ERP_PAD_R >> getComponentScaleX(ch, fmt)))))
						fWeight = 0;
					else
						fWeight = (!chan) ? m_fErpWeight_Y[y] : m_fErpWeight_C[y];
#else
					if (!chan)
					{
						if ((x < SVIDEO_ERP_PAD_L) || (x >= (iWidth - SVIDEO_ERP_PAD_R)))
							fWeight = 0;
						else
							fWeight = m_fErpWeight_Y[y];
					}
					else
					{
						ComponentID chId = ComponentID(chan);
						ChromaFormat fmt = pcPicD->getChromaFormat();

						if ((x < (SVIDEO_ERP_PAD_L >> getComponentScaleX(chId, fmt))) || (x >= (iWidth - (SVIDEO_ERP_PAD_R >> getComponentScaleX(chId, fmt)))))
							fWeight = 0;
						else
							fWeight = m_fErpWeight_C[y];
					}
#endif
				}
#endif	
				if (fWeight > 0) {
					fWeight = (fWeight - m_min) / (m_max - m_min);
					trueWeight += 0.5 * fWeight;
				}
				// read w2	
				fWeightSum += trueWeight;
				assert(iDiff <= 255);
				SSDwpsnr += iDiff * iDiff*trueWeight;
		}
		const Int maxval = 255 << (iBitDepthForPSNRCalc[toChannelType(ch)] - 8);
		m_dSMPSNRV[ch] = (SSDwpsnr ? 10.0 * log10((maxval * maxval*fWeightSum) / (Double)SSDwpsnr) : 999.99);
	}
	// calculate w2
	//for (Int chan = 0; chan < pcPicD->getNumberValidComponents(); chan++)
	//{
	//	const ComponentID ch = ComponentID(chan);
	//	const Pel*  pOrg = pcOrgPicYuv->getAddr(ch);
	//	const Int   iOrgStride = pcOrgPicYuv->getStride(ch);
	//	const Pel*  pRec = picd.getAddr(ch);
	//	const Int   iRecStride = picd.getStride(ch);

	//	const Int   iWidth = pcPicD->getWidth(ch);
	//	const Int   iHeight = pcPicD->getHeight(ch);

	//	for (Int np = 0; np < m_iFeaturePoints; np++)
	//	{
	//		Int x_loc;
	//		Int y_loc;
	//		if (!chan)
	//		{
	//			x_loc = (Int)(m_ffTable[np].x);
	//			y_loc = (Int)(m_ffTable[np].y);
	//		}
	//		else
	//		{
	//			x_loc = Int(m_ffTable[np].x >> pcPicD->getComponentScaleX(COMPONENT_Cb));
	//			y_loc = Int(m_ffTable[np].y >> pcPicD->getComponentScaleY(COMPONENT_Cb));
	//		}

	//		int x, y;
	//		x = x_loc;
	//		y = y_loc;

	//		
	//		if (((pOrg + y * iOrgStride)[x]) > 255 || ((pOrg + y * iOrgStride)[x]) < 0 || ((pRec + y * iRecStride)[x]) > 255 || ((pRec + y * iRecStride)[x]) < 0)
	//			printf("x:%d ,y:%d\norg:%d, rec:%d\n", x, y, (pOrg + y * iOrgStride)[x], (pRec + y * iRecStride)[x]);
	//		//int test = iReferenceBitShift[toChannelType(ch)];
	//		int org = ((pOrg + y * iOrgStride)[x]);
	//		int rec = ((pRec + y * iRecStride)[x]);
	//		//int refShift = iReferenceBitShift[toChannelType(ch)];
	//		//int outputShift = iOutputBitShift[toChannelType(ch)];
	//		Intermediate_Int iDiff = (Intermediate_Int)((((pOrg + y * iOrgStride)[x]) << iReferenceBitShift[toChannelType(ch)]) - (((pRec + y * iRecStride)[x]) << iOutputBitShift[toChannelType(ch)]));
	//		//Intermediate_Int iDiff = (Intermediate_Int)((pOrg[x] << iReferenceBitShift[toChannelType(ch)]) - (pRec[x] << iOutputBitShift[toChannelType(ch)]));
	//		if (m_response[np].z > 0)
	//			w2Weight[chan] += m_response[np].z;
	//		assert(iDiff <= 255);
	//		w2[chan] += iDiff * iDiff*m_response[np].z;
	//	}
	//}

	//// calculate 3 path individualy
	//for (Int ch_indx = 0; ch_indx < pcPicD->getNumberValidComponents(); ch_indx++)
	//{
	//	const ComponentID ch = ComponentID(ch_indx);
	//	const Int maxval = 255 << (iBitDepthForPSNRCalc[toChannelType(ch)] - 8);

	//	Double fReflpsnr = Double(iNumPoints)*maxval*maxval;
	//	w2PSNR[ch_indx] = (w2[ch_indx] ? 10.0 * log10(w2Weight[ch_indx] * fReflpsnr / (Double)w2[ch_indx]) : 999.99);
	//}

	/*for (int i = 0; i < 3; ++i) {
		std::cout << w1PSNR[i] << ", " << w2PSNR[i] << ", " << w3PSNR[i];
	}*/
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

#if SVIDEO_E2E_METRICS
Void TSMPSNRVMetric::xCalculateE2ESMPSNRV(TComPicYuv* pcPicYuv, TComPicYuv* pcOrigPicYuv)
#else
Void TWSPSNRMetric::xCalculateE2EWSPSNR(TComPicYuv* pcPicYuv, Int iPOC)
#endif
{
#if SVIDEO_E2E_METRICS
	xCalculateSMPSNRV(pcOrigPicYuv, pcPicYuv);
#else
	xCalculateWSPSNR(m_pcOrgPicYuv, m_pcRecPicYuv);
#endif
}
