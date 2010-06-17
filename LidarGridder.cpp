/***********************************************************************
LidarGridder - Experimental program to resample LiDAR data onto a
regular grid.
Copyright (c) 2009-2010 Oliver Kreylos

This file is part of the LiDAR processing and analysis package.

The LiDAR processing and analysis package is free software; you can
redistribute it and/or modify it under the terms of the GNU General
Public License as published by the Free Software Foundation; either
version 2 of the License, or (at your option) any later version.

The LiDAR processing and analysis package is distributed in the hope
that it will be useful, but WITHOUT ANY WARRANTY; without even the
implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with the LiDAR processing and analysis package; if not, write to the
Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
02111-1307 USA
***********************************************************************/

#include <stdlib.h>
#include <string.h>
#include <string>
#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <Misc/Array.h>
#include <Misc/ThrowStdErr.h>
#include <Misc/File.h>
#include <Misc/Timer.h>
#include <Math/Math.h>
#include <Math/Constants.h>
#include <Geometry/Box.h>
#include <GL/gl.h>
#include <GL/GLColor.h>
#include <GL/GLVertexArrayParts.h>
#include <GL/GLMaterial.h>
#include <GL/Extensions/GLARBVertexBufferObject.h>
#include <GL/GLObject.h>
#include <GL/GLContextData.h>
#include <GL/GLGeometryWrappers.h>
#include <GL/GLGeometryVertex.h>
#include <Vrui/Vrui.h>
#include <Vrui/Application.h>

#include "LidarTypes.h"
#include "LidarProcessOctree.h"

class PointNormalSampler // Class to sample the implicit surface defined by a LiDAR point cloud and its normal vector at an (x, y) position
	{
	/* Elements: */
	private:
	double pos[2]; // Sample point's (x, y) position
	double gridCellSize[2]; // Cell size of target grid in x and y
	double numLobes; // Number of lobes in Lanczos reconstruction filter
	double accumulator; // Accumulated convolution between point cloud and filter
	double weightSum; // Sum of weights for all accumulated points
	double absWeightSum; // Sum of absolute weights for all accumulated points
	double diffSum[2]; // Accumulated partial derivatives of elevation function
	double diffWeightSum[2]; // Sum of weights for accumulated partial derivatives
	
	/* Constructors and destructors: */
	public:
	PointNormalSampler(const double sPos[2],const double sGridCellSize[2],int sNumLobes) // Creates an empty sampler
		:numLobes(double(sNumLobes)),
		 accumulator(0.0),weightSum(0.0),absWeightSum(0.0)
		{
		for(int i=0;i<2;++i)
			{
			pos[i]=sPos[i];
			gridCellSize[i]=sGridCellSize[i];
			}
		for(int i=0;i<2;++i)
			{
			diffSum[i]=0.0;
			diffWeightSum[i]=0.0;
			}
		}
	
	/* Methods: */
	void dilate(void)
		{
		for(int i=0;i<2;++i)
			gridCellSize[i]*=2.0;
		accumulator=0.0;
		weightSum=0.0;
		absWeightSum=0.0;
		for(int i=0;i<2;++i)
			{
			diffSum[i]=0.0;
			diffWeightSum[i]=0.0;
			}
		}
	Box calcSampleBox(void) const
		{
		/* Calculate the sampler's bounding box: */
		Box sampleBox;
		for(int i=0;i<2;++i)
			{
			sampleBox.min[i]=Scalar(pos[i]-gridCellSize[i]*numLobes);
			sampleBox.max[i]=Scalar(pos[i]+gridCellSize[i]*numLobes);
			}
		sampleBox.min[2]=Box::full.min[2];
		sampleBox.max[2]=Box::full.max[2];
		
		return sampleBox;
		}
	void operator()(const LidarPoint& p)
		{
		/* Calculate the Lanczos filter weights and its partial derivatives for the LiDAR point: */
		double fx=Math::Constants<double>::pi/gridCellSize[0];
		double sx=(double(p[0])-pos[0])*fx;
		double ssx=Math::sin(sx);
		double sxn=sx/numLobes;
		double ssxn=Math::sin(sxn);
		double lancsx=sx!=0.0?(ssx*ssxn)/(sx*sxn):1.0;
		
		double fy=Math::Constants<double>::pi/gridCellSize[1];
		double sy=(double(p[1])-pos[1])*fy;
		double ssy=Math::sin(sy);
		double syn=sy/numLobes;
		double ssyn=Math::sin(syn);
		double lancsy=sy!=0.0?(ssy*ssyn)/(sy*syn):1.0;
		
		double pointWeight=lancsx*lancsy;
		double diffWeight[2];
		
		diffWeight[0]=sx!=0.0?lancsy*((numLobes*Math::cos(sx)*ssxn+Math::cos(sxn)*ssx)/sx-2.0*lancsx)/fx/sx:0.0;
		diffWeight[1]=sy!=0.0?lancsx*((numLobes*Math::cos(sy)*ssyn+Math::cos(syn)*ssy)/sy-2.0*lancsy)/fy/sy:0.0;
		
		/* Accumulate the weighted point and partial derivatives: */
		double z=double(p[2]);
		accumulator+=z*pointWeight;
		weightSum+=pointWeight;
		absWeightSum+=Math::abs(pointWeight);
		for(int i=0;i<2;++i)
			{
			diffSum[i]+=z*diffWeight[i];
			diffWeightSum[i]+=diffWeight[i];
			}
		}
	double getWeightSum(void) const // Returns the sum of point weights
		{
		return weightSum;
		}
	double getAbsWeightSum(void) const // Returns the sum of absolute point weights
		{
		return absWeightSum;
		}
	double getValue(void) const // Returns the convolution result
		{
		/* Return the weighted average: */
		return accumulator/weightSum;
		}
	double getDiff(int direction) const // Returns one of the partial derivatives
		{
		/* Return the weighted average of the partial derivative: */
		return (diffSum[direction]*weightSum-accumulator*diffWeightSum[direction])/Math::sqr(weightSum);
		}
	};

class DifferenceCalculator
	{
	/* Elements: */
	private:
	double gridOrigin[2]; // Grid's origin position
	int gridSize[2]; // Grid size in number of vertices
	double gridCellSize[2]; // Cell size of target grid in x and y
	const float* grid; // Array of grid elevations
	float* diffs; // Array of differences for each grid cell
	int* diffWeights;
	
	/* Constructors and destructors: */
	public:
	DifferenceCalculator(const double sGridOrigin[2],const int sGridSize[2],const double sGridCellSize[2],const float* sGrid)
		:grid(sGrid),
		 diffs(0),diffWeights(0)
		{
		for(int i=0;i<2;++i)
			{
			gridOrigin[i]=sGridOrigin[i];
			gridSize[i]=sGridSize[i];
			gridCellSize[i]=sGridCellSize[i];
			}
		diffs=new float[gridSize[1]*gridSize[0]];
		diffWeights=new int[gridSize[1]*gridSize[0]];
		for(int y=0;y<gridSize[1];++y)
			for(int x=0;x<gridSize[0];++x)
				{
				diffs[y*gridSize[0]+x]=0.0f;
				diffWeights[y*gridSize[0]+x]=0;
				}
		}
	~DifferenceCalculator(void)
		{
		delete[] diffs;
		delete[] diffWeights;
		}
	
	/* Methods: */
	void operator()(const LidarPoint& p)
		{
		/* Find the grid cell containing the point: */
		int pi[2];
		Scalar pd[2];
		for(int i=0;i<2;++i)
			{
			pd[i]=p[i]-gridOrigin[i];
			pi[i]=int(Math::floor(pd[i]));
			if(pi[i]<0||pi[i]>=gridSize[i]-1) // Bail out if the point is outside the grid's interior
				return;
			pd[i]-=Scalar(pi[i]);
			}
		
		/* Interpolate the grid elevation: */
		const float* gb=&grid[pi[1]*gridSize[0]+pi[0]];
		float h01=gb[0]*(Scalar(1)-pd[0])+gb[1]*pd[0];
		float h23=gb[gridSize[0]+0]*(Scalar(1)-pd[0])+gb[gridSize[0]+1]*pd[0];
		float h=h01*(Scalar(1)-pd[1])+h23*pd[1];
		
		/* Accumulate the difference: */
		diffs[pi[1]*gridSize[0]+pi[0]]+=Math::sqr(float(p[2])-h);
		++diffWeights[pi[1]*gridSize[0]+pi[0]];
		}
	const float* getDiffs(void) const
		{
		return diffs;
		}
	const int* getDiffWeights(void) const
		{
		return diffWeights;
		}
	float getDiffRMS(int x,int y) const
		{
		return Math::sqrt(diffs[y*gridSize[0]+x]/float(diffWeights[y*gridSize[0]+x]));
		}
	};

class LidarGridder:public Vrui::Application,public GLObject
	{
	/* Embedded classes: */
	private:
	typedef struct Geometry::Point<float,3> Point;
	typedef GLGeometry::Vertex<GLfloat,2,void,0,GLfloat,GLfloat,3> DemVertex;
	
	struct DataItem:public GLObject::DataItem
		{
		/* Elements: */
		public:
		GLuint demBufferIds[2]; // Buffer IDs for vertex and index buffer
		
		/* Constructors and destructors: */
		public:
		DataItem(void)
			{
			demBufferIds[0]=demBufferIds[1]=0;
			if(GLARBVertexBufferObject::isSupported())
				{
				GLARBVertexBufferObject::initExtension();
				glGenBuffersARB(2,demBufferIds);
				}
			}
		virtual ~DataItem(void)
			{
			if(demBufferIds[0]!=0)
				glDeleteBuffersARB(2,demBufferIds);
			}
		};
	
	/* Elements: */
	LidarProcessOctree* lpo; // The out-of-core LiDAR octree
	Cube lpoDomain; // Domain of the LiDAR octree
	Geometry::Box<double,2> gridBox; // Box limiting the extent of the extracted DEM
	double gridCellSize[2]; // Cell size for the extracted DEM
	int numLobes; // Number of lobes in the Lanczos reconstruction filter
	Misc::Array<DemVertex,2> dem; // The current DEM
	
	/* Private methods: */
	void createDem(void);
	
	/* Constructors and destructors: */
	public:
	LidarGridder(int& argc,char**& argv,char**& appDefaults);
	virtual ~LidarGridder(void);
	
	/* Methods from Vrui::Application: */
	virtual void display(GLContextData& contextData) const;
	
	/* Methods from GLObject: */
	virtual void initContext(GLContextData& contextData) const;
	};

/*****************************
Methods of class LidarGridder:
*****************************/

void LidarGridder::createDem(void)
	{
	/* Calculate the exact size of the DEM grid: */
	double gridOrigin[2];
	int gridSize[2];
	for(int i=0;i<2;++i)
		{
		gridSize[i]=int(Math::ceil((gridBox.max[i]-gridBox.min[i])/gridCellSize[i]))+1;
		double overshoot=double(gridSize[i]-1)*gridCellSize[i]-(gridBox.max[i]-gridBox.min[i]);
		gridOrigin[i]=gridBox.min[i]-0.5*overshoot;
		}
	
	std::ios::fmtflags flags=std::cout.flags();
	std::cout.setf(std::ios::fixed);
	std::streamsize precision=std::cout.precision(3);
	std::cout<<"Generating DEM:"<<std::endl;
	std::cout<<"Grid cell size: "<<gridCellSize[0]<<" x "<<gridCellSize[1]<<std::endl;
	std::cout<<"Grid size: "<<gridSize[0]<<" x "<<gridSize[1]<<std::endl;
	std::cout<<"Lower-left corner : "<<gridOrigin[0]<<", "<<gridOrigin[1]<<std::endl;
	std::cout<<"Upper-right corner: "<<gridOrigin[0]+double(gridSize[0]-1)*gridCellSize[0]<<", "<<gridOrigin[1]+double(gridSize[1]-1)*gridCellSize[1]<<std::endl;
	std::cout.precision(precision);
	std::cout.flags(flags);
	
	/* Create the DEM structure: */
	dem.resize(gridSize);
	
	/* Sample the LiDAR data set: */
	Misc::Timer timer;
	Vector offset=Point::origin-lpoDomain.getCenter();
	for(int y=0;y<gridSize[1];++y)
		{
		for(int x=0;x<gridSize[0];++x)
			{
			/* Calculate the grid vertex position: */
			double pos[2];
			pos[0]=gridOrigin[0]+double(x)*gridCellSize[0];
			pos[1]=gridOrigin[1]+double(y)*gridCellSize[1];
			
			/* Create a sampler: */
			PointNormalSampler sampler(pos,gridCellSize,numLobes);
			
			/* Repeatedly dilate the sampler until there is enough confidence in the returned value: */
			double yTotal=0.0;
			double wTotal=1.0;
			int numDilations=0;
			double confidence;
			while(true)
				{
				/* Sample the LiDAR data set: */
				lpo->processPointsInBox(sampler.calcSampleBox(),sampler);
				
				/* Calculate the sampling confidence: */
				confidence=sampler.getWeightSum();
				if(confidence>=1.0)
					{
					/* Store the step's contribution and stop: */
					yTotal+=wTotal*sampler.getValue();
					break;
					}
				
				/* Calculate the pre-dilation weight for this step: */
				double pdw=(confidence-0.25)/(1.0-0.25);
				if(pdw>0.0)
					{
					/* Store this step's partial contribution: */
					yTotal+=wTotal*pdw*sampler.getValue();
					wTotal*=(1.0-pdw);
					}
				
				/* Dilate the sampler and try again: */
				sampler.dilate();
				++numDilations;
				}
			
			/* Create the DEM vertex: */
			dem(x,y).texCoord=DemVertex::TexCoord((double(numDilations)+0.5)/8,confidence);
			DemVertex::Normal n(-sampler.getDiff(0),-sampler.getDiff(1),1);
			n.normalize();
			dem(x,y).normal=n;
			dem(x,y).position=DemVertex::Position(pos[0]+offset[0],pos[1]+offset[1],yTotal+offset[2]);
			}
		std::cout<<"\b\b\b\b"<<std::setw(3)<<(100*(y+1)+gridSize[1]/2)/gridSize[1]<<"%"<<std::flush;
		}
	std::cout<<std::endl;
	
	#if 1
	/* Calculate normal vectors: */
	for(int y=0;y<dem.getSize(1);++y)
		{
		dem(0,y).normal[0]=-(dem(1,y).position[2]-dem(0,y).position[2])/gridCellSize[0];
		for(int x=1;x<dem.getSize(0)-1;++x)
			dem(x,y).normal[0]=-(dem(x+1,y).position[2]-dem(x-1,y).position[2])/(2.0*gridCellSize[0]);
		dem(dem.getSize(0)-1,y).normal[0]=-(dem(dem.getSize(0)-1,y).position[2]-dem(dem.getSize(0)-2,y).position[2])/gridCellSize[0];
		}
	for(int x=0;x<dem.getSize(0);++x)
		{
		dem(x,0).normal[1]=-(dem(x,1).position[2]-dem(x,0).position[2])/gridCellSize[1];
		for(int y=1;y<dem.getSize(1)-1;++y)
			dem(x,y).normal[1]=-(dem(x,y+1).position[2]-dem(x,y-1).position[2])/(2.0*gridCellSize[1]);
		dem(x,dem.getSize(1)-1).normal[1]=-(dem(x,dem.getSize(1)-1).position[2]-dem(x,dem.getSize(1)-2).position[2])/gridCellSize[1];
		}
	for(int x=0;x<dem.getSize(0);++x)
		for(int y=0;y<dem.getSize(1);++y)
			{
			dem(x,y).normal[2]=1.0;
			dem(x,y).normal.normalize();
			}
	#endif
	timer.elapse();
	std::cout<<"Grid generation time: "<<timer.getTime()*1000.0<<" ms"<<std::endl;
	}

LidarGridder::LidarGridder(int& argc,char**& argv,char**& appDefaults)
	:Vrui::Application(argc,argv,appDefaults),
	 lpo(0)
	{
	/* Parse the command line: */
	const char* lidarFileName=0;
	size_t cacheSize=512;
	gridBox=Geometry::Box<double,2>::full;
	gridCellSize[0]=gridCellSize[1]=1.0;
	numLobes=3;
	const char* gridFileName=0;
	for(int i=1;i<argc;++i)
		{
		if(argv[i][0]=='-')
			{
			if(strcasecmp(argv[i]+1,"cache")==0)
				{
				++i;
				cacheSize=atoi(argv[i]);
				}
			else if(strcasecmp(argv[i]+1,"gridMin")==0)
				{
				for(int j=0;j<2;++j)
					{
					++i;
					gridBox.min[j]=atof(argv[i]);
					}
				}
			else if(strcasecmp(argv[i]+1,"gridMax")==0)
				{
				for(int j=0;j<2;++j)
					{
					++i;
					gridBox.max[j]=atof(argv[i]);
					}
				}
			else if(strcasecmp(argv[i]+1,"gridCellSize")==0)
				{
				for(int j=0;j<2;++j)
					{
					++i;
					gridCellSize[j]=atof(argv[i]);
					}
				}
			else if(strcasecmp(argv[i]+1,"numLobes")==0)
				{
				++i;
				numLobes=atoi(argv[i]);
				}
			else
				std::cerr<<"Unrecognized switch "<<argv[i]<<std::endl;
			}
		else if(lidarFileName==0)
			lidarFileName=argv[i];
		else if(gridFileName==0)
			gridFileName=argv[i];
		else
			std::cerr<<"Unrecognized argument "<<argv[i]<<std::endl;
		}
	if(lidarFileName==0)
		Misc::throwStdErr("No LiDAR file name provided");
	
	/* Open the LiDAR data set: */
	lpo=new LidarProcessOctree(lidarFileName,cacheSize*1024*1024);
	lpoDomain=lpo->getDomain();
	
	/* Limit the grid box to the LiDAR data set's domain: */
	for(int i=0;i<2;++i)
		{
		if(gridBox.min[i]<lpoDomain.getMin()[i])
			gridBox.min[i]=lpoDomain.getMin()[i];
		if(gridBox.max[i]>lpoDomain.getMax()[i])
			gridBox.max[i]=lpoDomain.getMax()[i];
		}
	
	/* Sample the initial DEM: */
	createDem();
	
	if(gridFileName!=0)
		{
		/* Save the dem as a grid file: */
		Misc::File gridFile(gridFileName,"wb",Misc::File::LittleEndian);
		gridFile.write<int>(dem.getSize(0));
		gridFile.write<int>(dem.getSize(1));
		gridFile.write<float>(dem(0,0).position[0]);
		gridFile.write<float>(dem(0,0).position[1]);
		gridFile.write<float>(dem(dem.getSize(0)-1,dem.getSize(1)-1).position[0]);
		gridFile.write<float>(dem(dem.getSize(0)-1,dem.getSize(1)-1).position[1]);
		for(int y=0;y<dem.getSize(1);++y)
			for(int x=0;x<dem.getSize(0);++x)
				gridFile.write<float>(dem(x,y).position[2]);
		}
	
	/* Initialize the navigation transformation: */
	Vrui::Point center(Math::mid(gridBox.min[0],gridBox.max[0]),Math::mid(gridBox.min[1],gridBox.max[1]),0);
	Vrui::Scalar size=Math::div2(Math::sqrt(Math::sqr(gridBox.max[0]-gridBox.min[0])+Math::sqr(gridBox.max[1]-gridBox.min[1])));
	Vrui::setNavigationTransformation(Vrui::Point::origin,size,Vrui::Vector(0,1,0));
	}

LidarGridder::~LidarGridder(void)
	{
	delete lpo;
	}

void LidarGridder::display(GLContextData& contextData) const
	{
	/* Save and set up OpenGL state: */
	glPushAttrib(GL_ENABLE_BIT|GL_LIGHTING_BIT);
	#if 0
	glDisable(GL_CULL_FACE);
	glLightModeli(GL_LIGHT_MODEL_TWO_SIDE,GL_TRUE);
	#endif
	glDisable(GL_COLOR_MATERIAL);
	glMaterial(GLMaterialEnums::FRONT_AND_BACK,GLMaterial(GLMaterial::Color(0.75f,0.75f,0.75f),GLMaterial::Color(0.5f,0.5f,0.5f),25.0f));
	
	glEnable(GL_TEXTURE_1D);
	glTexEnvi(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,GL_MODULATE);
	glLightModeli(GL_LIGHT_MODEL_COLOR_CONTROL,GL_SEPARATE_SPECULAR_COLOR);
	
	#if 0
	GLColor<GLfloat,3> spikinessColorMap[4]=
		{
		GLColor<GLfloat,3>(1.0f,0.0f,0.0f),GLColor<GLfloat,3>(1.0f,1.0f,0.0f),GLColor<GLfloat,3>(0.0f,1.0f,0.0f),
		GLColor<GLfloat,3>(0.0f,1.0f,0.0f)
		};
	glTexParameteri(GL_TEXTURE_1D,GL_TEXTURE_MIN_FILTER,GL_LINEAR);
	glTexParameteri(GL_TEXTURE_1D,GL_TEXTURE_MAG_FILTER,GL_LINEAR);
	glTexParameteri(GL_TEXTURE_1D,GL_TEXTURE_WRAP_S,GL_CLAMP);
	glTexImage1D(GL_TEXTURE_1D,0,GL_RGB,4,0,GL_RGB,GL_FLOAT,spikinessColorMap);
	#else
	GLColor<GLfloat,3> dilationColorMap[8]=
		{
		GLColor<GLfloat,3>(1.0f,0.0f,0.0f),GLColor<GLfloat,3>(1.0f,1.0f,0.0f),GLColor<GLfloat,3>(0.0f,1.0f,0.0f),
		GLColor<GLfloat,3>(0.0f,1.0f,1.0f),GLColor<GLfloat,3>(0.0f,0.0f,1.0f),GLColor<GLfloat,3>(1.0f,0.0f,1.0f),
		GLColor<GLfloat,3>(1.0f,1.0f,1.0f),GLColor<GLfloat,3>(1.0f,1.0f,1.0f)
		};
	glTexParameteri(GL_TEXTURE_1D,GL_TEXTURE_MIN_FILTER,GL_NEAREST);
	glTexParameteri(GL_TEXTURE_1D,GL_TEXTURE_MAG_FILTER,GL_NEAREST);
	glTexParameteri(GL_TEXTURE_1D,GL_TEXTURE_WRAP_S,GL_CLAMP);
	glTexImage1D(GL_TEXTURE_1D,0,GL_RGB,8,0,GL_RGB,GL_FLOAT,dilationColorMap);
	#endif
	
	/* Get the data item: */
	DataItem* dataItem=contextData.retrieveDataItem<DataItem>(this);
	if(dataItem->demBufferIds[0]!=0)
		{
		/* Bind the vertex buffer and index buffer: */
		glBindBufferARB(GL_ARRAY_BUFFER_ARB,dataItem->demBufferIds[0]);
		glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER_ARB,dataItem->demBufferIds[1]);
		
		/* Install the vertex array: */
		GLVertexArrayParts::enable(DemVertex::getPartsMask());
		glVertexPointer(static_cast<const DemVertex*>(0));
		
		/* Draw the DEM: */
		const GLuint* indexPtr=0;
		for(int y=1;y<dem.getSize(1);++y)
			{
			glDrawElements(GL_QUAD_STRIP,dem.getSize(0)*2,GL_UNSIGNED_INT,indexPtr);
			indexPtr+=dem.getSize(0)*2;
			}
		
		/* Disable the vertex array: */
		GLVertexArrayParts::disable(DemVertex::getPartsMask());
		
		/* Protect the vertex buffer and index buffer: */
		glBindBufferARB(GL_ARRAY_BUFFER_ARB,0);
		glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER_ARB,0);
		}
	else
		{
		/* Render the DEM: */
		for(int y=1;y<dem.getSize(1);++y)
			{
			glBegin(GL_QUAD_STRIP);
			for(int x=0;x<dem.getSize(0);++x)
				{
				glTexCoord1f(dem(x,y).texCoord[0]);
				glNormal(dem(x,y).normal);
				glVertex(dem(x,y).position);
				glTexCoord1f(dem(x,y-1).texCoord[0]);
				glNormal(dem(x,y-1).normal);
				glVertex(dem(x,y-1).position);
				}
			glEnd();
			}
		}
	
	/* Restore OpenGL state: */
	glPopAttrib();
	}

void LidarGridder::initContext(GLContextData& contextData) const
	{
	/* Create a context data item: */
	DataItem* dataItem=new DataItem;
	contextData.addDataItem(this,dataItem);
	
	if(dataItem->demBufferIds[0]!=0)
		{
		/* Upload the DEM to the vertex buffer object: */
		glBindBufferARB(GL_ARRAY_BUFFER_ARB,dataItem->demBufferIds[0]);
		glBufferDataARB(GL_ARRAY_BUFFER_ARB,size_t(dem.getSize(0))*size_t(dem.getSize(1))*sizeof(DemVertex),dem.getArray(),GL_STATIC_DRAW_ARB);
		glBindBufferARB(GL_ARRAY_BUFFER_ARB,0);
		
		/* Upload an index set into the index buffer object: */
		glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER_ARB,dataItem->demBufferIds[1]);
		glBufferDataARB(GL_ELEMENT_ARRAY_BUFFER_ARB,size_t(dem.getSize(1)-1)*size_t(dem.getSize(0)*2)*sizeof(GLuint),0,GL_STATIC_DRAW_ARB);
		GLuint* indexPtr=static_cast<GLuint*>(glMapBufferARB(GL_ELEMENT_ARRAY_BUFFER_ARB,GL_WRITE_ONLY_ARB));
		for(int y=1;y<dem.getSize(1);++y)
			for(int x=0;x<dem.getSize(0);++x,indexPtr+=2)
				{
				indexPtr[0]=dem.calcLinearIndex(x,y);
				indexPtr[1]=dem.calcLinearIndex(x,y-1);
				}
		glUnmapBufferARB(GL_ELEMENT_ARRAY_BUFFER_ARB);
		glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER_ARB,0);
		}
	}

int main(int argc,char* argv[])
	{
	try
		{
		char** appDefaults=0;
		LidarGridder app(argc,argv,appDefaults);
		app.run();
		}
	catch(std::runtime_error err)
		{
		std::cerr<<"Caught exception "<<err.what()<<std::endl;
		return 1;
		}
	
	return 0;
	}

#if 0

int main(int argc,char* argv[])
	{
	const char* lidarFileName=0;
	const char* gridFileName=0;
	double gridOrigin[2]={0.0,0.0};
	int gridSize[2]={128,128};
	double gridCellSize[2]={1.0,1.0};
	int numLobes=3;
	int cacheSize=512;
	for(int i=1;i<argc;++i)
		{
		if(argv[i][0]=='-')
			{
			if(strcasecmp(argv[i]+1,"gridOrigin")==0)
				{
				for(int j=0;j<2;++j)
					{
					++i;
					gridOrigin[j]=atof(argv[i]);
					}
				}
			else if(strcasecmp(argv[i]+1,"gridSize")==0)
				{
				for(int j=0;j<2;++j)
					{
					++i;
					gridSize[j]=atoi(argv[i]);
					}
				}
			else if(strcasecmp(argv[i]+1,"gridCellSize")==0)
				{
				for(int j=0;j<2;++j)
					{
					++i;
					gridCellSize[j]=atof(argv[i]);
					}
				}
			else if(strcasecmp(argv[i]+1,"numLobes")==0)
				{
				++i;
				numLobes=atoi(argv[i]);
				}
			else if(strcasecmp(argv[i]+1,"cache")==0)
				{
				++i;
				cacheSize=atoi(argv[i]);
				}
			}
		else if(lidarFileName==0)
			lidarFileName=argv[i];
		else if(gridFileName==0)
			gridFileName=argv[i];
		}
	if(lidarFileName==0)
		{
		std::cerr<<"No LiDAR file name provided"<<std::endl;
		return 1;
		}
	if(gridFileName==0)
		{
		std::cerr<<"No grid file name provided"<<std::endl;
		return 1;
		}
	
	/* Open the LiDAR data set: */
	LidarProcessOctree lpo(lidarFileName,cacheSize*1024*1024);
	const Cube& domain=lpo.getDomain();
	
	/* Sample the grid: */
	std::cout<<"Sampling grid...   0%"<<std::flush;
	float* grid=new float[gridSize[1]*gridSize[0]];
	float* weights=new float[gridSize[1]*gridSize[0]];
	for(int y=0;y<gridSize[1];++y)
		{
		for(int x=0;x<gridSize[0];++x)
			{
			/* Calculate the grid vertex position: */
			double pos[2];
			pos[0]=gridOrigin[0]+double(x)*gridCellSize[0];
			pos[1]=gridOrigin[1]+double(y)*gridCellSize[1];
			
			/* Create a sampler: */
			Sampler sampler(pos,gridCellSize,numLobes);
			
			/* Calculate the sampler's bounding box: */
			Box sampleBox;
			for(int i=0;i<2;++i)
				{
				sampleBox.min[i]=Scalar(pos[i]-gridCellSize[i]*double(numLobes));
				sampleBox.max[i]=Scalar(pos[i]+gridCellSize[i]*double(numLobes));
				}
			sampleBox.min[2]=domain.getMin()[2];
			sampleBox.max[2]=domain.getMax()[2];
			
			/* Sample the LiDAR data set: */
			lpo.processPointsInBox(sampleBox,sampler);
			grid[y*gridSize[0]+x]=float(sampler.getValue());
			weights[y*gridSize[0]+x]=float(Math::abs(sampler.getWeightSum()));
			}
		std::cout<<"\b\b\b\b"<<std::setw(3)<<(100*(y+1)+gridSize[1]/2)/gridSize[1]<<"%"<<std::flush;
		}
	std::cout<<std::endl;
	
	/* Save the grid: */
	{
	Misc::File headerFile((std::string(gridFileName)+".hdr").c_str(),"wt");
	fprintf(headerFile.getFilePtr(),"BANDGAPBYTES 0\n");
	fprintf(headerFile.getFilePtr(),"BANDROWBYTES %u\n",size_t(gridSize[0])*sizeof(float));
	fprintf(headerFile.getFilePtr(),"BYTEORDER I\n");
	fprintf(headerFile.getFilePtr(),"LAYOUT BIL\n");
	fprintf(headerFile.getFilePtr(),"NBANDS 1\n");
	fprintf(headerFile.getFilePtr(),"NBITS %u\n",sizeof(float)*8);
	fprintf(headerFile.getFilePtr(),"NCOLS %u\n",size_t(gridSize[0]));
	fprintf(headerFile.getFilePtr(),"NROWS %u\n",size_t(gridSize[1]));
	fprintf(headerFile.getFilePtr(),"TOTALROWBYTES %u\n",size_t(gridSize[0])*sizeof(float));
	fprintf(headerFile.getFilePtr(),"XDIM %lf\n",gridCellSize[0]);
	fprintf(headerFile.getFilePtr(),"YDIM %lf\n",gridCellSize[1]);
	}
	
	{
	Misc::File gridFile((std::string(gridFileName)+".bil").c_str(),"wb",Misc::File::LittleEndian);
	for(int y=gridSize[1]-1;y>=0;--y)
		gridFile.write(&grid[y*gridSize[0]],gridSize[0]);
	}
	
	/* Save the weight sum image: */
	float weightMin,weightMax;
	weightMin=weightMax=weights[0];
	for(int i=1;i<gridSize[1]*gridSize[0];++i)
		{
		if(weightMin>weights[i])
			weightMin=weights[i];
		if(weightMax<weights[i])
			weightMax=weights[i];
		}
	std::cout<<"Weight sum range: ["<<weightMin<<", "<<weightMax<<"]"<<std::endl;
	
	{
	Misc::File weightFile("Weights.ppm","wb",Misc::File::DontCare);
	fprintf(weightFile.getFilePtr(),"P5\n%d %d\n255\n",gridSize[0],gridSize[1]);
	unsigned char* row=new unsigned char[gridSize[0]];
	for(int y=gridSize[1]-1;y>=0;--y)
		{
		for(int x=0;x<gridSize[0];++x)
			row[x]=(unsigned char)((weights[y*gridSize[0]+x]-weightMin)/(weightMax-weightMin)*255.0f+0.5f);
		weightFile.write(row,gridSize[0]);
		}
	delete[] row;
	}
	
	/* Calculate the difference image between the reconstructed grid and the LiDAR point cloud: */
	DifferenceCalculator dc(gridOrigin,gridSize,gridCellSize,grid);
	Box gridBox;
	for(int i=0;i<2;++i)
		{
		gridBox.min[i]=Scalar(gridOrigin[i]);
		gridBox.max[i]=Scalar(gridOrigin[i]+double(gridSize[i]-1)*gridCellSize[i]);
		}
	gridBox.min[2]=domain.getMin()[2];
	gridBox.max[2]=domain.getMax()[2];
	lpo.processPointsInBox(gridBox,dc);
	float* diffs=new float[gridSize[1]*gridSize[0]];
	float diffMax=0.0f;
	for(int y=0;y<gridSize[1];++y)
		for(int x=0;x<gridSize[0];++x)
			{
			if(dc.getDiffWeights()[y*gridSize[0]+x]>0)
				{
				float d=dc.getDiffRMS(x,y);
				diffs[y*gridSize[0]+x]=d;
				if(diffMax<d)
					diffMax=d;
				}
			else
				diffs[y*gridSize[0]+x]=0.0f;
			}
	std::cout<<"Max difference RMS per grid cell: "<<diffMax<<std::endl;
	
	{
	Misc::File diffFile("Diffs.ppm","wb",Misc::File::DontCare);
	fprintf(diffFile.getFilePtr(),"P5\n%d %d\n255\n",gridSize[0],gridSize[1]);
	unsigned char* row=new unsigned char[gridSize[0]];
	for(int y=gridSize[1]-1;y>=0;--y)
		{
		for(int x=0;x<gridSize[0];++x)
			row[x]=(unsigned char)(diffs[y*gridSize[0]+x]/diffMax*255.0f+0.5f);
		diffFile.write(row,gridSize[0]);
		}
	delete[] row;
	}
	
	/* Clean up and exit: */
	delete[] grid;
	delete[] weights;
	delete[] diffs;
	return 0;
	}

#endif
