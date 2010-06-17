/***********************************************************************
LoadPointSet - Helper function to load a 2D point set into LiDAR Viewer
by elevating the points to the implicit LiDAR point cloud surface.
Copyright (c) 2010 Oliver Kreylos

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

#include "LoadPointSet.h"

#include <string>
#include <Misc/FileCharacterSource.h>
#include <Misc/ValueSource.h>
#include <SceneGraph/ColorNode.h>
#include <SceneGraph/CoordinateNode.h>
#include <SceneGraph/PointSetNode.h>
#include <SceneGraph/ShapeNode.h>

#include "LidarProcessOctree.h"
#include "LidarElevationSampler.h"
#include "SceneGraph.h"

void loadPointSet(const char* lidarFileName,unsigned int memCacheSize,const char* pointFileName,Scalar filterRadius,int numLobes,const Vector& pointOffset)
	{
	/* Create a color and coordinate node to receive points: */
	SceneGraph::ColorNode* color=new SceneGraph::ColorNode;
	SceneGraph::CoordinateNode* coord=new SceneGraph::CoordinateNode;
	
	/* Create the point set node: */
	SceneGraph::PointSetNode* pointSet=new SceneGraph::PointSetNode;
	pointSet->color.setValue(color);
	pointSet->coord.setValue(coord);
	pointSet->pointSize.setValue(5);
	pointSet->update();
	SceneGraph::ShapeNode* shape=new SceneGraph::ShapeNode;
	shape->geometry.setValue(pointSet);
	shape->update();
	getSceneGraphRoot().children.appendValue(shape);
	getSceneGraphRoot().update();
	
	/* Create a temporary LiDAR processing octree to calculate point elevations: */
	LidarProcessOctree lpo(lidarFileName,size_t(memCacheSize)*size_t(1024*1024));
	
	/* Open the point file: */
	Misc::FileCharacterSource pointFile(pointFileName);
	Misc::ValueSource pointSource(pointFile);
	pointSource.setWhitespace("");
	pointSource.setPunctuation(",\n");
	pointSource.setQuotes("\"");
	
	/* Read the header line: */
	int valueCol=-1;
	int posCol[2]={-1,-1};
	for(int col=0;!pointSource.eof()&&pointSource.peekc()!='\n';++col)
		{
		std::string columnHeader=pointSource.readString();
		if(columnHeader=="DamageID1")
			valueCol=col;
		else if(columnHeader=="POINT_X")
			posCol[0]=col;
		else if(columnHeader=="POINT_Y")
			posCol[1]=col;
		if(pointSource.peekc()==',')
			pointSource.skipString();
		}
	
	/* Go to the next line: */
	pointSource.skipString();
	
	/* Read all points: */
	while(!pointSource.eof())
		{
		/* Read the current line: */
		double value=0.0;
		double pos[2]={0.0,0.0};
		for(int col=0;!pointSource.eof()&&pointSource.peekc()!='\n';++col)
			{
			if(col==valueCol)
				value=pointSource.readNumber();
			else if(col==posCol[0])
				pos[0]=pointSource.readNumber();
			else if(col==posCol[1])
				pos[1]=pointSource.readNumber();
			else
				pointSource.skipString();
			if(pointSource.peekc()==',')
				pointSource.skipString();
			}
		
		/* Calculate the current point's elevation: */
		LidarRadialElevationSampler les(pos,filterRadius,numLobes);
		lpo.processPointsInBox(les.getBox(),les);
		if(les.getAbsWeightSum()>0.0)
			{
			/* Add the point to the coordinate node: */
			SceneGraph::Point p;
			for(int i=0;i<2;++i)
				p[i]=SceneGraph::Scalar(pos[i]-double(pointOffset[i]));
			p[2]=SceneGraph::Scalar(les.getValue()-double(pointOffset[2]));
			coord->point.appendValue(p);
			
			if(value<2.5)
				color->color.appendValue(SceneGraph::Color(1.0f,0.0f,0.0f));
			else if(value<3.5)
				color->color.appendValue(SceneGraph::Color(1.0f,0.5f,0.0f));
			else if(value<4.5)
				color->color.appendValue(SceneGraph::Color(0.5f,1.0f,0.0f));
			else if(value<5.5)
				color->color.appendValue(SceneGraph::Color(0.0f,1.0f,0.0f));
			else
				color->color.appendValue(SceneGraph::Color(0.0f,1.0f,0.0f));
			}
		
		/* Go to the next line: */
		pointSource.skipString();
		}
	
	color->update();
	coord->update();
	}
