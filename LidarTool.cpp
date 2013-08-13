/***********************************************************************
LidarTool - Vrui tool class to position a virtual input device at the
intersection of a ray and a LiDAR octree.
Copyright (c) 2008-2010 Oliver Kreylos

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

#include "LidarTool.h"

#include <Cluster/MulticastPipe.h>
#include <Geometry/Ray.h>
#include <Geometry/OrthogonalTransformation.h>
#include <GL/gl.h>
#include <GL/GLGeometryWrappers.h>
#include <Vrui/InputDevice.h>
#include <Vrui/ToolManager.h>
#include <Vrui/InputGraphManager.h>
#include <Vrui/Vrui.h>

#include "LidarTypes.h"
#include "LidarOctree.h"
#include "LidarViewer.h"

/*********************************
Methods of class LidarToolFactory:
*********************************/

LidarToolFactory::LidarToolFactory(Vrui::ToolManager& toolManager)
	:Vrui::ToolFactory("LidarTool",toolManager)
	{
	/* Insert class into class hierarchy: */
	Vrui::TransformToolFactory* transformToolFactory=dynamic_cast<Vrui::TransformToolFactory*>(toolManager.loadClass("TransformTool"));
	transformToolFactory->addChildClass(this);
	addParentClass(transformToolFactory);
	
	/* Initialize tool layout: */
	layout.setNumButtons(0,true);
	layout.setNumValuators(0,true);
	
	/* Set the custom tool class' factory pointer: */
	LidarTool::factory=this;
	}

LidarToolFactory::~LidarToolFactory(void)
	{
	/* Reset the custom tool class' factory pointer: */
	LidarTool::factory=0;
	}

const char* LidarToolFactory::getName(void) const
	{
	return "Point Cloud Projector";
	}

Vrui::Tool* LidarToolFactory::createTool(const Vrui::ToolInputAssignment& inputAssignment) const
	{
	/* Create a new object of the custom tool class: */
	LidarTool* newTool=new LidarTool(this,inputAssignment);
	
	return newTool;
	}

void LidarToolFactory::destroyTool(Vrui::Tool* tool) const
	{
	/* Cast the tool pointer to the Lidar tool class (not really necessary): */
	LidarTool* lidarTool=dynamic_cast<LidarTool*>(tool);
	
	/* Destroy the tool: */
	delete lidarTool;
	}

/**********************************
Static elements of class LidarTool:
**********************************/

LidarToolFactory* LidarTool::factory=0;

/**************************
Methods of class LidarTool:
**************************/

LidarTool::LidarTool(const Vrui::ToolFactory* factory,const Vrui::ToolInputAssignment& inputAssignment)
	:Vrui::TransformTool(factory,inputAssignment)
	{
	/* Set the source device: */
	if(input.getNumButtonSlots()>0)
		sourceDevice=getButtonDevice(0);
	else
		sourceDevice=getValuatorDevice(0);
	}

void LidarTool::initialize(void)
	{
	/* Initialize the base tool: */
	TransformTool::initialize();
	
	/* Disable the transformed device's glyph: */
	Vrui::getInputGraphManager()->getInputDeviceGlyph(transformedDevice).disable();
	}

const Vrui::ToolFactory* LidarTool::getFactory(void) const
	{
	return factory;
	}

void LidarTool::frame(void)
	{
	/* Calculate ray equation in navigation coordinates: */
	LidarOctree::Ray deviceRay=sourceDevice->getRay();
	LidarOctree::Ray modelRay=deviceRay;
	modelRay.transform(Vrui::getInverseNavigationTransformation());
	
	/* Intersect the ray with the LiDAR data set: */
	Scalar rayParameter;
	if(Vrui::isMaster())
		{
		/* Calculate the intersection: */
		rayParameter=Scalar(-1);
		LidarOctree::ConeIntersection cone(modelRay,Vrui::getRayPickCosine());
		for(int i=0;i<application->numOctrees;++i)
			{
			application->octrees[i]->intersectCone(cone);
			if(cone.isValid())
				{
				rayParameter=cone.getParameter();
				cone.testLambda2=cone.testLambdaMin;
				}
			}
		
		if(Vrui::getMainPipe()!=0)
			{
			/* Send the intersection to the slaves: */
			Vrui::getMainPipe()->write<Scalar>(rayParameter);
			// Vrui::getMainPipe()->flush();
			}
		}
	else
		{
		/* Receive the intersection from the master: */
		Vrui::getMainPipe()->read<Scalar>(rayParameter);
		}
	
	transformedDevice->setDeviceRay(sourceDevice->getDeviceRayDirection(),sourceDevice->getDeviceRayStart());
	if(rayParameter>=Scalar(0))
		{
		/* Set the device position to the intersection point: */
		Vrui::TrackerState ts(Vrui::getNavigationTransformation().transform(modelRay(rayParameter))-Vrui::Point::origin,sourceDevice->getOrientation());
		transformedDevice->setTransformation(ts);
		}
	else
		{
		/* Move the device in the plane it currently inhabits: */
		rayParameter=(deviceRay.getDirection()*Vector(transformedDevice->getPosition()-sourceDevice->getPosition()))/Geometry::sqr(deviceRay.getDirection());
		Vrui::TrackerState ts(deviceRay(rayParameter)-Point::origin,sourceDevice->getOrientation());
		transformedDevice->setTransformation(ts);
		}
	}

void LidarTool::display(GLContextData& contextData) const
	{
	/* Draw a line from the source device's position to the transformed device's position: */
	glPushAttrib(GL_ENABLE_BIT|GL_LINE_BIT);
	glDisable(GL_LIGHTING);
	glLineWidth(1.0f);
	
	glColor3f(0.0f,1.0f,0.0f);
	glBegin(GL_LINES);
	glVertex(sourceDevice->getPosition());
	glVertex(transformedDevice->getPosition());
	glEnd();
	
	glPopAttrib();
	}
