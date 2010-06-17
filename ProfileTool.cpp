/***********************************************************************
ProfileTool - Vrui tool class to extract profile curves from 2.5D LiDAR
data sets.
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

#include "ProfileTool.h"

#include <GL/gl.h>
#include <GL/GLGeometryWrappers.h>
#include <GL/GLTransformationWrappers.h>
#include <GLMotif/StyleSheet.h>
#include <GLMotif/WidgetManager.h>
#include <GLMotif/PopupWindow.h>
#include <GLMotif/RowColumn.h>
#include <GLMotif/Label.h>
#include <Vrui/DisplayState.h>
#include <Vrui/Vrui.h>

#include "LidarOctree.h"
#include "ProfileExtractor.h"

/***********************************
Methods of class ProfileToolFactory:
***********************************/

ProfileToolFactory::ProfileToolFactory(Vrui::ToolManager& toolManager,const LidarOctree* sOctree)
	:Vrui::ToolFactory("ProfileTool",toolManager),
	 octree(sOctree)
	{
	#if 0
	/* Insert class into class hierarchy: */
	Vrui::TransformToolFactory* transformToolFactory=dynamic_cast<Vrui::TransformToolFactory*>(toolManager.loadClass("TransformTool"));
	transformToolFactory->addChildClass(this);
	addParentClass(transformToolFactory);
	#endif
	
	/* Initialize tool layout: */
	layout.setNumDevices(1);
	layout.setNumButtons(0,1);
	
	/* Set the custom tool class' factory pointer: */
	ProfileTool::factory=this;
	}

ProfileToolFactory::~ProfileToolFactory(void)
	{
	/* Reset the custom tool class' factory pointer: */
	ProfileTool::factory=0;
	}

const char* ProfileToolFactory::getName(void) const
	{
	return "Profile Extractor";
	}

Vrui::Tool* ProfileToolFactory::createTool(const Vrui::ToolInputAssignment& inputAssignment) const
	{
	/* Create a new object of the custom tool class: */
	ProfileTool* newTool=new ProfileTool(this,inputAssignment);
	
	return newTool;
	}

void ProfileToolFactory::destroyTool(Vrui::Tool* tool) const
	{
	/* Cast the tool pointer to the Lidar tool class (not really necessary): */
	ProfileTool* lidarTool=dynamic_cast<ProfileTool*>(tool);
	
	/* Destroy the tool: */
	delete lidarTool;
	}

/************************************
Static elements of class ProfileTool:
************************************/

ProfileToolFactory* ProfileTool::factory=0;

/****************************
Methods of class ProfileTool:
****************************/

void ProfileTool::segmentLengthCallback(GLMotif::TextFieldSlider::ValueChangedCallbackData* cbData)
	{
	segmentLength=cbData->value;
	}

void ProfileTool::oversamplingCallback(GLMotif::TextFieldSlider::ValueChangedCallbackData* cbData)
	{
	oversampling=int(cbData->value+0.5);
	}

void ProfileTool::filterWidthCallback(GLMotif::TextFieldSlider::ValueChangedCallbackData* cbData)
	{
	filterWidth=cbData->value;
	}

void ProfileTool::numLobesCallback(GLMotif::TextFieldSlider::ValueChangedCallbackData* cbData)
	{
	numLobes=int(cbData->value+0.5);
	}

ProfileTool::ProfileTool(const Vrui::ToolFactory* factory,const Vrui::ToolInputAssignment& inputAssignment)
	:Vrui::Tool(factory,inputAssignment),
	 segmentLength(1.0),oversampling(2),filterWidth(1.0),numLobes(3),
	 settingsDialogPopup(0),
	 state(0)
	{
	/* Get the style sheet: */
	const GLMotif::StyleSheet* ss=Vrui::getWidgetManager()->getStyleSheet();
	
	/* Create the settings dialog window: */
	settingsDialogPopup=new GLMotif::PopupWindow("SettingsDialogPopup",Vrui::getWidgetManager(),"Profile Extraction Settings");
	settingsDialogPopup->setResizableFlags(true,false);
	
	GLMotif::RowColumn* settingsDialog=new GLMotif::RowColumn("Settings",settingsDialogPopup,false);
	settingsDialog->setNumMinorWidgets(2);
	
	/* Create a slider to adjust the segment length: */
	new GLMotif::Label("SegmentLengthLabel",settingsDialog,"Segment Length");
	
	GLMotif::TextFieldSlider* segmentLengthSlider=new GLMotif::TextFieldSlider("SegmentLengthSlider",settingsDialog,8,ss->fontHeight*10.0f);
	segmentLengthSlider->setSliderMapping(GLMotif::TextFieldSlider::EXP10);
	segmentLengthSlider->setValueRange(1.0e-3,1.0e3,0.1);
	segmentLengthSlider->setValue(segmentLength);
	segmentLengthSlider->getValueChangedCallbacks().add(this,&ProfileTool::segmentLengthCallback);
	
	/* Create a slider to adjust the segment length: */
	new GLMotif::Label("OversamplingLabel",settingsDialog,"Oversampling");
	
	GLMotif::TextFieldSlider* oversamplingSlider=new GLMotif::TextFieldSlider("OversamplingSlider",settingsDialog,8,ss->fontHeight*10.0f);
	oversamplingSlider->setSliderMapping(GLMotif::TextFieldSlider::LINEAR);
	oversamplingSlider->setValueRange(1.0,10.0,1.0);
	oversamplingSlider->setValue(oversampling);
	oversamplingSlider->getValueChangedCallbacks().add(this,&ProfileTool::oversamplingCallback);
	
	/* Create a slider to adjust the filter width: */
	new GLMotif::Label("FilterWidthLabel",settingsDialog,"Filter Width");
	
	GLMotif::TextFieldSlider* filterWidthSlider=new GLMotif::TextFieldSlider("FilterWidthSlider",settingsDialog,8,ss->fontHeight*10.0f);
	filterWidthSlider->setSliderMapping(GLMotif::TextFieldSlider::EXP10);
	filterWidthSlider->setValueRange(1.0e-3,1.0e3,0.1);
	filterWidthSlider->setValue(filterWidth);
	filterWidthSlider->getValueChangedCallbacks().add(this,&ProfileTool::filterWidthCallback);
	
	/* Create a slider to adjust the number of Lanczos filter lobes: */
	new GLMotif::Label("NumLobesLabel",settingsDialog,"Lanczos Filter Size");
	
	GLMotif::TextFieldSlider* numLobesSlider=new GLMotif::TextFieldSlider("NumLobesSlider",settingsDialog,8,ss->fontHeight*10.0f);
	numLobesSlider->setSliderMapping(GLMotif::TextFieldSlider::LINEAR);
	numLobesSlider->setValueRange(1.0,11.0,1.0);
	numLobesSlider->setValue(numLobes);
	numLobesSlider->getValueChangedCallbacks().add(this,&ProfileTool::numLobesCallback);
	
	settingsDialog->manageChild();
	
	Vrui::popupPrimaryWidget(settingsDialogPopup);
	}

ProfileTool::~ProfileTool(void)
	{
	delete settingsDialogPopup;
	}

const Vrui::ToolFactory* ProfileTool::getFactory(void) const
	{
	return factory;
	}

void ProfileTool::buttonCallback(int deviceIndex,int deviceButtonIndex,Vrui::InputDevice::ButtonCallbackData* cbData)
	{
	/* Change state depending on the button event: */
	switch(state)
		{
		case 0:
			if(cbData->newButtonState)
				state=1;
			break;
		
		case 1:
			if(!cbData->newButtonState)
				state=2;
			break;
		
		case 2:
			if(cbData->newButtonState)
				state=3;
			break;
		
		case 3:
			if(!cbData->newButtonState)
				state=4;
			break;
		}
	}

void ProfileTool::frame(void)
	{
	/* Get pointer to input device: */
	Vrui::InputDevice* iDevice=input.getDevice(0);
	
	/* Act depending on current state: */
	switch(state)
		{
		case 1:
			/* Update the first point: */
			p0=Point(Vrui::getInverseNavigationTransformation().transform(iDevice->getPosition()));
			break;
		
		case 3:
			/* Update the second point: */
			p1=Point(Vrui::getInverseNavigationTransformation().transform(iDevice->getPosition()));
			break;
		
		case 4:
			/* Extract the profile: */
			extractProfile(factory->octree,p0,p1,segmentLength,oversampling,filterWidth,numLobes,Vrui::getMainPipe());
			
			/* Go back to idle state: */
			state=0;
			break;
		}
	}

void ProfileTool::display(GLContextData& contextData) const
	{
	if(state>=1)
		{
		/* Save OpenGL state: */
		glPushAttrib(GL_ENABLE_BIT|GL_LINE_BIT|GL_POINT_BIT);
		glDisable(GL_LIGHTING);
		
		/* Go to navigation coordinates: */
		glPushMatrix();
		glLoadIdentity();
		glMultMatrix(Vrui::getDisplayState(contextData).modelviewNavigational);
		
		/* Draw stuff: */
		if(state>=1)
			{
			glPointSize(5.0f);
			glColor3f(0.5f,0.0f,0.0f);
			glBegin(GL_POINTS);
			glVertex(p0);
			if(state>=3)
				glVertex(p1);
			glEnd();
			
			if(state>=3)
				{
				glLineWidth(3.0f);
				glColor3f(0.0f,0.5f,0.0f);
				glBegin(GL_LINES);
				glVertex(p0);
				glVertex(p1);
				glEnd();
				}
			}
		
		/* Restore OpenGL state: */
		glPopMatrix();
		glPopAttrib();
		}
	}
