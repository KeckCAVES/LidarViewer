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

#ifndef PROFILETOOL_INCLUDED
#define PROFILETOOL_INCLUDED

#include <GLMotif/TextFieldSlider.h>
#include <Vrui/Tool.h>
#include <Vrui/Application.h>

#include "LidarTypes.h"

/* Forward declarations: */
namespace GLMotif {
class PopupWindow;
}
class LidarViewer;

class ProfileTool; // Forward declaration of the profile tool class

class ProfileToolFactory:public Vrui::ToolFactory // Class for factories that create/destroy LiDAR tool objects
	{
	friend class ProfileTool;
	
	/* Constructors and destructors: */
	public:
	ProfileToolFactory(Vrui::ToolManager& toolManager);
	static void factoryDestructor(Vrui::ToolFactory* factory)
		{
		delete factory;
		}
	virtual ~ProfileToolFactory(void);
	
	/* Methods from ToolFactory: */
	virtual const char* getName(void) const;
	virtual const char* getButtonFunction(int buttonSlotIndex) const;
	virtual Vrui::Tool* createTool(const Vrui::ToolInputAssignment& inputAssignment) const;
	virtual void destroyTool(Vrui::Tool* tool) const;
	};

class ProfileTool:public Vrui::Tool,public Vrui::Application::Tool<LidarViewer> // The profile tool class
	{
	friend class ProfileToolFactory;
	
	/* Elements: */
	private:
	static ProfileToolFactory* factory; // Pointer to the factory object for this class
	double segmentLength;
	int oversampling;
	double filterWidth;
	int numLobes;
	GLMotif::PopupWindow* settingsDialogPopup; // Dialog to adjust the profile extraction settings
	Point p0,p1; // First and second profile points
	int state; // Tool state (0: idle, 1: dragging first point, 2: waiting, 3: dragging second point, 4: extracting profile
	
	/* Private methods: */
	void segmentLengthCallback(GLMotif::TextFieldSlider::ValueChangedCallbackData* cbData);
	void oversamplingCallback(GLMotif::TextFieldSlider::ValueChangedCallbackData* cbData);
	void filterWidthCallback(GLMotif::TextFieldSlider::ValueChangedCallbackData* cbData);
	void numLobesCallback(GLMotif::TextFieldSlider::ValueChangedCallbackData* cbData);
	
	/* Constructors and destructors: */
	public:
	ProfileTool(const Vrui::ToolFactory* factory,const Vrui::ToolInputAssignment& inputAssignment);
	virtual ~ProfileTool(void);
	
	/* Methods from Tool: */
	virtual const Vrui::ToolFactory* getFactory(void) const;
	virtual void buttonCallback(int buttonSlotIndex,Vrui::InputDevice::ButtonCallbackData* cbData);
	virtual void frame(void);
	virtual void display(GLContextData& contextData) const;
	};

#endif
