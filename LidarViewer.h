/***********************************************************************
LidarViewer - Viewer program for multiresolution LiDAR data.
Copyright (c) 2005-2008 Oliver Kreylos

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

#ifndef LIDARVIEWER_INCLUDED
#define LIDARVIEWER_INCLUDED

#include <vector>
#include <Misc/CallbackData.h>
#include <Geometry/Plane.h>
#include <GL/gl.h>
#include <GL/GLColor.h>
#include <GL/GLObject.h>
#include <GLMotif/RadioBox.h>
#include <GLMotif/ToggleButton.h>
#include <GLMotif/Slider.h>
#include <Vrui/Tools/LocatorTool.h>
#include <Vrui/LocatorToolAdapter.h>
#include <Vrui/Tools/DraggingTool.h>
#include <Vrui/Vrui.h>
#include <Vrui/Application.h>

#include "LidarTypes.h"
#include "PointBasedLightingShader.h"

/* Forward declarations: */
namespace Comm {
class MulticastPipe;
}
namespace GLMotif {
class Popup;
class PopupMenu;
class PopupWindow;
class TextField;
}
namespace Vrui {
class InputDevice;
class Lightsource;
}
class LidarOctree;
class Primitive;

class LidarViewer:public Vrui::Application,public GLObject
	{
	/* Embedded classes: */
	private:
	typedef Geometry::Plane<double,3> GPlane;
	
	class Locator:public Vrui::LocatorToolAdapter // Base class for application locators
		{
		/* Elements: */
		protected:
		LidarViewer* application; // Pointer to the application object
		
		/* Constructors and destructors: */
		Locator(Vrui::LocatorTool* sTool,LidarViewer* sApplication)
			:LocatorToolAdapter(sTool),
			 application(sApplication)
			{
			};
		public:
		virtual ~Locator(void)
			{
			};
		
		/* Methods: */
		virtual void updateSettings(void)
			{
			};
		virtual void glRenderAction(GLContextData& contextData) const
			{
			};
		};
	
	class SelectorLocator:public Locator // Class for point selection locators
		{
		friend class LidarViewer;
		
		/* Embedded classes: */
		public:
		enum SelectorMode // Enumerated type for selection modes
			{
			ADD,SUBTRACT
			};
		
		/* Elements: */
		private:
		Vrui::Scalar influenceRadius; // Radius of selector's influence in physical coordinates
		SelectorMode selectorMode; // Selection mode of this locator
		Vrui::NavTrackerState locatorTransform; // Current transformation of the locator
		Point modelCenter; // Locator's position in model coordinates
		Scalar modelRadius; // Locator's influence radius in model coordinates
		bool ready; // Glag if the selector has been properly initialized (modelCenter and modelRadius are valid)
		bool active; // Flag whether the selector is active (selecting points)
		
		/* Constructors and destructors: */
		public:
		SelectorLocator(Vrui::LocatorTool* sTool,LidarViewer* sApplication);
		virtual ~SelectorLocator(void);
		
		/* Methods: */
		virtual void motionCallback(Vrui::LocatorTool::MotionCallbackData* cbData);
		virtual void buttonPressCallback(Vrui::LocatorTool::ButtonPressCallbackData* cbData);
		virtual void buttonReleaseCallback(Vrui::LocatorTool::ButtonReleaseCallbackData* cbData);
		virtual void updateSettings();
		virtual void glRenderAction(GLContextData& contextData) const;
		};
	
	typedef std::vector<Locator*> LocatorList;
	typedef std::vector<Primitive*> PrimitiveList;
	
	struct DataItem:public GLObject::DataItem
		{
		/* Elements: */
		public:
		GLuint influenceSphereDisplayListId; // ID of display list to render transparent spheres
		GLuint planeColorMapTextureId; // Texture object ID of texture plane color map
		PointBasedLightingShader pbls; // Shader for point-based lighting
		
		/* Constructors and destructors: */
		DataItem(void);
		virtual ~DataItem(void);
		};
	
	friend class SelectorLocator;
	
	/* Elements: */
	LidarOctree* octree; // The LiDAR data representation
	Scalar renderQuality; // The current rendering quality (adapted to achieve optimal frame rate)
	Scalar fncWeight; // Weight factor for focus+context LOD adjustment
	float pointSize; // The pixel size used to render LiDAR points
	bool pointBasedLighting; // Flag whether points are rendered with illumination
	bool usePointColors; // Flag whether to use points' colors during illuminated rendering
	bool enableSun; // Flag whether to use a sun light source instead of all viewer's headlights
	bool* viewerHeadlightStates; // Enable states of all viewers' headlights at the last time the sun light source was turned on
	Vrui::Scalar sunAzimuth,sunElevation; // Azimuth and elevation angles of sun light source in degrees
	Vrui::Lightsource* sun; // Light source representing the sun
	bool useTexturePlane; // Flag whether to use automatically generated texture coordinates to visualize point distance from a plane
	GPlane texturePlane; // Plane equation of the texture-generating plane
	double texturePlaneScale; // Scale factor for texture plane distances
	bool updateTree; // Flag if the tree is continuously updated
	double lastFrameTime; // Application time of last frame; used to calculate render performance
	
	/* Interaction state: */
	bool overrideTools; // Flag whether interaction settings changes influence existing tools
	Vrui::Scalar brushSize; // Default physical-coordinate size for new interaction brushes
	GLColor<GLfloat,4> brushColor; // Color to render selection brush, with transparency
	SelectorLocator::SelectorMode defaultSelectorMode; // Selection mode for new selector locators
	Scalar neighborhoodSize; // Size of neighborhood for point classification
	LocatorList locators; // List of currently existing locators
	Comm::MulticastPipe* extractorPipe; // Pipe to synchronize feature extraction on a distributed rendering cluster
	GLColor<GLfloat,4> primitiveColor; // Color to render primitives, with transparency
	GLColor<GLfloat,4> selectedPrimitiveColor; // Color to render selected primitives, with transparency
	PrimitiveList primitives; // List of extracted primitives
	int lastPickedPrimitive; // Index of the most recently picked primitive
	std::vector<bool> primitiveSelectedFlags; // List of selected flags for extracted primitives
	
	/* Vrui state: */
	GLMotif::PopupMenu* mainMenu; // The program's main menu
	GLMotif::RadioBox* mainMenuSelectorModes;
	GLMotif::PopupWindow* renderDialog; // The rendering settings dialog
	GLMotif::TextField* renderQualityValue;
	GLMotif::Slider* renderQualitySlider;
	GLMotif::TextField* fncWeightValue;
	GLMotif::Slider* fncWeightSlider;
	GLMotif::TextField* pointSizeValue;
	GLMotif::Slider* pointSizeSlider;
	GLMotif::TextField* sunAzimuthValue;
	GLMotif::Slider* sunAzimuthSlider;
	GLMotif::TextField* sunElevationValue;
	GLMotif::Slider* sunElevationSlider;
	GLMotif::TextField* texturePlaneScaleValue;
	GLMotif::Slider* texturePlaneScaleSlider;
	GLMotif::PopupWindow* interactionDialog; // The interaction settings dialog
	GLMotif::RadioBox* interactionDialogSelectorModes;
	GLMotif::TextField* brushSizeValue;
	GLMotif::Slider* brushSizeSlider;
	GLMotif::TextField* neighborhoodSizeValue;
	GLMotif::Slider* neighborhoodSizeSlider;
	
	/* Private methods: */
	GLMotif::Popup* createSelectorModesMenu(void);
	GLMotif::Popup* createSelectionMenu(void);
	GLMotif::Popup* createExtractionMenu(void);
	GLMotif::Popup* createDialogMenu(void);
	GLMotif::PopupMenu* createMainMenu(void);
	GLMotif::PopupWindow* createRenderDialog(void);
	GLMotif::PopupWindow* createInteractionDialog(void);
	static void treeUpdateNotificationCB(void* userData)
		{
		Vrui::requestUpdate();
		};
	void draggingToolCallback(Vrui::DraggingTool::DragStartCallbackData* cbData); // Callback for selection events from dragging tools
	int addPrimitive(Primitive* newPrimitive); // Adds a primitive to the list; returns the index of the newly added primitive
	void selectPrimitive(int primitiveIndex); // Selects the given primitive
	void deselectPrimitive(int primitiveIndex); // Deselects the given primitive
	void deletePrimitive(int primitiveIndex); // Deletes the given primitive from the list
	void updateSun(void); // Updates the state of the sun light source
	void setEnableSun(bool newEnableSun); // Enables or disables the sun light source
	
	/* Constructors and destructors: */
	public:
	LidarViewer(int& argc,char**& argv,char**& appDefaults);
	virtual ~LidarViewer(void);
	
	/* Methods: */
	virtual void initContext(GLContextData& contextData) const;
	virtual void toolCreationCallback(Vrui::ToolManager::ToolCreationCallbackData* cbData);
	virtual void toolDestructionCallback(Vrui::ToolManager::ToolDestructionCallbackData* cbData);
	virtual void frame(void);
	virtual void display(GLContextData& contextData) const;
	void centerDisplayCallback(Misc::CallbackData* cbData);
	void changeSelectorModeCallback(GLMotif::RadioBox::ValueChangedCallbackData* cbData);
	void classifySelectionCallback(Misc::CallbackData* cbData);
	void saveSelectionCallback(Misc::CallbackData* cbData);
	void clearSelectionCallback(Misc::CallbackData* cbData);
	void extractPlaneCallback(Misc::CallbackData* cbData);
	void extractBruntonCallback(Misc::CallbackData* cbData);
	void extractSphereCallback(Misc::CallbackData* cbData);
	void extractCylinderCallback(Misc::CallbackData* cbData);
	void intersectPrimitivesCallback(Misc::CallbackData* cbData);
	void loadPrimitivesCallback(Misc::CallbackData* cbData);
	void savePrimitivesCallback(Misc::CallbackData* cbData);
	void deleteSelectedPrimitivesCallback(Misc::CallbackData* cbData);
	void clearPrimitivesCallback(Misc::CallbackData* cbData);
	void showRenderDialogCallback(GLMotif::ToggleButton::ValueChangedCallbackData* cbData);
	void renderQualitySliderCallback(GLMotif::Slider::ValueChangedCallbackData* cbData);
	void fncWeightSliderCallback(GLMotif::Slider::ValueChangedCallbackData* cbData);
	void pointSizeSliderCallback(GLMotif::Slider::ValueChangedCallbackData* cbData);
	void enableLightingCallback(GLMotif::ToggleButton::ValueChangedCallbackData* cbData);
	void usePointColorsCallback(GLMotif::ToggleButton::ValueChangedCallbackData* cbData);
	void enableSunCallback(GLMotif::ToggleButton::ValueChangedCallbackData* cbData);
	void sunAzimuthSliderCallback(GLMotif::Slider::ValueChangedCallbackData* cbData);
	void sunElevationSliderCallback(GLMotif::Slider::ValueChangedCallbackData* cbData);
	void enableTexturePlaneCallback(GLMotif::ToggleButton::ValueChangedCallbackData* cbData);
	void texturePlaneScaleSliderCallback(GLMotif::Slider::ValueChangedCallbackData* cbData);
	void showInteractionDialogCallback(GLMotif::ToggleButton::ValueChangedCallbackData* cbData);
	void overrideToolsCallback(GLMotif::ToggleButton::ValueChangedCallbackData* cbData);
	void brushSizeSliderCallback(GLMotif::Slider::ValueChangedCallbackData* cbData);
	void neighborhoodSizeSliderCallback(GLMotif::Slider::ValueChangedCallbackData* cbData);
	void updateTreeCallback(GLMotif::ToggleButton::ValueChangedCallbackData* cbData);
	};

#endif
