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

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <stdexcept>
#include <Misc/ThrowStdErr.h>
#include <Misc/File.h>
#include <Comm/MulticastPipe.h>
#include <Geometry/Vector.h>
#include <Geometry/TranslationTransformation.h>
#include <Geometry/OrthogonalTransformation.h>
#include <GL/gl.h>
#include <GL/GLContextData.h>
#include <GL/GLGeometryWrappers.h>
#include <GL/GLTransformationWrappers.h>
#include <GL/GLModels.h>
#include <GL/GLColorMap.h>
#include <GL/GLFrustum.h>
#include <GLMotif/StyleSheet.h>
#include <GLMotif/Button.h>
#include <GLMotif/CascadeButton.h>
#include <GLMotif/Label.h>
#include <GLMotif/TextField.h>
#include <GLMotif/RowColumn.h>
#include <GLMotif/Menu.h>
#include <GLMotif/SubMenu.h>
#include <GLMotif/Popup.h>
#include <GLMotif/PopupMenu.h>
#include <GLMotif/PopupWindow.h>
#include <GLMotif/WidgetManager.h>
#include <Vrui/GlyphRenderer.h>
#include <Vrui/OrthogonalCoordinateTransform.h>
#include <Vrui/CoordinateManager.h>
#include <Vrui/ToolManager.h>
#include <Vrui/ClusterSupport.h>

#include "LidarOctree.h"
#include "LidarTool.h"
#include "LidarSelectionSaver.h"
#include "Primitive.h"
#include "PlanePrimitive.h"
#include "LinePrimitive.h"
#include "PointPrimitive.h"
#include "SpherePrimitive.h"
#include "CylinderPrimitive.h"
#include "ArcInfoExportFile.h"

#include "LidarViewer.h"

/*********************************************
Methods of class LidarViewer::SelectorLocator:
*********************************************/

LidarViewer::SelectorLocator::SelectorLocator(Vrui::LocatorTool* sTool,LidarViewer* sApplication)
	:Locator(sTool,sApplication),
	 influenceRadius(application->brushSize),
	 selectorMode(application->defaultSelectorMode),
	 ready(false),
	 active(false)
	{
	}

LidarViewer::SelectorLocator::~SelectorLocator(void)
	{
	}

void LidarViewer::SelectorLocator::motionCallback(Vrui::LocatorTool::MotionCallbackData* cbData)
	{
	/* Update the locator's position and radius in model coordinates: */
	locatorTransform=cbData->currentTransformation;
	modelCenter=Point(locatorTransform.getOrigin());
	modelRadius=Scalar(influenceRadius*locatorTransform.getScaling());
	ready=true;
	
	if(active)
		{
		switch(selectorMode)
			{
			case ADD:
				application->octree->selectPoints(LidarOctree::Interactor(modelCenter,modelRadius));
				break;
			
			case SUBTRACT:
				application->octree->deselectPoints(LidarOctree::Interactor(modelCenter,modelRadius));
				break;
			}
		}
	}

void LidarViewer::SelectorLocator::buttonPressCallback(Vrui::LocatorTool::ButtonPressCallbackData* cbData)
	{
	active=true;
	}

void LidarViewer::SelectorLocator::buttonReleaseCallback(Vrui::LocatorTool::ButtonReleaseCallbackData* cbData)
	{
	active=false;
	}

void LidarViewer::SelectorLocator::updateSettings(void)
	{
	influenceRadius=application->brushSize;
	selectorMode=application->defaultSelectorMode;
	}

void LidarViewer::SelectorLocator::glRenderAction(GLContextData& contextData) const
	{
	glPushAttrib(GL_COLOR_BUFFER_BIT|GL_ENABLE_BIT|GL_LINE_BIT|GL_POLYGON_BIT);
	
	/* Retrieve context entry: */
	DataItem* dataItem=contextData.retrieveDataItem<DataItem>(application);
	
	/* Render the influence sphere: */
	glDisable(GL_LIGHTING);
	
	glPushMatrix();
	glMultMatrix(locatorTransform);
	glScale(influenceRadius);
	glCallList(dataItem->influenceSphereDisplayListId);
	glPopMatrix();
	
	glPopAttrib();
	}

/**************************************
Methods of class LidarViewer::DataItem:
**************************************/

LidarViewer::DataItem::DataItem(void)
	:influenceSphereDisplayListId(glGenLists(1))
	{
	glGenTextures(1,&planeColorMapTextureId);
	}

LidarViewer::DataItem::~DataItem(void)
	{
	glDeleteLists(influenceSphereDisplayListId,1);
	glDeleteTextures(1,&planeColorMapTextureId);
	}

/****************************
Methods of class LidarViewer:
****************************/

GLMotif::Popup* LidarViewer::createSelectorModesMenu(void)
	{
	GLMotif::Popup* selectorModesMenuPopup=new GLMotif::Popup("SelectorModesMenuPopup",Vrui::getWidgetManager());
	
	GLMotif::RadioBox* selectorModes=new GLMotif::RadioBox("SelectorModes",selectorModesMenuPopup,false);
	selectorModes->setSelectionMode(GLMotif::RadioBox::ALWAYS_ONE);
	
	selectorModes->addToggle("Add");
	selectorModes->addToggle("Subtract");
	
	selectorModes->manageChild();
	switch(defaultSelectorMode)
		{
		case SelectorLocator::ADD:
			selectorModes->setSelectedToggle(0);
			break;
		
		case SelectorLocator::SUBTRACT:
			selectorModes->setSelectedToggle(1);
			break;
		}
	selectorModes->getValueChangedCallbacks().add(this,&LidarViewer::changeSelectorModeCallback);
	
	return selectorModesMenuPopup;
	}

GLMotif::Popup* LidarViewer::createSelectionMenu(void)
	{
	GLMotif::Popup* selectionMenuPopup=new GLMotif::Popup("SelectionMenuPopup",Vrui::getWidgetManager());
	
	GLMotif::SubMenu* selectionMenu=new GLMotif::SubMenu("SelectionMenu",selectionMenuPopup,false);
	
	GLMotif::Button* clearSelectionButton=new GLMotif::Button("ClearSelectionButton",selectionMenu,"Clear Selection");
	clearSelectionButton->getSelectCallbacks().add(this,&LidarViewer::clearSelectionCallback);
	
	GLMotif::Button* saveSelectionButton=new GLMotif::Button("SaveSelectionButton",selectionMenu,"Save Selection");
	saveSelectionButton->getSelectCallbacks().add(this,&LidarViewer::saveSelectionCallback);
	
	selectionMenu->manageChild();
	
	return selectionMenuPopup;
	}

GLMotif::Popup* LidarViewer::createExtractionMenu(void)
	{
	GLMotif::Popup* extractionMenuPopup=new GLMotif::Popup("ExtractionMenuPopup",Vrui::getWidgetManager());
	
	GLMotif::SubMenu* extractionMenu=new GLMotif::SubMenu("ExtractionMenu",extractionMenuPopup,false);
	
	GLMotif::Button* extractPlaneButton=new GLMotif::Button("ExtractPlaneButton",extractionMenu,"Extract Plane");
	extractPlaneButton->getSelectCallbacks().add(this,&LidarViewer::extractPlaneCallback);
	
	GLMotif::Button* extractSphereButton=new GLMotif::Button("ExtractSphereButton",extractionMenu,"Extract Sphere");
	extractSphereButton->getSelectCallbacks().add(this,&LidarViewer::extractSphereCallback);
	
	GLMotif::Button* extractCylinderButton=new GLMotif::Button("ExtractCylinderButton",extractionMenu,"Extract Cylinder");
	extractCylinderButton->getSelectCallbacks().add(this,&LidarViewer::extractCylinderCallback);
	
	GLMotif::Button* intersectPrimitivesButton=new GLMotif::Button("IntersectPrimitivesButton",extractionMenu,"Intersect Primitives");
	intersectPrimitivesButton->getSelectCallbacks().add(this,&LidarViewer::intersectPrimitivesCallback);
	
	GLMotif::Button* loadPrimitivesButton=new GLMotif::Button("LoadPrimitivesButton",extractionMenu,"Load Primitives");
	loadPrimitivesButton->getSelectCallbacks().add(this,&LidarViewer::loadPrimitivesCallback);
	
	GLMotif::Button* savePrimitivesButton=new GLMotif::Button("SavePrimitivesButton",extractionMenu,"Save Primitives");
	savePrimitivesButton->getSelectCallbacks().add(this,&LidarViewer::savePrimitivesCallback);
	
	GLMotif::Button* deleteSelectedPrimitivesButton=new GLMotif::Button("DeleteSelectedPrimitivesButton",extractionMenu,"Delete Selected Primitives");
	deleteSelectedPrimitivesButton->getSelectCallbacks().add(this,&LidarViewer::deleteSelectedPrimitivesCallback);
	
	GLMotif::Button* clearPrimitivesButton=new GLMotif::Button("ClearPrimitivesButton",extractionMenu,"Clear Primitives");
	clearPrimitivesButton->getSelectCallbacks().add(this,&LidarViewer::clearPrimitivesCallback);
	
	extractionMenu->manageChild();
	
	return extractionMenuPopup;
	}

GLMotif::Popup* LidarViewer::createViewMenu(void)
	{
	GLMotif::Popup* viewMenuPopup=new GLMotif::Popup("ViewMenuPopup",Vrui::getWidgetManager());
	
	GLMotif::SubMenu* viewMenu=new GLMotif::SubMenu("ViewMenu",viewMenuPopup,false);
	
	GLMotif::Button* centerDisplayButton=new GLMotif::Button("CenterDisplayButton",viewMenu,"Center Display");
	centerDisplayButton->getSelectCallbacks().add(this,&LidarViewer::centerDisplayCallback);
	
	GLMotif::Button* saveViewpointButton=new GLMotif::Button("SaveViewpointButton",viewMenu,"Save Viewpoint");
	saveViewpointButton->getSelectCallbacks().add(this,&LidarViewer::saveViewpointCallback);
	
	GLMotif::Button* loadViewpointButton=new GLMotif::Button("LoadViewpointButton",viewMenu,"Load Viewpoint");
	loadViewpointButton->getSelectCallbacks().add(this,&LidarViewer::loadViewpointCallback);
	
	viewMenu->manageChild();
	
	return viewMenuPopup;
	}

GLMotif::Popup* LidarViewer::createDialogMenu(void)
	{
	GLMotif::Popup* dialogMenuPopup=new GLMotif::Popup("DialogMenuPopup",Vrui::getWidgetManager());
	
	GLMotif::SubMenu* dialogMenu=new GLMotif::SubMenu("DialogMenu",dialogMenuPopup,false);
	
	GLMotif::ToggleButton* showRenderDialogToggle=new GLMotif::ToggleButton("ShowRenderDialogToggle",dialogMenu,"Show Render Dialog");
	showRenderDialogToggle->setToggle(false);
	showRenderDialogToggle->getValueChangedCallbacks().add(this,&LidarViewer::showRenderDialogCallback);
	
	GLMotif::ToggleButton* showInteractionDialogToggle=new GLMotif::ToggleButton("ShowInteractionDialogToggle",dialogMenu,"Show Interaction Dialog");
	showInteractionDialogToggle->setToggle(false);
	showInteractionDialogToggle->getValueChangedCallbacks().add(this,&LidarViewer::showInteractionDialogCallback);
	
	dialogMenu->manageChild();
	
	return dialogMenuPopup;
	}

GLMotif::PopupMenu* LidarViewer::createMainMenu(void)
	{
	GLMotif::PopupMenu* mainMenuPopup=new GLMotif::PopupMenu("MainMenuPopup",Vrui::getWidgetManager());
	mainMenuPopup->setTitle("LiDAR Viewer");
	
	GLMotif::Menu* mainMenu=new GLMotif::Menu("MainMenu",mainMenuPopup,false);
	
	GLMotif::CascadeButton* selectorModesCascade=new GLMotif::CascadeButton("SelectorModesCascade",mainMenu,"Selector Modes");
	selectorModesCascade->setPopup(createSelectorModesMenu());
	
	GLMotif::CascadeButton* selectionCascade=new GLMotif::CascadeButton("SelectionCascade",mainMenu,"Selection");
	selectionCascade->setPopup(createSelectionMenu());
	
	GLMotif::CascadeButton* extractionCascade=new GLMotif::CascadeButton("ExtractionCascade",mainMenu,"Primitives");
	extractionCascade->setPopup(createExtractionMenu());
	
	GLMotif::CascadeButton* viewCascade=new GLMotif::CascadeButton("ViewCascade",mainMenu,"View");
	viewCascade->setPopup(createViewMenu());
	
	GLMotif::CascadeButton* dialogCascade=new GLMotif::CascadeButton("DialogCascade",mainMenu,"Dialogs");
	dialogCascade->setPopup(createDialogMenu());
	
	#if 0
	GLMotif::ToggleButton* updateTreeToggle=new GLMotif::ToggleButton("UpdateTreeToggle",mainMenu,"Update Tree");
	updateTreeToggle->setToggle(updateTree);
	updateTreeToggle->getValueChangedCallbacks().add(this,&LidarViewer::updateTreeCallback);
	#endif
	
	mainMenu->manageChild();
	
	return mainMenuPopup;
	}

GLMotif::PopupWindow* LidarViewer::createRenderDialog(void)
	{
	const GLMotif::StyleSheet& ss=*Vrui::getWidgetManager()->getStyleSheet();
	
	GLMotif::PopupWindow* renderDialog=new GLMotif::PopupWindow("RenderDialog",Vrui::getWidgetManager(),"Render Settings");
	
	GLMotif::RowColumn* renderSettings=new GLMotif::RowColumn("RenderSettings",renderDialog,false);
	renderSettings->setNumMinorWidgets(3);
	
	/* Create a slider/textfield combo to change the rendering quality: */
	GLMotif::Label* renderQualityLabel=new GLMotif::Label("RenderQualityLabel",renderSettings,"Render Quality");
	
	renderQualityValue=new GLMotif::TextField("RenderQualityValue",renderSettings,5);
	renderQualityValue->setFieldWidth(5);
	renderQualityValue->setPrecision(2);
	renderQualityValue->setValue(double(renderQuality));
	
	renderQualitySlider=new GLMotif::Slider("RenderQualitySlider",renderSettings,GLMotif::Slider::HORIZONTAL,ss.fontHeight*10.0f);
	renderQualitySlider->setValueRange(-3.0,3.0,0.01);
	renderQualitySlider->setValue(double(renderQuality));
	renderQualitySlider->getValueChangedCallbacks().add(this,&LidarViewer::renderQualitySliderCallback);
	
	/* Create a slider/textfield combo to change the focus+context weight: */
	GLMotif::Label* fncWeightLabel=new GLMotif::Label("FncWeightLabel",renderSettings,"Focus + Context");
	
	fncWeightValue=new GLMotif::TextField("FncWeightValue",renderSettings,5);
	fncWeightValue->setFieldWidth(5);
	fncWeightValue->setPrecision(2);
	fncWeightValue->setValue(double(fncWeight));
	
	fncWeightSlider=new GLMotif::Slider("FncWeightSlider",renderSettings,GLMotif::Slider::HORIZONTAL,ss.fontHeight*10.0f);
	fncWeightSlider->setValueRange(0.0,2.0,0.01);
	fncWeightSlider->setValue(double(fncWeight));
	fncWeightSlider->getValueChangedCallbacks().add(this,&LidarViewer::fncWeightSliderCallback);
	
	/* Create a slider/textfield combo to change the point size: */
	GLMotif::Label* pointSizeLabel=new GLMotif::Label("PointSizeLabel",renderSettings,"Point Size");
	
	pointSizeValue=new GLMotif::TextField("PointSizeValue",renderSettings,3);
	pointSizeValue->setFieldWidth(3);
	pointSizeValue->setPrecision(1);
	pointSizeValue->setValue(double(pointSize));
	
	pointSizeSlider=new GLMotif::Slider("PointSizeSlider",renderSettings,GLMotif::Slider::HORIZONTAL,ss.fontHeight*10.0f);
	pointSizeSlider->setValueRange(1.0,10.0,0.5);
	pointSizeSlider->setValue(double(pointSize));
	pointSizeSlider->getValueChangedCallbacks().add(this,&LidarViewer::pointSizeSliderCallback);
	
	/* Create a toggle button to enable plane distance visualization: */
	GLMotif::ToggleButton* enableTextureToggle=new GLMotif::ToggleButton("EnableTextureToggle",renderSettings,"Show Plane Distance");
	enableTextureToggle->setToggle(useTexturePlane);
	enableTextureToggle->getValueChangedCallbacks().add(this,&LidarViewer::enableTexturePlaneCallback);

	texturePlaneScaleValue=new GLMotif::TextField("TexturePlaneScaleValue",renderSettings,8);
	texturePlaneScaleValue->setFieldWidth(8);
	texturePlaneScaleValue->setPrecision(3);
	texturePlaneScaleValue->setValue(double(texturePlaneScale));
	
	texturePlaneScaleSlider=new GLMotif::Slider("TexturePlaneScaleSlider",renderSettings,GLMotif::Slider::HORIZONTAL,ss.fontHeight*10.0f);
	texturePlaneScaleSlider->setValueRange(-1.0,4.0,0.1);
	texturePlaneScaleSlider->setValue(Math::log10(double(texturePlaneScale)));
	texturePlaneScaleSlider->getValueChangedCallbacks().add(this,&LidarViewer::texturePlaneScaleSliderCallback);
	
	renderSettings->manageChild();
	
	return renderDialog;
	}

GLMotif::PopupWindow* LidarViewer::createInteractionDialog(void)
	{
	const GLMotif::StyleSheet& ss=*Vrui::getWidgetManager()->getStyleSheet();
	
	GLMotif::PopupWindow* interactionDialog=new GLMotif::PopupWindow("InteractionDialog",Vrui::getWidgetManager(),"Interaction Settings");
	
	GLMotif::RowColumn* interactionSettings=new GLMotif::RowColumn("InteractionSettings",interactionDialog,false);
	
	GLMotif::ToggleButton* overrideToolsToggle=new GLMotif::ToggleButton("OverrideToolsToggle",interactionSettings,"Override Tools");
	overrideToolsToggle->setToggle(overrideTools);
	overrideToolsToggle->getValueChangedCallbacks().add(this,&LidarViewer::overrideToolsCallback);
	
	/* Create a slider/textfield combo to change the interaction sphere size: */
	GLMotif::RowColumn* brushSizeBox=new GLMotif::RowColumn("BrushSize",interactionSettings,false);
	brushSizeBox->setOrientation(GLMotif::RowColumn::HORIZONTAL);
	
	GLMotif::Label* brushSizeLabel=new GLMotif::Label("BrushSizeLabel",brushSizeBox,"Brush Size");
	
	brushSizeValue=new GLMotif::TextField("BrushSizeValue",brushSizeBox,7);
	brushSizeValue->setFieldWidth(7);
	brushSizeValue->setPrecision(3);
	brushSizeValue->setValue(double(brushSize));
	
	GLMotif::Slider* brushSizeSlider=new GLMotif::Slider("BrushSizeSlider",brushSizeBox,GLMotif::Slider::HORIZONTAL,ss.fontHeight*10.0f);
	brushSizeSlider->setValueRange(double(brushSize)*0.1,double(brushSize)*5.0,double(brushSize)*0.01);
	brushSizeSlider->setValue(double(brushSize));
	brushSizeSlider->getValueChangedCallbacks().add(this,&LidarViewer::brushSizeSliderCallback);
	
	brushSizeBox->manageChild();
	
	/* Create a radio box to select interaction modes: */
	GLMotif::RadioBox* selectorModes=new GLMotif::RadioBox("SelectorModes",interactionSettings,false);
	selectorModes->setOrientation(GLMotif::RowColumn::HORIZONTAL);
	selectorModes->setPacking(GLMotif::RowColumn::PACK_GRID);
	selectorModes->setSelectionMode(GLMotif::RadioBox::ALWAYS_ONE);
	
	selectorModes->addToggle("Add");
	selectorModes->addToggle("Subtract");
	
	switch(defaultSelectorMode)
		{
		case SelectorLocator::ADD:
			selectorModes->setSelectedToggle(0);
			break;
		
		case SelectorLocator::SUBTRACT:
			selectorModes->setSelectedToggle(1);
			break;
		}
	selectorModes->getValueChangedCallbacks().add(this,&LidarViewer::changeSelectorModeCallback);
	selectorModes->manageChild();
	
	interactionSettings->manageChild();
	
	return interactionDialog;
	}

void LidarViewer::draggingToolCallback(Vrui::DraggingTool::DragStartCallbackData* cbData)
	{
	/* Get the dragging tool's selection point in model coordinates: */
	Point pickPoint=cbData->startTransformation.getOrigin();
	
	/* Pick against all primitives: */
	Scalar minDistance=Vrui::getInchFactor()*Vrui::getInverseNavigationTransformation().getScaling();
	int pickedPrimitiveIndex=-1;
	for(int i=0;i<int(primitives.size());++i)
		{
		Scalar dist=primitives[i]->pick(pickPoint,minDistance);
		if(minDistance>dist)
			{
			minDistance=dist;
			pickedPrimitiveIndex=i;
			}
		}
	
	if(pickedPrimitiveIndex>=0)
		{
		/* Select or deselect the primitive: */
		if(primitiveSelectedFlags[pickedPrimitiveIndex])
			deselectPrimitive(pickedPrimitiveIndex);
		else
			selectPrimitive(pickedPrimitiveIndex);
		
		lastPickedPrimitive=pickedPrimitiveIndex;
		}
	}

int LidarViewer::addPrimitive(Primitive* newPrimitive)
	{
	primitives.push_back(newPrimitive);
	primitiveSelectedFlags.push_back(false);
	return int(primitives.size())-1;
	}

void LidarViewer::selectPrimitive(int primitiveIndex)
	{
	if(!primitiveSelectedFlags[primitiveIndex])
		{
		primitives[primitiveIndex]->setSurfaceColor(Primitive::Color(0.1f,0.5f,0.5f,0.5f));
		primitiveSelectedFlags[primitiveIndex]=true;
		}
	}

void LidarViewer::deselectPrimitive(int primitiveIndex)
	{
	if(primitiveSelectedFlags[primitiveIndex])
		{
		primitives[primitiveIndex]->setSurfaceColor(Primitive::Color(0.5f,0.5f,0.1f,0.5f));
		primitiveSelectedFlags[primitiveIndex]=false;
		}
	}

void LidarViewer::deletePrimitive(int primitiveIndex)
	{
	delete primitives[primitiveIndex];
	primitives.erase(primitives.begin()+primitiveIndex);
	primitiveSelectedFlags.erase(primitiveSelectedFlags.begin()+primitiveIndex);
	}

LidarViewer::LidarViewer(int& argc,char**& argv,char**& appDefaults)
	:Vrui::Application(argc,argv,appDefaults),
	 octree(0),renderQuality(0),fncWeight(0.5),pointSize(3.0f),
	 useTexturePlane(false),texturePlane(GPlane::Vector(0.0,0.0,1.0),0.0),texturePlaneScale(100.0),
	 updateTree(true),
	 lastFrameTime(Vrui::getApplicationTime()),
	 overrideTools(true),
	 brushSize(Vrui::getGlyphRenderer()->getGlyphSize()*Vrui::Scalar(2.5)),
	 defaultSelectorMode(SelectorLocator::ADD),
	 extractorPipe(Vrui::openPipe()),
	 lastPickedPrimitive(-1),
	 mainMenu(0),renderDialog(0)
	{
	/* Parse the command line: */
	unsigned int memCacheSize=512*1024*1024;
	unsigned int gfxCacheSize=160*1024*1024;
	bool haveOctreeFile=false;
	const char* viewpointFileName=0;
	for(int i=1;i<argc;++i)
		{
		if(argv[i][0]=='-')
			{
			if(strcasecmp(argv[i]+1,"cache")==0)
				{
				if(i+1<argc)
					{
					++i;
					memCacheSize=atoi(argv[i])*1024*1024;
					}
				}
			else if(strcasecmp(argv[i]+1,"gfxcache")==0)
				{
				if(i+1<argc)
					{
					++i;
					gfxCacheSize=atoi(argv[i])*1024*1024;
					}
				}
			else if(strcasecmp(argv[i]+1,"pointsize")==0)
				{
				if(i+1<argc)
					{
					++i;
					pointSize=float(atof(argv[i]));
					}
				}
			else if(strcasecmp(argv[i]+1,"loadview")==0)
				{
				if(i+1<argc)
					{
					++i;
					viewpointFileName=argv[i];
					}
				}
			else if(strcasecmp(argv[i]+1,"quality")==0)
				{
				if(i+1<argc)
					{
					++i;
					renderQuality=Scalar(atof(argv[i]));
					}
				}
			}
		else if(!haveOctreeFile)
			{
			/* Load a LiDAR octree file: */
			octree=new LidarOctree(argv[i],memCacheSize,gfxCacheSize);
			
			haveOctreeFile=true;
			}
		else
			{
			/* Get the octree's point offset: */
			typedef Geometry::TranslationTransformation<double,3> ReadTransform;
			ReadTransform t=ReadTransform::translate(-octree->getPointOffset());
			
			/* Load an ARC/INFO export file: */
			arcInfoExportFiles.push_back(new ArcInfoExportFile(argv[i]));
			}
		}
	
	if(!haveOctreeFile)
		Misc::throwStdErr("LidarViewer::LidarViewer: No octree file name provided");
	
	/* Initialize the octree: */
	octree->setRenderQuality(renderQuality);
	octree->setTreeUpdateFunction(treeUpdateNotificationCB,0);
	
	/* Register a coordinate transform object to undo the coordinate offset done by the octree object: */
	Vrui::getCoordinateManager()->setCoordinateTransform(new Vrui::OrthogonalCoordinateTransform(Vrui::OGTransform::translate(Vrui::Vector(-octree->getPointOffset()))));
	
	/* Create the GUI: */
	mainMenu=createMainMenu();
	Vrui::setMainMenu(mainMenu);
	renderDialog=createRenderDialog();
	interactionDialog=createInteractionDialog();
	
	/* Initialize the navigation transformation: */
	if(viewpointFileName!=0)
		{
		/* Open viewpoint file: */
		Misc::File viewpointFile(viewpointFileName,"rb",Misc::File::LittleEndian);
		
		/* Read the navigation transformation: */
		Vrui::NavTransform::Vector translation;
		viewpointFile.read(translation.getComponents(),3);
		Vrui::NavTransform::Rotation::Scalar quaternion[4];
		viewpointFile.read(quaternion,4);
		Vrui::NavTransform::Scalar scaling=viewpointFile.read<Vrui::NavTransform::Scalar>();
		
		/* Set the navigation transformation: */
		Vrui::setNavigationTransformation(Vrui::NavTransform(translation,Vrui::NavTransform::Rotation::fromQuaternion(quaternion),scaling));
		}
	else
		centerDisplayCallback(0);
	
	/* Register the custom tool class with the Vrui tool manager: */
	LidarToolFactory* lidarToolFactory=new LidarToolFactory(*Vrui::getToolManager(),octree);
	Vrui::getToolManager()->addClass(lidarToolFactory,LidarToolFactory::factoryDestructor);
	}

LidarViewer::~LidarViewer(void)
	{
	/* Destroy all locators: */
	for(LocatorList::iterator lIt=locators.begin();lIt!=locators.end();++lIt)
		delete *lIt;
	
	/* Close the synchronization pipe: */
	delete extractorPipe;
	
	/* Destroy all primitives: */
	for(PrimitiveList::iterator pIt=primitives.begin();pIt!=primitives.end();++pIt)
		delete* pIt;
	
	/* Delete the GUI: */
	delete mainMenu;
	delete renderDialog;
	delete interactionDialog;
	
	/* Delete the octree: */
	delete octree;
	
	/* Delete all ARC/INFO export files: */
	for(std::vector<ArcInfoExportFile*>::iterator aIt=arcInfoExportFiles.begin();aIt!=arcInfoExportFiles.end();++aIt)
		delete *aIt;
	}

void LidarViewer::initContext(GLContextData& contextData) const
	{
	/* Create a new context entry: */
	DataItem* dataItem=new DataItem;
	contextData.addDataItem(this,dataItem);
	
	/* Create the influence sphere display list: */
	glNewList(dataItem->influenceSphereDisplayListId,GL_COMPILE);
	glDisable(GL_CULL_FACE);
	glEnable(GL_BLEND);
	glDepthMask(GL_FALSE);
	glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
	glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
	glColor4f(0.6f,0.6f,0.1f,0.5f); // Yellow sphere
	// glColor4f(0.6f,0.1f,0.6f,0.5f); // Purple sphere
	glDrawSphereIcosahedron(1.0,5);
	glBlendFunc(GL_ONE,GL_ONE);
	glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
	glLineWidth(1.0f);
	glColor3f(0.025f,0.025f,0.025f);
	glDrawSphereIcosahedron(1.0,5);
	glDepthMask(GL_TRUE);
	glDisable(GL_BLEND);
	glEndList();
	
	/* Create the texture plane color map: */
	GLColorMap colorMap(GLColorMap::RAINBOW,1.0f,1.0f,0.0,1.0);
	
	/* Create the color map texture image: */
	glBindTexture(GL_TEXTURE_1D,dataItem->planeColorMapTextureId);
	glTexParameteri(GL_TEXTURE_1D,GL_TEXTURE_MIN_FILTER,GL_LINEAR);
	glTexParameteri(GL_TEXTURE_1D,GL_TEXTURE_MAG_FILTER,GL_LINEAR);
	glTexParameteri(GL_TEXTURE_1D,GL_TEXTURE_BASE_LEVEL,0);
	glTexParameteri(GL_TEXTURE_1D,GL_TEXTURE_MAX_LEVEL,0);
	glTexParameteri(GL_TEXTURE_1D,GL_TEXTURE_WRAP_S,GL_CLAMP_TO_EDGE);
	glTexImage1D(GL_TEXTURE_1D,0,GL_RGB,colorMap.getNumEntries(),0,GL_RGBA,GL_FLOAT,colorMap.getColors());
	glBindTexture(GL_TEXTURE_1D,0);
	}

void LidarViewer::toolCreationCallback(Vrui::ToolManager::ToolCreationCallbackData* cbData)
	{
	/* Check if the new tool is a locator tool: */
	Vrui::LocatorTool* ltool=dynamic_cast<Vrui::LocatorTool*>(cbData->tool);
	if(ltool!=0)
		{
		/* Create new locator: */
		Locator* newLocator=new SelectorLocator(ltool,this);
		
		/* Add new locator to list: */
		locators.push_back(newLocator);
		}
	
	/* Check if the new tool is a dragging tool: */
	Vrui::DraggingTool* dtool=dynamic_cast<Vrui::DraggingTool*>(cbData->tool);
	if(dtool!=0)
		{
		/* Register a dragging tool callback with the new dragging tool: */
		dtool->getDragStartCallbacks().add(this,&LidarViewer::draggingToolCallback);
		}
	}

void LidarViewer::toolDestructionCallback(Vrui::ToolManager::ToolDestructionCallbackData* cbData)
	{
	/* Check if the to-be-destroyed tool is a locator tool: */
	Vrui::LocatorTool* ltool=dynamic_cast<Vrui::LocatorTool*>(cbData->tool);
	if(ltool!=0)
		{
		/* Find the locator associated with the tool in the list: */
		LocatorList::iterator lIt;
		for(lIt=locators.begin();lIt!=locators.end()&&(*lIt)->getTool()!=ltool;++lIt)
			;
		if(lIt!=locators.end())
			{
			/* Remove the locator: */
			delete *lIt;
			locators.erase(lIt);
			}
		}
	}

void LidarViewer::frame(void)
	{
	#if 0
	/* Get the current frame time and adapt the rendering quality: */
	double fr=Vrui::getCurrentFrameTime();
	if(fr>0.0)
		{
		fr=1.0/fr;
		if(fr<45.0)
			{
			renderQuality-=0.01f*(45.0-fr);
			octree->setRenderQuality(renderQuality);
			renderQualitySlider->setValue(double(renderQuality));
			}
		else if(fr>55.0)
			{
			renderQuality+=0.01f;
			octree->setRenderQuality(renderQuality);
			renderQualitySlider->setValue(double(renderQuality));
			}
		}
	#endif
	
	/* Prepare the next rendering pass: */
	octree->startRenderPass();
	
	/* Prepare for interaction with all interactors: */
	for(LocatorList::iterator lIt=locators.begin();lIt!=locators.end();++lIt)
		{
		SelectorLocator* sl=dynamic_cast<SelectorLocator*>(*lIt);
		if(sl!=0&&sl->ready)
			octree->interact(LidarOctree::Interactor(sl->modelCenter,sl->modelRadius));
		}
	}

void LidarViewer::display(GLContextData& contextData) const
	{
	/* Retrieve context entry: */
	DataItem* dataItem=contextData.retrieveDataItem<DataItem>(this);
	
	/* Set up basic OpenGL state: */
	glPushAttrib(GL_ENABLE_BIT|GL_LINE_BIT|GL_POINT_BIT|GL_TEXTURE_BIT);
	glDisable(GL_LIGHTING);
	glPointSize(pointSize);
	
	if(useTexturePlane)
		{
		/* Set up automatic texture coordinate generation: */
		glTexGeni(GL_S,GL_TEXTURE_GEN_MODE,GL_OBJECT_LINEAR);
		GLdouble planeCoeff[4];
		for(int i=0;i<3;++i)
			planeCoeff[i]=texturePlane.getNormal()[i]/texturePlaneScale;
		planeCoeff[3]=0.5-texturePlane.getOffset()/texturePlaneScale;
		glTexGendv(GL_S,GL_OBJECT_PLANE,planeCoeff);
		glEnable(GL_TEXTURE_GEN_S);
		
		/* Enable 1D texture mapping: */
		glEnable(GL_TEXTURE_1D);
		glDisable(GL_TEXTURE_2D);
		glDisable(GL_TEXTURE_3D);
		glBindTexture(GL_TEXTURE_1D,dataItem->planeColorMapTextureId);
		glTexEnvi(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,GL_REPLACE);
		}
	
	/* Render the LiDAR point tree: */
	Point displayCenter=Point(Vrui::getInverseNavigationTransformation().transform(Vrui::getDisplayCenter()));
	Scalar displaySize=Scalar(Vrui::getInverseNavigationTransformation().getScaling()*Vrui::getDisplaySize());
	octree->setFocusAndContext(displayCenter,displaySize*Scalar(0.5),fncWeight);
	LidarOctree::Frustum frustum;
	frustum.setFromGL();
	octree->glRenderAction(frustum,contextData);
	
	if(useTexturePlane)
		{
		/* Disable 1D texture mapping: */
		glBindTexture(GL_TEXTURE_1D,0);
		glDisable(GL_TEXTURE_1D);
		
		/* Disable automatic texture coordinate generation: */
		glDisable(GL_TEXTURE_GEN_S);
		}
	
	/* Draw all ARC/INFO export files: */
	glLineWidth(1.0f);
	glColor3f(1.0f,1.0f,0.0f);
	for(std::vector<ArcInfoExportFile*>::const_iterator aIt=arcInfoExportFiles.begin();aIt!=arcInfoExportFiles.end();++aIt)
		(*aIt)->glRenderAction(contextData);
	
	/* Reset OpenGL state: */
	glPopAttrib();
	
	/* Render all locators: */
	for(LocatorList::const_iterator lIt=locators.begin();lIt!=locators.end();++lIt)
		(*lIt)->glRenderAction(contextData);
	
	/* Render all extracted primitives: */
	for(PrimitiveList::const_iterator pIt=primitives.begin();pIt!=primitives.end();++pIt)
		(*pIt)->glRenderAction(contextData);
	}

void LidarViewer::centerDisplayCallback(Misc::CallbackData* cbData)
	{
	/* Initialize the navigation transformation: */
	Vrui::setNavigationTransformation(octree->getDomainCenter(),octree->getDomainRadius(),Vrui::Vector(0,0,1));
	}

void LidarViewer::changeSelectorModeCallback(GLMotif::RadioBox::ValueChangedCallbackData* cbData)
	{
	switch(cbData->radioBox->getToggleIndex(cbData->newSelectedToggle))
		{
		case 0:
			defaultSelectorMode=SelectorLocator::ADD;
			break;
		
		case 1:
			defaultSelectorMode=SelectorLocator::SUBTRACT;
			break;
		}
	
	if(overrideTools)
		{
		/* Apply current tool settings to all existing tools: */
		for(LocatorList::iterator lIt=locators.begin();lIt!=locators.end();++lIt)
			(*lIt)->updateSettings();
		}
	}

void LidarViewer::saveSelectionCallback(Misc::CallbackData* cbData)
	{
	if(Vrui::isMaster())
		{
		/* Create a selection saver functor: */
		LidarSelectionSaver lss("SelectedPoints.xyzrgb",octree->getPointOffset());
		
		/* Save all selected points: */
		octree->processSelectedPoints(lss);
		}
	}

void LidarViewer::clearSelectionCallback(Misc::CallbackData* cbData)
	{
	/* Clear the point selection: */
	octree->clearSelection();
	}

void LidarViewer::extractPlaneCallback(Misc::CallbackData* cbData)
	{
	try
		{
		PlanePrimitive* primitive;
		if(Vrui::isMaster())
			primitive=new PlanePrimitive(octree,extractorPipe);
		else
			primitive=new PlanePrimitive(extractorPipe);
		
		/* Store the primitive: */
		lastPickedPrimitive=addPrimitive(primitive);
		
		/* Update the texture generation plane: */
		texturePlane=primitive->getPlane();
		texturePlane.normalize();
		}
	catch(std::runtime_error err)
		{
		if(Vrui::isMaster())
			std::cerr<<"Did not extract plane due to exception "<<err.what()<<std::endl;
		}
	}

void LidarViewer::extractSphereCallback(Misc::CallbackData* cbData)
	{
	try
		{
		SpherePrimitive* primitive;
		if(Vrui::isMaster())
			primitive=new SpherePrimitive(octree,extractorPipe);
		else
			primitive=new SpherePrimitive(extractorPipe);
		
		/* Store the primitive: */
		lastPickedPrimitive=addPrimitive(primitive);
		}
	catch(std::runtime_error err)
		{
		if(Vrui::isMaster())
			std::cerr<<"Did not extract sphere due to exception "<<err.what()<<std::endl;
		}
	}

void LidarViewer::extractCylinderCallback(Misc::CallbackData* cbData)
	{
	try
		{
		CylinderPrimitive* primitive;
		if(Vrui::isMaster())
			primitive=new CylinderPrimitive(octree,extractorPipe);
		else
			primitive=new CylinderPrimitive(extractorPipe);
		
		/* Store the primitive: */
		lastPickedPrimitive=addPrimitive(primitive);
		}
	catch(std::runtime_error err)
		{
		if(Vrui::isMaster())
			std::cerr<<"Did not extract cylinder due to exception "<<err.what()<<std::endl;
		}
	}

void LidarViewer::intersectPrimitivesCallback(Misc::CallbackData* cbData)
	{
	try
		{
		/* Get all selected primitives: */
		std::vector<PlanePrimitive*> planes;
		std::vector<LinePrimitive*> lines;
		std::vector<PointPrimitive*> points;
		for(unsigned int i=0;i<primitives.size();++i)
			if(primitiveSelectedFlags[i])
				{
				if(dynamic_cast<PlanePrimitive*>(primitives[i])!=0)
					planes.push_back(dynamic_cast<PlanePrimitive*>(primitives[i]));
				else if(dynamic_cast<LinePrimitive*>(primitives[i])!=0)
					lines.push_back(dynamic_cast<LinePrimitive*>(primitives[i]));
				else if(dynamic_cast<PointPrimitive*>(primitives[i])!=0)
					points.push_back(dynamic_cast<PointPrimitive*>(primitives[i]));
				}
		
		/* Determine what to do based on the selection set: */
		Primitive* primitive;
		if(planes.size()==2&&lines.empty()&&points.empty())
			{
			/* Create a line by intersecting two planes: */
			if(Vrui::isMaster())
				primitive=new LinePrimitive(planes[0],planes[1],extractorPipe);
			else
				primitive=new LinePrimitive(extractorPipe);
			}
		else if(planes.size()==3&&lines.empty()&&points.empty())
			{
			/* Create a point by intersecting three planes: */
			if(Vrui::isMaster())
				primitive=new PointPrimitive(planes[0],planes[1],planes[2],extractorPipe);
			else
				primitive=new PointPrimitive(extractorPipe);
			}
		else if(planes.size()==1&&lines.size()==1&&points.empty())
			{
			/* Create a point by intersecting a plane and a line: */
			if(Vrui::isMaster())
				primitive=new PointPrimitive(planes[0],lines[0],extractorPipe);
			else
				primitive=new PointPrimitive(extractorPipe);
			}
		else
			Misc::throwStdErr("mismatching selected primitives");
		
		/* Store the primitive: */
		lastPickedPrimitive=addPrimitive(primitive);
		
		/* Unselect all primitives: */
		for(int i=0;i<int(primitives.size());++i)
			deselectPrimitive(i);
		}
	catch(std::runtime_error err)
		{
		if(Vrui::isMaster())
			std::cerr<<"Did not intersect primitives due to exception "<<err.what()<<std::endl;
		}
	}

void LidarViewer::loadPrimitivesCallback(Misc::CallbackData* cbData)
	{
	/* Open the primitive file: */
	Misc::File primitiveFile("Primitives.dat","rb",Misc::File::LittleEndian);
	
	/* Read the file header: */
	char header[40];
	primitiveFile.read<char>(header,sizeof(header));
	if(strcmp(header,"LidarViewer primitive file v1.2       \n")!=0)
		Misc::throwStdErr("LidarViewer::loadPrimitivesCallback: File %s is not a valid version 1.2 primitive file","Primitives.dat");
	
	/* Read all primitives in the file: */
	while(true)
		{
		/* Read the primitive type: */
		int primitiveType;
		try
			{
			primitiveType=primitiveFile.read<int>();
			}
		catch(Misc::File::ReadError err)
			{
			/* Stop reading: */
			break;
			}
		
		/* Create a primitive of the appropriate type: */
		Primitive* newPrimitive=0;
		switch(primitiveType)
			{
			case 0:
				newPrimitive=new PlanePrimitive(primitiveFile,Primitive::Vector(-octree->getPointOffset()));
				break;
			
			case 1:
				newPrimitive=new SpherePrimitive(primitiveFile,Primitive::Vector(-octree->getPointOffset()));
				break;
			
			case 2:
				newPrimitive=new CylinderPrimitive(primitiveFile,Primitive::Vector(-octree->getPointOffset()));
				break;
			
			case 3:
				newPrimitive=new LinePrimitive(primitiveFile,Primitive::Vector(-octree->getPointOffset()));
				break;
			
			case 4:
				newPrimitive=new PointPrimitive(primitiveFile,Primitive::Vector(-octree->getPointOffset()));
				break;
			
			default:
				Misc::throwStdErr("LidarViewer::loadPrimitivesCallback: Unknown primitive type %d",primitiveType);
			}
		
		/* Store the primitive: */
		lastPickedPrimitive=addPrimitive(newPrimitive);
		}
	}

void LidarViewer::savePrimitivesCallback(Misc::CallbackData* cbData)
	{
	/* Open the primitive file: */
	Misc::File primitiveFile("Primitives.dat","wb",Misc::File::LittleEndian);
	
	/* Write the file header: */
	char header[40];
	snprintf(header,sizeof(header),"LidarViewer primitive file v1.2       \n");
	primitiveFile.write<char>(header,sizeof(header));
	
	/* Write all primitives to the file: */
	for(std::vector<Primitive*>::const_iterator pIt=primitives.begin();pIt!=primitives.end();++pIt)
		{
		/* Determine and write the primitive type: */
		if(dynamic_cast<const PlanePrimitive*>(*pIt)!=0)
			primitiveFile.write<int>(0);
		else if(dynamic_cast<const SpherePrimitive*>(*pIt)!=0)
			primitiveFile.write<int>(1);
		else if(dynamic_cast<const CylinderPrimitive*>(*pIt)!=0)
			primitiveFile.write<int>(2);
		else if(dynamic_cast<const LinePrimitive*>(*pIt)!=0)
			primitiveFile.write<int>(3);
		else if(dynamic_cast<const PointPrimitive*>(*pIt)!=0)
			primitiveFile.write<int>(4);
		
		/* Write the primitive: */
		(*pIt)->write(primitiveFile,Primitive::Vector(octree->getPointOffset()));
		}
	}

void LidarViewer::deleteSelectedPrimitivesCallback(Misc::CallbackData* cbData)
	{
	for(int i=primitives.size()-1;i>=0;--i)
		if(primitiveSelectedFlags[i])
			deletePrimitive(i);
	lastPickedPrimitive=-1;
	}

void LidarViewer::clearPrimitivesCallback(Misc::CallbackData* cbData)
	{
	/* Destroy all primitives: */
	for(PrimitiveList::iterator pIt=primitives.begin();pIt!=primitives.end();++pIt)
		delete* pIt;
	primitives.clear();
	primitiveSelectedFlags.clear();
	lastPickedPrimitive=-1;
	}

void LidarViewer::showRenderDialogCallback(GLMotif::ToggleButton::ValueChangedCallbackData* cbData)
	{
	if(cbData->set)
		{
		/* Open the render dialog at the same position as the main menu: */
		Vrui::getWidgetManager()->popupPrimaryWidget(renderDialog,Vrui::getWidgetManager()->calcWidgetTransformation(mainMenu));
		}
	else
		{
		/* Close the render dialog: */
		Vrui::popdownPrimaryWidget(renderDialog);
		}
	}

void LidarViewer::renderQualitySliderCallback(GLMotif::Slider::ValueChangedCallbackData* cbData)
	{
	/* Get the new render quality: */
	renderQuality=Scalar(cbData->value);
	
	/* Update the render quality value label: */
	renderQualityValue->setValue(double(cbData->value));
	
	/* Set the octree's render quality: */
	octree->setRenderQuality(renderQuality);
	}

void LidarViewer::fncWeightSliderCallback(GLMotif::Slider::ValueChangedCallbackData* cbData)
	{
	/* Get the new focus+context weight: */
	fncWeight=Scalar(cbData->value);
	
	/* Update the focus+context weight value label: */
	fncWeightValue->setValue(double(cbData->value));
	}

void LidarViewer::pointSizeSliderCallback(GLMotif::Slider::ValueChangedCallbackData* cbData)
	{
	/* Get the new point size: */
	pointSize=cbData->value;
	
	/* Update the point size value label: */
	pointSizeValue->setValue(double(cbData->value));
	}

void LidarViewer::enableTexturePlaneCallback(GLMotif::ToggleButton::ValueChangedCallbackData* cbData)
	{
	useTexturePlane=cbData->set;
	}

void LidarViewer::texturePlaneScaleSliderCallback(GLMotif::Slider::ValueChangedCallbackData* cbData)
	{
	/* Get the new texture plane scaling factor: */
	texturePlaneScale=Math::pow(10.0,double(cbData->value));
	
	/* Update the texture plane scale value label: */
	texturePlaneScaleValue->setValue(double(texturePlaneScale));
	}

void LidarViewer::showInteractionDialogCallback(GLMotif::ToggleButton::ValueChangedCallbackData* cbData)
	{
	if(cbData->set)
		{
		/* Open the interaction dialog at the same position as the main menu: */
		Vrui::getWidgetManager()->popupPrimaryWidget(interactionDialog,Vrui::getWidgetManager()->calcWidgetTransformation(mainMenu));
		}
	else
		{
		/* Close the interaction dialog: */
		Vrui::popdownPrimaryWidget(interactionDialog);
		}
	}

void LidarViewer::overrideToolsCallback(GLMotif::ToggleButton::ValueChangedCallbackData* cbData)
	{
	overrideTools=cbData->set;
	if(overrideTools)
		{
		/* Apply current tool settings to all existing tools: */
		for(LocatorList::iterator lIt=locators.begin();lIt!=locators.end();++lIt)
			(*lIt)->updateSettings();
		}
	}

void LidarViewer::brushSizeSliderCallback(GLMotif::Slider::ValueChangedCallbackData* cbData)
	{
	/* Get the new brush size: */
	brushSize=cbData->value;
	
	/* Update the brush size value label: */
	brushSizeValue->setValue(double(cbData->value));
	
	if(overrideTools)
		{
		/* Apply current tool settings to all existing tools: */
		for(LocatorList::iterator lIt=locators.begin();lIt!=locators.end();++lIt)
			(*lIt)->updateSettings();
		}
	}

void LidarViewer::loadViewpointCallback(Misc::CallbackData* cbData)
	{
	/* Open viewpoint file: */
	Misc::File viewpointFile("Viewpoint.dat","rb",Misc::File::LittleEndian);
	
	/* Read the navigation transformation: */
	Vrui::NavTransform::Vector translation;
	viewpointFile.read(translation.getComponents(),3);
	Vrui::NavTransform::Rotation::Scalar quaternion[4];
	viewpointFile.read(quaternion,4);
	Vrui::NavTransform::Scalar scaling=viewpointFile.read<Vrui::NavTransform::Scalar>();
	
	/* Set the navigation transformation: */
	Vrui::setNavigationTransformation(Vrui::NavTransform(translation,Vrui::NavTransform::Rotation::fromQuaternion(quaternion),scaling));
	}

void LidarViewer::saveViewpointCallback(Misc::CallbackData* cbData)
	{
	/* Open viewpoint file: */
	Misc::File viewpointFile("Viewpoint.dat","wb",Misc::File::LittleEndian);
	
	/* Get the current navigation transformation: */
	const Vrui::NavTransform& nt=Vrui::getNavigationTransformation();
	
	/* Write the navigation transformation: */
	viewpointFile.write(nt.getTranslation().getComponents(),3);
	viewpointFile.write(nt.getRotation().getQuaternion(),4);
	viewpointFile.write(nt.getScaling());
	}


void LidarViewer::updateTreeCallback(GLMotif::ToggleButton::ValueChangedCallbackData* cbData)
	{
	updateTree=cbData->set;
	}

int main(int argc,char* argv[])
	{
	try
		{
		char** appDefaults=0;
		LidarViewer lv(argc,argv,appDefaults);
		lv.run();
		return 0;
		}
	catch(std::runtime_error err)
		{
		std::cerr<<"Caught exception "<<err.what()<<std::endl;
		return 1;
		}
	}

