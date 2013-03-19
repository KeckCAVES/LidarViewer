/***********************************************************************
LidarViewer - Viewer program for multiresolution LiDAR data.
Copyright (c) 2005-2013 Oliver Kreylos

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

#include "LidarViewer.h"

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <stdexcept>
#include <Misc/SelfDestructPointer.h>
#include <Misc/FunctionCalls.h>
#include <Misc/ThrowStdErr.h>
#include <Misc/PrintInteger.h>
#include <Misc/CreateNumberedFileName.h>
#include <Misc/FileTests.h>
#include <Misc/StringMarshaller.h>
#include <Misc/StandardValueCoders.h>
#include <Misc/ConfigurationFile.h>
#include <IO/File.h>
#include <IO/ValueSource.h>
#include <Cluster/MulticastPipe.h>
#include <Geometry/Vector.h>
#include <Geometry/AffineTransformation.h>
#include <Geometry/LinearUnit.h>
#include <GL/gl.h>
#include <GL/GLColorTemplates.h>
#include <GL/GLMaterial.h>
#include <GL/GLContextData.h>
#include <GL/GLValueCoders.h>
#include <GL/GLGeometryWrappers.h>
#include <GL/GLTransformationWrappers.h>
#include <GL/GLModels.h>
#include <GL/GLColorMap.h>
#include <GL/GLFrustum.h>
#include <GL/GLPrintError.h>
#include <GLMotif/StyleSheet.h>
#include <GLMotif/WidgetManager.h>
#include <GLMotif/Margin.h>
#include <GLMotif/Separator.h>
#include <GLMotif/Button.h>
#include <GLMotif/CascadeButton.h>
#include <GLMotif/Label.h>
#include <GLMotif/RowColumn.h>
#include <GLMotif/Menu.h>
#include <GLMotif/SubMenu.h>
#include <GLMotif/Popup.h>
#include <GLMotif/PopupMenu.h>
#include <GLMotif/Pager.h>
#include <GLMotif/PopupWindow.h>
#include <SceneGraph/NodeCreator.h>
#include <SceneGraph/VRMLFile.h>
#include <SceneGraph/TransformNode.h>
#include <Vrui/GlyphRenderer.h>
#include <Vrui/AffineCoordinateTransform.h>
#include <Vrui/Lightsource.h>
#include <Vrui/LightsourceManager.h>
#include <Vrui/Viewer.h>
#include <Vrui/CoordinateManager.h>
#include <Vrui/ToolManager.h>
#include <Vrui/ClusterSupport.h>
#include <Vrui/SurfaceNavigationTool.h>
#include <Vrui/OpenFile.h>
#include <Vrui/SceneGraphSupport.h>

#include "LidarOctree.h"
#include "LidarTool.h"
#include "LidarSelectionSaver.h"
#include "Primitive.h"
#include "PlanePrimitive.h"
#include "BruntonPrimitive.h"
#include "LinePrimitive.h"
#include "PointPrimitive.h"
#include "SpherePrimitive.h"
#include "CylinderPrimitive.h"
#include "PointClassifier.h"
#include "RidgeFinder.h"
#include "ProfileTool.h"
#include "SceneGraph.h"
#include "LoadPointSet.h"
#include "FallingSphereProcessor.h"

/**************
Helper classes:
**************/

namespace Misc {

template <>
class ValueCoder<LidarViewer::SelectorLocator::SelectorMode>
	{
	/* Methods: */
	public:
	static std::string encode(const LidarViewer::SelectorLocator::SelectorMode& value)
		{
		if(value==LidarViewer::SelectorLocator::ADD)
			return "Add";
		else
			return "Subtract";
		}
	static LidarViewer::SelectorLocator::SelectorMode decode(const char* start,const char* end,const char** decodeEnd =0)
		{
		if(strncasecmp(start,"Add",end-start)==0)
			{
			if(decodeEnd!=0)
				*decodeEnd=start+3;
			return LidarViewer::SelectorLocator::ADD;
			}
		else if(strncasecmp(start,"Subtract",end-start)==0)
			{
			if(decodeEnd!=0)
				*decodeEnd=start+8;
			return LidarViewer::SelectorLocator::SUBTRACT;
			}
		else
			throw DecodingError(Misc::printStdErrMsg("Could not convert %s to LidarViewer::SelectorLocator::SelectorMode",std::string(start,end).c_str()));
		}
	};

}

/*********************************************
Methods of class LidarViewer::SelectorLocator:
*********************************************/

LidarViewer::SelectorLocator::SelectorLocator(Vrui::LocatorTool* sTool,LidarViewer* sApplication,Misc::ConfigurationFileSection* cfg)
	:Locator(sTool,sApplication),
	 influenceRadius(application->brushSize),
	 selectorMode(application->defaultSelectorMode),
	 ready(false),
	 active(false)
	{
	if(cfg!=0)
		{
		/* Load tool settings from configuration file section: */
		influenceRadius=cfg->retrieveValue<Vrui::Scalar>("./influenceRadius",influenceRadius);
		selectorMode=cfg->retrieveValue<SelectorMode>("./selectorMode",selectorMode);
		}
	}

LidarViewer::SelectorLocator::~SelectorLocator(void)
	{
	}

void LidarViewer::SelectorLocator::storeState(Misc::ConfigurationFileSection& configFileSection) const
	{
	configFileSection.storeValue<Vrui::Scalar>("./influenceRadius",influenceRadius);
	configFileSection.storeValue<SelectorMode>("./selectorMode",selectorMode);
	}

void LidarViewer::SelectorLocator::getName(std::string& name) const
	{
	name="Select Points";
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
				for(int i=0;i<application->numOctrees;++i)
					if(application->showOctrees[i])
						application->octrees[i]->selectPoints(LidarOctree::Interactor(modelCenter,modelRadius));
				break;
			
			case SUBTRACT:
				for(int i=0;i<application->numOctrees;++i)
					if(application->showOctrees[i])
						application->octrees[i]->deselectPoints(LidarOctree::Interactor(modelCenter,modelRadius));
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

LidarViewer::DataItem::DataItem(GLContextData& contextData)
	:influenceSphereDisplayListId(glGenLists(1)),
	 pbls(contextData)
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
	
	mainMenuSelectorModes=selectorModes;
	
	return selectorModesMenuPopup;
	}

GLMotif::Popup* LidarViewer::createSelectionMenu(void)
	{
	GLMotif::Popup* selectionMenuPopup=new GLMotif::Popup("SelectionMenuPopup",Vrui::getWidgetManager());
	
	GLMotif::SubMenu* selectionMenu=new GLMotif::SubMenu("SelectionMenu",selectionMenuPopup,false);
	
	GLMotif::Button* classifySelectionButton=new GLMotif::Button("ClassifySelectionButton",selectionMenu,"Classify Selection");
	classifySelectionButton->getSelectCallbacks().add(this,&LidarViewer::classifySelectionCallback);
	
	GLMotif::Button* saveSelectionButton=new GLMotif::Button("SaveSelectionButton",selectionMenu,"Save Selection...");
	saveSelectionButton->getSelectCallbacks().add(this,&LidarViewer::saveSelectionCallback);
	
	new GLMotif::Separator("Separator1",selectionMenu,GLMotif::Separator::HORIZONTAL,0.0f,GLMotif::Separator::LOWERED);
	
	GLMotif::Button* clearSelectionButton=new GLMotif::Button("ClearSelectionButton",selectionMenu,"Clear Selection");
	clearSelectionButton->getSelectCallbacks().add(this,&LidarViewer::clearSelectionCallback);
	
	selectionMenu->manageChild();
	
	return selectionMenuPopup;
	}

GLMotif::Popup* LidarViewer::createExtractionMenu(void)
	{
	GLMotif::Popup* extractionMenuPopup=new GLMotif::Popup("ExtractionMenuPopup",Vrui::getWidgetManager());
	
	GLMotif::SubMenu* extractionMenu=new GLMotif::SubMenu("ExtractionMenu",extractionMenuPopup,false);
	
	GLMotif::Button* extractPlaneButton=new GLMotif::Button("ExtractPlaneButton",extractionMenu,"Extract Plane");
	extractPlaneButton->getSelectCallbacks().add(this,&LidarViewer::extractPlaneCallback);
	
	GLMotif::Button* extractBruntonButton=new GLMotif::Button("ExtractBruntonButton",extractionMenu,"Indicate Strike+Dip");
	extractBruntonButton->getSelectCallbacks().add(this,&LidarViewer::extractBruntonCallback);
	
	GLMotif::Button* extractLineButton=new GLMotif::Button("ExtractLineButton",extractionMenu,"Extract Line");
	extractLineButton->getSelectCallbacks().add(this,&LidarViewer::extractLineCallback);
	
	GLMotif::Button* extractSphereButton=new GLMotif::Button("ExtractSphereButton",extractionMenu,"Extract Sphere");
	extractSphereButton->getSelectCallbacks().add(this,&LidarViewer::extractSphereCallback);
	
	GLMotif::Button* extractCylinderButton=new GLMotif::Button("ExtractCylinderButton",extractionMenu,"Extract Cylinder");
	extractCylinderButton->getSelectCallbacks().add(this,&LidarViewer::extractCylinderCallback);
	
	GLMotif::Button* intersectPrimitivesButton=new GLMotif::Button("IntersectPrimitivesButton",extractionMenu,"Intersect Primitives");
	intersectPrimitivesButton->getSelectCallbacks().add(this,&LidarViewer::intersectPrimitivesCallback);
	
	GLMotif::Button* loadPrimitivesButton=new GLMotif::Button("LoadPrimitivesButton",extractionMenu,"Load Primitives...");
	loadPrimitivesButton->getSelectCallbacks().add(this,&LidarViewer::loadPrimitivesCallback);
	
	GLMotif::Button* savePrimitivesButton=new GLMotif::Button("SavePrimitivesButton",extractionMenu,"Save Primitives...");
	savePrimitivesButton->getSelectCallbacks().add(this,&LidarViewer::savePrimitivesCallback);
	
	new GLMotif::Separator("Separator1",extractionMenu,GLMotif::Separator::HORIZONTAL,0.0f,GLMotif::Separator::LOWERED);
	
	GLMotif::Button* deleteSelectedPrimitivesButton=new GLMotif::Button("DeleteSelectedPrimitivesButton",extractionMenu,"Delete Selected Primitives");
	deleteSelectedPrimitivesButton->getSelectCallbacks().add(this,&LidarViewer::deleteSelectedPrimitivesCallback);
	
	GLMotif::Button* clearPrimitivesButton=new GLMotif::Button("ClearPrimitivesButton",extractionMenu,"Clear Primitives");
	clearPrimitivesButton->getSelectCallbacks().add(this,&LidarViewer::clearPrimitivesCallback);
	
	extractionMenu->manageChild();
	
	return extractionMenuPopup;
	}

GLMotif::Popup* LidarViewer::createDialogMenu(void)
	{
	GLMotif::Popup* dialogMenuPopup=new GLMotif::Popup("DialogMenuPopup",Vrui::getWidgetManager());
	
	GLMotif::SubMenu* dialogMenu=new GLMotif::SubMenu("DialogMenu",dialogMenuPopup,false);
	
	if(numOctrees>1)
		{
		GLMotif::Button* showOctreeDialogButton=new GLMotif::Button("ShowOctreeDialogButton",dialogMenu,"Show Octree Dialog");
		showOctreeDialogButton->getSelectCallbacks().add(this,&LidarViewer::showOctreeDialogCallback);
		}
	
	GLMotif::Button* showRenderDialogButton=new GLMotif::Button("ShowRenderDialogButton",dialogMenu,"Show Render Dialog");
	showRenderDialogButton->getSelectCallbacks().add(this,&LidarViewer::showRenderDialogCallback);
	
	GLMotif::Button* showInteractionDialogButton=new GLMotif::Button("ShowInteractionDialogButton",dialogMenu,"Show Interaction Dialog");
	showInteractionDialogButton->getSelectCallbacks().add(this,&LidarViewer::showInteractionDialogCallback);
	
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
	
	GLMotif::Button* centerDisplayButton=new GLMotif::Button("CenterDisplayButton",mainMenu,"Center Display");
	centerDisplayButton->getSelectCallbacks().add(this,&LidarViewer::centerDisplayCallback);
	
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

GLMotif::PopupWindow* LidarViewer::createOctreeDialog(void)
	{
	GLMotif::PopupWindow* octreeDialog=new GLMotif::PopupWindow("OctreeDialog",Vrui::getWidgetManager(),"Octree Selection");
	octreeDialog->setCloseButton(true);
	octreeDialog->setResizableFlags(false,false);
	octreeDialog->popDownOnClose();
	
	GLMotif::RowColumn* octreeSelection=new GLMotif::RowColumn("OctreeSelection",octreeDialog,false);
	octreeSelection->setOrientation(GLMotif::RowColumn::HORIZONTAL);
	octreeSelection->setPacking(GLMotif::RowColumn::PACK_GRID);
	octreeSelection->setNumMinorWidgets(1);
	
	char toggleName[20]="OctreeToggle";
	char toggleLabel[10];
	for(int i=0;i<numOctrees;++i)
		{
		char* tnPtr=Misc::print(i,toggleName+12+4);
		while(tnPtr>toggleName+12)
			*(--tnPtr)='0';
		GLMotif::ToggleButton* octreeToggle=new GLMotif::ToggleButton(toggleName,octreeSelection,Misc::print(i,toggleLabel+9));
		octreeToggle->setToggle(showOctrees[i]);
		octreeToggle->getValueChangedCallbacks().add(this,&LidarViewer::octreeSelectionCallback,i);
		}
	
	octreeSelection->manageChild();
	
	return octreeDialog;
	}

GLMotif::PopupWindow* LidarViewer::createRenderDialog(void)
	{
	const GLMotif::StyleSheet& ss=*Vrui::getWidgetManager()->getStyleSheet();
	
	GLMotif::PopupWindow* renderDialog=new GLMotif::PopupWindow("RenderDialog",Vrui::getWidgetManager(),"Render Settings");
	renderDialog->setCloseButton(true);
	renderDialog->setResizableFlags(true,false);
	renderDialog->popDownOnClose();
	
	GLMotif::Pager* renderPager=new GLMotif::Pager("RenderPager",renderDialog,false);
	
	/* Create a page with LOD settings: */
	renderPager->setNextPageName("LOD");
	
	GLMotif::RowColumn* lodBox=new GLMotif::RowColumn("LODBox",renderPager,false);
	lodBox->setOrientation(GLMotif::RowColumn::VERTICAL);
	lodBox->setPacking(GLMotif::RowColumn::PACK_TIGHT);
	lodBox->setNumMinorWidgets(2);
	
	/* Create a slider/textfield combo to change the rendering quality: */
	new GLMotif::Label("RenderQualityLabel",lodBox,"Render Quality");
	
	GLMotif::TextFieldSlider* renderQualitySlider=new GLMotif::TextFieldSlider("RenderQualitySlider",lodBox,6,ss.fontHeight*10.0f);
	renderQualitySlider->getTextField()->setFloatFormat(GLMotif::TextField::FIXED);
	renderQualitySlider->getTextField()->setFieldWidth(5);
	renderQualitySlider->getTextField()->setPrecision(2);
	renderQualitySlider->setValueRange(-3.0,3.0,0.01);
	renderQualitySlider->getSlider()->addNotch(0.0f);
	renderQualitySlider->setValue(double(renderQuality));
	renderQualitySlider->getValueChangedCallbacks().add(this,&LidarViewer::renderQualitySliderCallback);
	
	/* Create a slider/textfield combo to change the focus+context weight: */
	new GLMotif::Label("FncWeightLabel",lodBox,"Focus + Context");
	
	GLMotif::TextFieldSlider* fncWeightSlider=new GLMotif::TextFieldSlider("FncWeightSlider",lodBox,6,ss.fontHeight*10.0f);
	fncWeightSlider->getTextField()->setFloatFormat(GLMotif::TextField::FIXED);
	fncWeightSlider->getTextField()->setFieldWidth(5);
	fncWeightSlider->getTextField()->setPrecision(2);
	fncWeightSlider->setValueRange(0.0,2.0,0.01);
	fncWeightSlider->setValue(double(fncWeight));
	fncWeightSlider->getValueChangedCallbacks().add(this,&LidarViewer::fncWeightSliderCallback);
	
	/* Create a slider/textfield combo to change the point size: */
	new GLMotif::Label("PointSizeLabel",lodBox,"Point Size");
	
	GLMotif::TextFieldSlider* pointSizeSlider=new GLMotif::TextFieldSlider("PointSizeSlider",lodBox,6,ss.fontHeight*10.0f);
	pointSizeSlider->getTextField()->setFloatFormat(GLMotif::TextField::FIXED);
	pointSizeSlider->getTextField()->setFieldWidth(4);
	pointSizeSlider->getTextField()->setPrecision(1);
	pointSizeSlider->setValueRange(1.0,10.0,0.5);
	pointSizeSlider->setValue(double(pointSize));
	pointSizeSlider->getValueChangedCallbacks().add(this,&LidarViewer::pointSizeSliderCallback);
	
	for(int i=0;i<3;++i)
		lodBox->setRowWeight(i,1.0f);
	
	lodBox->manageChild();
	
	/* Check if any of the octrees have normal vectors: */
	bool haveNormalVectors=false;
	for(int i=0;i<numOctrees;++i)
		haveNormalVectors=haveNormalVectors||octrees[i]->hasNormalVectors();
	if(haveNormalVectors)
		{
		/* Create a page with lighting settings: */
		renderPager->setNextPageName("Lighting");
		
		GLMotif::RowColumn* lightingBox=new GLMotif::RowColumn("LightingBox",renderPager,false);
		lightingBox->setOrientation(GLMotif::RowColumn::VERTICAL);
		lightingBox->setPacking(GLMotif::RowColumn::PACK_TIGHT);
		lightingBox->setNumMinorWidgets(2);
		
		/* Create a toggle button to enable point-based lighting: */
		GLMotif::Margin* enableLightingMargin=new GLMotif::Margin("EnableLightingMargin",lightingBox,false);
		enableLightingMargin->setAlignment(GLMotif::Alignment(GLMotif::Alignment::LEFT,GLMotif::Alignment::VCENTER));
		
		GLMotif::ToggleButton* enableLightingToggle=new GLMotif::ToggleButton("EnableLightingToggle",enableLightingMargin,"Lighting");
		enableLightingToggle->setBorderWidth(0.0f);
		enableLightingToggle->setHAlignment(GLFont::Left);
		enableLightingToggle->setToggle(pointBasedLighting);
		enableLightingToggle->getValueChangedCallbacks().add(this,&LidarViewer::enableLightingCallback);
		
		enableLightingMargin->manageChild();
		
		/* Create a toggle button to enable point colors: */
		GLMotif::Margin* usePointColorsMargin=new GLMotif::Margin("UsePointColorsMargin",lightingBox,false);
		usePointColorsMargin->setAlignment(GLMotif::Alignment(GLMotif::Alignment::LEFT,GLMotif::Alignment::VCENTER));
		
		GLMotif::ToggleButton* usePointColorsToggle=new GLMotif::ToggleButton("UsePointColorsToggle",usePointColorsMargin,"Use Point Colors");
		usePointColorsToggle->setBorderWidth(0.0f);
		usePointColorsToggle->setHAlignment(GLFont::Left);
		usePointColorsToggle->setToggle(usePointColors);
		usePointColorsToggle->getValueChangedCallbacks().add(this,&LidarViewer::usePointColorsCallback);
		
		usePointColorsMargin->manageChild();
		
		/* Create a button to enable point splatting: */
		GLMotif::Margin* useSplattingMargin=new GLMotif::Margin("UseSplattingMargin",lightingBox,false);
		useSplattingMargin->setAlignment(GLMotif::Alignment(GLMotif::Alignment::LEFT,GLMotif::Alignment::VCENTER));
		
		GLMotif::ToggleButton* useSplattingToggle=new GLMotif::ToggleButton("UseSplattingToggle",useSplattingMargin,"Splat Size");
		useSplattingToggle->setBorderWidth(0.0f);
		useSplattingToggle->setHAlignment(GLFont::Left);
		useSplattingToggle->setToggle(useSplatting);
		useSplattingToggle->getValueChangedCallbacks().add(this,&LidarViewer::useSplattingCallback);
		
		useSplattingMargin->manageChild();
		
		GLMotif::TextFieldSlider* splatSizeSlider=new GLMotif::TextFieldSlider("SplatSizeSlider",lightingBox,6,ss.fontHeight*10.0f);
		splatSizeSlider->getTextField()->setFloatFormat(GLMotif::TextField::FIXED);
		splatSizeSlider->getTextField()->setFieldWidth(6);
		splatSizeSlider->getTextField()->setPrecision(3);
		splatSizeSlider->setSliderMapping(GLMotif::TextFieldSlider::EXP10);
		splatSizeSlider->setValueRange(0.001,1000.0,0.01);
		splatSizeSlider->setValue(double(splatSize));
		splatSizeSlider->getValueChangedCallbacks().add(this,&LidarViewer::splatSizeSliderCallback);
		
		/* Create a toggle button to enable a fixed-position light source: */
		GLMotif::Margin* enableSunMargin=new GLMotif::Margin("EnableSunMargin",lightingBox,false);
		enableSunMargin->setAlignment(GLMotif::Alignment(GLMotif::Alignment::LEFT,GLMotif::Alignment::VCENTER));
		
		GLMotif::ToggleButton* enableSunToggle=new GLMotif::ToggleButton("SunToggle",enableSunMargin,"Sun Light Source");
		enableSunToggle->setBorderWidth(0.0f);
		enableSunToggle->setHAlignment(GLFont::Left);
		enableSunToggle->setToggle(enableSun);
		enableSunToggle->getValueChangedCallbacks().add(this,&LidarViewer::enableSunCallback);
		
		enableSunMargin->manageChild();
		
		/* Create a box with a pair of sliders to control the sun light source: */
		GLMotif::RowColumn* sunBox=new GLMotif::RowColumn("SunBox",lightingBox,false);
		sunBox->setOrientation(GLMotif::RowColumn::VERTICAL);
		sunBox->setNumMinorWidgets(2);
		
		new GLMotif::Label("SunAzimuthLabel",sunBox,"Azimuth");
		
		GLMotif::TextFieldSlider* sunAzimuthSlider=new GLMotif::TextFieldSlider("SunAzimuthSlider",sunBox,6,ss.fontHeight*10.0f);
		sunAzimuthSlider->getTextField()->setFloatFormat(GLMotif::TextField::FIXED);
		sunAzimuthSlider->getTextField()->setFieldWidth(3);
		sunAzimuthSlider->getTextField()->setPrecision(0);
		sunAzimuthSlider->setValueRange(0.0,360.0,1.0);
		sunAzimuthSlider->setValue(double(sunAzimuth));
		sunAzimuthSlider->getValueChangedCallbacks().add(this,&LidarViewer::sunAzimuthSliderCallback);
		
		new GLMotif::Label("SunElevationLabel",sunBox,"Elevation");
		
		GLMotif::TextFieldSlider* sunElevationSlider=new GLMotif::TextFieldSlider("SunElevationSlider",sunBox,6,ss.fontHeight*10.0f);
		sunElevationSlider->getTextField()->setFloatFormat(GLMotif::TextField::FIXED);
		sunElevationSlider->getTextField()->setFieldWidth(2);
		sunElevationSlider->getTextField()->setPrecision(0);
		sunElevationSlider->setValueRange(0.0,90.0,1.0);
		sunElevationSlider->setValue(double(sunElevation));
		sunElevationSlider->getValueChangedCallbacks().add(this,&LidarViewer::sunElevationSliderCallback);
		
		sunBox->manageChild();
		
		for(int i=0;i<3;++i)
			lightingBox->setRowWeight(i,1.0f);
		
		lightingBox->manageChild();
		
		/* Create a page to edit material properties: */
		renderPager->setNextPageName("Material");
		
		GLMotif::Margin* materialMargin=new GLMotif::Margin("MaterialMargin",renderPager,false);
		materialMargin->setAlignment(GLMotif::Alignment(GLMotif::Alignment::HCENTER,GLMotif::Alignment::VCENTER));
		
		GLMotif::MaterialEditor* materialEditor=new GLMotif::MaterialEditor("MaterialEditor",materialMargin);
		materialEditor->setMaterial(surfaceMaterial);
		materialEditor->getValueChangedCallbacks().add(this,&LidarViewer::materialCallback);
		
		materialMargin->manageChild();
		}
	
	/* Create a page with plane setting: */
	renderPager->setNextPageName("Plane");
	
	GLMotif::RowColumn* planeBox=new GLMotif::RowColumn("PlaneBox",renderPager,false);
	planeBox->setOrientation(GLMotif::RowColumn::VERTICAL);
	planeBox->setPacking(GLMotif::RowColumn::PACK_TIGHT);
	planeBox->setNumMinorWidgets(2);
	
	new GLMotif::Label("ExaggerationLabel",planeBox,"Exaggeration");
	
	GLMotif::TextFieldSlider* exaggerationSlider=new GLMotif::TextFieldSlider("ExaggerationSlider",planeBox,8,ss.fontHeight*10.0f);
	exaggerationSlider->getTextField()->setFieldWidth(8);
	exaggerationSlider->getTextField()->setPrecision(3);
	exaggerationSlider->setSliderMapping(GLMotif::TextFieldSlider::EXP10);
	exaggerationSlider->setValueRange(0.05,20.0,0.02);
	exaggerationSlider->getSlider()->addNotch(0.0f);
	exaggerationSlider->setValue(double(planeDistanceExaggeration));
	exaggerationSlider->getValueChangedCallbacks().add(this,&LidarViewer::distanceExaggerationSliderCallback);
	
	/* Create a toggle button to enable plane distance visualization: */
	GLMotif::Margin* enableTexturePlaneMargin=new GLMotif::Margin("EnableTexturePlaneMargin",planeBox,false);
	enableTexturePlaneMargin->setAlignment(GLMotif::Alignment(GLMotif::Alignment::LEFT,GLMotif::Alignment::VCENTER));
	
	GLMotif::ToggleButton* enableTexturePlaneToggle=new GLMotif::ToggleButton("EnableTexturePlaneToggle",enableTexturePlaneMargin,"Show Plane Distance");
	enableTexturePlaneToggle->setBorderWidth(0.0f);
	enableTexturePlaneToggle->setHAlignment(GLFont::Left);
	enableTexturePlaneToggle->setToggle(useTexturePlane);
	enableTexturePlaneToggle->getValueChangedCallbacks().add(this,&LidarViewer::enableTexturePlaneCallback);
	
	enableTexturePlaneMargin->manageChild();
	
	/* Create a slider to select the plane distance visualization scale: */
	GLMotif::TextFieldSlider* texturePlaneScaleSlider=new GLMotif::TextFieldSlider("TexturePlaneScaleSlider",planeBox,8,ss.fontHeight*10.0f);
	texturePlaneScaleSlider->getTextField()->setFieldWidth(8);
	texturePlaneScaleSlider->getTextField()->setPrecision(3);
	texturePlaneScaleSlider->setSliderMapping(GLMotif::TextFieldSlider::EXP10);
	texturePlaneScaleSlider->setValueRange(0.01,10000.0,0.1);
	texturePlaneScaleSlider->setValue(double(texturePlaneScale));
	texturePlaneScaleSlider->getValueChangedCallbacks().add(this,&LidarViewer::texturePlaneScaleSliderCallback);
	
	#ifdef LIDARVIEWER_VISUALIZE_WATER
	
	/* Create a slider to adjust the vertical plane position: */
	new GLMotif::Label("PlaneOffsetLabel",planeBox,"Plane Offset");
	
	GLMotif::TextFieldSlider* texturePlaneOffsetSlider=new GLMotif::TextFieldSlider("TexturePlaneOffsetSlider",planeBox,8,ss.fontHeight*10.0f);
	texturePlaneOffsetSlider->getTextField()->setFieldWidth(8);
	texturePlaneOffsetSlider->getTextField()->setPrecision(3);
	texturePlaneOffsetSlider->setSliderMapping(GLMotif::TextFieldSlider::LINEAR);
	texturePlaneOffsetSlider->setValueRange(-20.0,20.0,0.01);
	texturePlaneOffsetSlider->getSlider()->addNotch(0.0f);
	texturePlaneOffsetSlider->setValue(double(texturePlaneOffset));
	texturePlaneOffsetSlider->getValueChangedCallbacks().add(this,&LidarViewer::texturePlaneOffsetSliderCallback);
	
	for(int i=0;i<3;++i)
		planeBox->setRowWeight(i,1.0f);
	
	#else
	
	for(int i=0;i<2;++i)
		planeBox->setRowWeight(i,1.0f);
	
	#endif
	
	planeBox->manageChild();
	
	renderPager->manageChild();
	
	return renderDialog;
	}

GLMotif::PopupWindow* LidarViewer::createInteractionDialog(void)
	{
	const GLMotif::StyleSheet& ss=*Vrui::getWidgetManager()->getStyleSheet();
	
	GLMotif::PopupWindow* interactionDialog=new GLMotif::PopupWindow("InteractionDialog",Vrui::getWidgetManager(),"Interaction Settings");
	interactionDialog->setCloseButton(true);
	interactionDialog->setResizableFlags(true,false);
	interactionDialog->popDownOnClose();
	
	GLMotif::RowColumn* interactionSettings=new GLMotif::RowColumn("InteractionSettings",interactionDialog,false);
	interactionSettings->setOrientation(GLMotif::RowColumn::VERTICAL);
	interactionSettings->setPacking(GLMotif::RowColumn::PACK_TIGHT);
	interactionSettings->setNumMinorWidgets(1);
	
	/* Create a box with tool settings: */
	GLMotif::Margin* toolSettingsMargin=new GLMotif::Margin("ToolSettingsMargin",interactionSettings,false);
	toolSettingsMargin->setAlignment(GLMotif::Alignment(GLMotif::Alignment::LEFT,GLMotif::Alignment::VCENTER));
	
	GLMotif::RowColumn* toolSettingsBox=new GLMotif::RowColumn("ToolSettingsBox",toolSettingsMargin,false);
	toolSettingsBox->setOrientation(GLMotif::RowColumn::HORIZONTAL);
	toolSettingsBox->setPacking(GLMotif::RowColumn::PACK_TIGHT);
	toolSettingsBox->setNumMinorWidgets(1);
	
	/* Create a toggle button to override existing tools' settings: */
	GLMotif::ToggleButton* overrideToolsToggle=new GLMotif::ToggleButton("OverrideToolsToggle",toolSettingsBox,"Override Tools");
	overrideToolsToggle->setBorderWidth(0.0f);
	overrideToolsToggle->setHAlignment(GLFont::Left);
	overrideToolsToggle->setToggle(overrideTools);
	overrideToolsToggle->getValueChangedCallbacks().add(this,&LidarViewer::overrideToolsCallback);
	
	new GLMotif::Separator("Separator1",toolSettingsBox,GLMotif::Separator::VERTICAL,0.0f,GLMotif::Separator::LOWERED);
	
	/* Create a radio box to select selection modes: */
	interactionDialogSelectorModes=new GLMotif::RadioBox("InteractionDialogSelectorModes",toolSettingsBox,false);
	interactionDialogSelectorModes->setOrientation(GLMotif::RowColumn::HORIZONTAL);
	interactionDialogSelectorModes->setPacking(GLMotif::RowColumn::PACK_TIGHT);
	interactionDialogSelectorModes->setSelectionMode(GLMotif::RadioBox::ALWAYS_ONE);
	
	interactionDialogSelectorModes->addToggle("Add");
	interactionDialogSelectorModes->addToggle("Subtract");
	
	switch(defaultSelectorMode)
		{
		case SelectorLocator::ADD:
			interactionDialogSelectorModes->setSelectedToggle(0);
			break;
		
		case SelectorLocator::SUBTRACT:
			interactionDialogSelectorModes->setSelectedToggle(1);
			break;
		}
	interactionDialogSelectorModes->getValueChangedCallbacks().add(this,&LidarViewer::changeSelectorModeCallback);
	interactionDialogSelectorModes->manageChild();
	
	toolSettingsBox->manageChild();
	toolSettingsMargin->manageChild();
	
	/* Create a box with sliders to adjust interaction sizes: */
	GLMotif::RowColumn* sliderBox=new GLMotif::RowColumn("SliderBox",interactionSettings,false);
	sliderBox->setOrientation(GLMotif::RowColumn::VERTICAL);
	sliderBox->setPacking(GLMotif::RowColumn::PACK_TIGHT);
	sliderBox->setNumMinorWidgets(2);
	
	/* Create a slider to change the size of the selection brush: */
	new GLMotif::Label("BrushSizeLabel",sliderBox,"Brush Size");
	
	GLMotif::TextFieldSlider* brushSizeSlider=new GLMotif::TextFieldSlider("BrushSizeSlider",sliderBox,8,ss.fontHeight*10.0f);
	brushSizeSlider->getTextField()->setFieldWidth(7);
	brushSizeSlider->getTextField()->setPrecision(4);
	brushSizeSlider->setValueRange(double(brushSize)*0.1,double(brushSize)*5.0,double(brushSize)*0.01);
	brushSizeSlider->setValue(double(brushSize));
	brushSizeSlider->getValueChangedCallbacks().add(this,&LidarViewer::brushSizeSliderCallback);
	
	/* Create a slider to change the size of the processing neighborhood: */
	new GLMotif::Label("NeighborhoodSizeLabel",sliderBox,"Neighborhood Size");
	
	GLMotif::TextFieldSlider* neighborhoodSizeSlider=new GLMotif::TextFieldSlider("NeighborhoodSizeSlider",sliderBox,8,ss.fontHeight*10.0f);
	neighborhoodSizeSlider->getTextField()->setFieldWidth(7);
	neighborhoodSizeSlider->getTextField()->setPrecision(4);
	neighborhoodSizeSlider->setSliderMapping(GLMotif::TextFieldSlider::EXP10);
	neighborhoodSizeSlider->setValueRange(10.0e-3,10.0e3,0.1);
	neighborhoodSizeSlider->setValue(double(neighborhoodSize));
	neighborhoodSizeSlider->getValueChangedCallbacks().add(this,&LidarViewer::neighborhoodSizeSliderCallback);
	
	sliderBox->manageChild();
	
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
	/* Set the new primitive's color: */
	newPrimitive->setSurfaceColor(primitiveColor);
	newPrimitive->setGridColor(Primitive::Color(0.2f,0.2f,0.2f));
	
	/* Store the new primitive: */
	primitives.push_back(newPrimitive);
	primitiveSelectedFlags.push_back(false);
	return int(primitives.size())-1;
	}

void LidarViewer::selectPrimitive(int primitiveIndex)
	{
	if(!primitiveSelectedFlags[primitiveIndex])
		{
		primitives[primitiveIndex]->setSurfaceColor(selectedPrimitiveColor);
		primitiveSelectedFlags[primitiveIndex]=true;
		}
	}

void LidarViewer::deselectPrimitive(int primitiveIndex)
	{
	if(primitiveSelectedFlags[primitiveIndex])
		{
		primitives[primitiveIndex]->setSurfaceColor(primitiveColor);
		primitiveSelectedFlags[primitiveIndex]=false;
		}
	}

void LidarViewer::deletePrimitive(int primitiveIndex)
	{
	delete primitives[primitiveIndex];
	primitives.erase(primitives.begin()+primitiveIndex);
	primitiveSelectedFlags.erase(primitiveSelectedFlags.begin()+primitiveIndex);
	}

void LidarViewer::updateSun(void)
	{
	/* Compute the light source's direction vector: */
	Vrui::Scalar z=Math::sin(Math::rad(sunElevation));
	Vrui::Scalar xy=Math::cos(Math::rad(sunElevation));
	Vrui::Scalar x=xy*Math::sin(Math::rad(sunAzimuth));
	Vrui::Scalar y=xy*Math::cos(Math::rad(sunAzimuth));
	sun->getLight().position=GLLight::Position(GLLight::Scalar(x),GLLight::Scalar(y),GLLight::Scalar(z),GLLight::Scalar(0));
	}

void LidarViewer::setEnableSun(bool newEnableSun)
	{
	if(enableSun!=newEnableSun)
		{
		enableSun=newEnableSun;
		
		if(enableSun)
			{
			/* Store the headlight enable states of all viewers and disable all headlights: */
			viewerHeadlightStates=new bool[Vrui::getNumViewers()];
			for(int i=0;i<Vrui::getNumViewers();++i)
				{
				viewerHeadlightStates[i]=Vrui::getViewer(i)->getHeadlight().isEnabled();
				Vrui::getViewer(i)->setHeadlightState(false);
				}
			
			/* Enable and update the sun light source: */
			sun->enable();
			updateSun();
			}
		else
			{
			/* Reset the headlight enable states of all viewers: */
			for(int i=0;i<Vrui::getNumViewers();++i)
				Vrui::getViewer(i)->setHeadlightState(viewerHeadlightStates[i]);
			delete[] viewerHeadlightStates;
			viewerHeadlightStates=0;
			
			/* Disable the sun light source: */
			sun->disable();
			}
		}
	}

LidarViewer::LidarViewer(int& argc,char**& argv,char**& appDefaults)
	:Vrui::Application(argc,argv,appDefaults),
	 numOctrees(0),octrees(0),showOctrees(0),
	 coordTransform(0),
	 renderQuality(0),fncWeight(0.5),
	 pointSize(3.0f),
	 pointBasedLighting(false),
	 surfaceMaterial(GLMaterial::Color(1.0f,1.0f,1.0f),GLMaterial::Color(0.6f,0.6f,0.6f),30.0f),
	 usePointColors(true),useSplatting(false),splatSize(0.1),
	 enableSun(false),viewerHeadlightStates(0),sunAzimuth(180),sunElevation(45),sun(0),
	 useTexturePlane(false),texturePlane(GPlane::Vector(0.0,0.0,1.0),0.0),texturePlaneScale(100.0),
	 #ifdef LIDARVIEWER_VISUALIZE_WATER
	 texturePlaneOffset(0.0),
	 #endif
	 planeDistanceExaggeration(1.0),
	 updateTree(true),
	 lastFrameTime(Vrui::getApplicationTime()),
	 overrideTools(true),
	 brushSize(Vrui::getGlyphRenderer()->getGlyphSize()*Vrui::Scalar(2.5)),
	 brushColor(0.6f,0.6f,0.1f,0.5f),
	 defaultSelectorMode(SelectorLocator::ADD),
	 neighborhoodSize(1),
	 extractorPipe(Vrui::openPipe()),
	 primitiveColor(0.5f,0.5f,0.1f,0.5f),
	 selectedPrimitiveColor(0.1f,0.5f,0.5f,0.5f),
	 lastPickedPrimitive(-1),
	 mainMenu(0),octreeDialog(0),renderDialog(0),interactionDialog(0),
	 sceneGraphRoot(0),
	 dataDirectory(0)
	{
	memCacheSize=512;
	unsigned int gfxCacheSize=128;
	try
		{
		/* Open LidarViewer's configuration file: */
		Misc::ConfigurationFile configFile(LIDARVIEWER_CONFIGFILENAME);
		Misc::ConfigurationFileSection cfg=configFile.getSection("/LidarViewer");
		
		/* Override program settings from configuration file: */
		renderQuality=cfg.retrieveValue<Scalar>("./renderQuality",renderQuality);
		fncWeight=cfg.retrieveValue<Scalar>("./focusAndContextWeight",fncWeight);
		pointSize=cfg.retrieveValue<float>("./pointSize",pointSize);
		pointBasedLighting=cfg.retrieveValue<bool>("./enableLighting",pointBasedLighting);
		surfaceMaterial=cfg.retrieveValue<GLMaterial>("./surfaceMaterial",surfaceMaterial);
		usePointColors=cfg.retrieveValue<bool>("./usePointColors",usePointColors);
		useSplatting=cfg.retrieveValue<bool>("./useSplatting",useSplatting);
		splatSize=cfg.retrieveValue<double>("./splatSize",splatSize);
		enableSun=cfg.retrieveValue<bool>("./enableSun",enableSun);
		sunAzimuth=cfg.retrieveValue<Scalar>("./sunAzimuth",sunAzimuth);
		sunElevation=cfg.retrieveValue<Scalar>("./sunElevation",sunElevation);
		overrideTools=cfg.retrieveValue<bool>("./overrideTools",overrideTools);
		brushSize=cfg.retrieveValue<Vrui::Scalar>("./brushSize",brushSize);
		brushColor=cfg.retrieveValue<GLColor<GLfloat,4> >("./brushColor",brushColor);
		primitiveColor=cfg.retrieveValue<GLColor<GLfloat,4> >("./primitiveColor",primitiveColor);
		selectedPrimitiveColor=cfg.retrieveValue<GLColor<GLfloat,4> >("./selectedPrimitiveColor",selectedPrimitiveColor);
		memCacheSize=cfg.retrieveValue<unsigned int>("./memoryCacheSize",memCacheSize);
		gfxCacheSize=cfg.retrieveValue<unsigned int>("./graphicsCacheSize",gfxCacheSize);
		}
	catch(std::runtime_error err)
		{
		/* Just ignore the error */
		}
	
	/* Parse the command line: */
	for(int i=1;i<argc;++i)
		{
		if(argv[i][0]=='-')
			{
			if(strcasecmp(argv[i]+1,"memoryCacheSize")==0)
				{
				if(i+1<argc)
					{
					++i;
					memCacheSize=atoi(argv[i]);
					}
				}
			else if(strcasecmp(argv[i]+1,"graphicsCacheSize")==0)
				{
				if(i+1<argc)
					{
					++i;
					gfxCacheSize=atoi(argv[i]);
					}
				}
			else if(strcasecmp(argv[i]+1,"renderQuality")==0)
				{
				if(i+1<argc)
					{
					++i;
					renderQuality=Scalar(atof(argv[i]));
					}
				}
			else if(strcasecmp(argv[i]+1,"focusAndContextWeight")==0)
				{
				if(i+1<argc)
					{
					++i;
					fncWeight=Scalar(atof(argv[i]));
					}
				}
			else if(strcasecmp(argv[i]+1,"pointSize")==0)
				{
				if(i+1<argc)
					{
					++i;
					pointSize=float(atof(argv[i]));
					}
				}
			else if(strcasecmp(argv[i]+1,"enableLighting")==0)
				pointBasedLighting=true;
			else if(strcasecmp(argv[i]+1,"usePointColors")==0)
				usePointColors=true;
			else if(strcasecmp(argv[i]+1,"sceneGraph")==0)
				{
				if(i+1<argc)
					{
					++i;
					
					/* Create a root node if there is none yet: */
					if(sceneGraphRoot==0)
						sceneGraphRoot=new SceneGraph::TransformNode;
					
					/* Load the scene graph file: */
					SceneGraph::NodeCreator nodeCreator;
					SceneGraph::VRMLFile vrmlFile(argv[i],Vrui::openFile(argv[i]),nodeCreator,Vrui::getClusterMultiplexer());
					vrmlFile.parse(sceneGraphRoot);
					}
				}
			}
		else
			{
			/* Load another LiDAR octree file: */
			lidarFileNames.push_back(argv[i]);
			LidarOctree** newOctrees=new LidarOctree*[numOctrees+1];
			for(int j=0;j<numOctrees;++j)
				newOctrees[j]=octrees[j];
			delete[] octrees;
			++numOctrees;
			octrees=newOctrees;
			
			/* Load a LiDAR octree file: */
			octrees[numOctrees-1]=new LidarOctree(argv[i],size_t(memCacheSize)*size_t(1024*1024),size_t(gfxCacheSize)*size_t(1024*1024));
			}
		}
	
	if(numOctrees==0)
		Misc::throwStdErr("LidarViewer::LidarViewer: No octree file name provided");
	
	/* Initialize all octrees: */
	showOctrees=new bool[numOctrees];
	for(int i=0;i<numOctrees;++i)
		{
		octrees[i]->setRenderQuality(renderQuality);
		octrees[i]->setTreeUpdateFunction(treeUpdateNotificationCB,0);
		showOctrees[i]=true;
		}
	
	/* Check if all the octrees have the same linear unit: */
	Geometry::LinearUnit linearUnit;
	for(int i=0;i<numOctrees;++i)
		{
		/* Check if the LiDAR file contains a unit file: */
		std::string unitFileName=lidarFileNames[i];
		unitFileName.append("/Unit");
		if(Misc::isFileReadable(unitFileName.c_str()))
			{
			/* Read the unit file: */
			IO::ValueSource unit(Vrui::openFile(unitFileName.c_str()));
			unit.skipWs();
			Vrui::Scalar unitFactor=Vrui::Scalar(unit.readNumber());
			std::string unitName=unit.readString();
			
			/* Create a linear unit: */
			Geometry::LinearUnit fileLinearUnit(unitName.c_str(),unitFactor);
			if(linearUnit.unit==Geometry::LinearUnit::UNKNOWN)
				linearUnit=fileLinearUnit;
			else if(linearUnit.unit!=fileLinearUnit.unit||linearUnit.factor!=fileLinearUnit.factor)
				Misc::throwStdErr("LidarViewer::LidarViewer: Octree file %s has mismatching units",lidarFileNames[i].c_str());
			}
		}
	
	/* Set the coordinate manager's linear unit: */
	Vrui::getCoordinateManager()->setUnit(linearUnit);
	
	/* Register a coordinate transform object to undo the coordinate offset done by the octree object and LiDAR preprocessor: */
	/* WARNING: This does not work properly for multiple octree files! */
	for(int i=0;i<3;++i)
		offsets[i]=double(octrees[0]->getPointOffset()[i]);
	std::string offsetFileName=lidarFileNames[0];
	offsetFileName.append("/Offset");
	if(Misc::isFileReadable(offsetFileName.c_str()))
		{
		/* Read the offset file: */
		IO::FilePtr offsetFile(Vrui::openFile(offsetFileName.c_str()));
		offsetFile->setEndianness(Misc::LittleEndian);
		for(int i=0;i<3;++i)
			offsets[i]-=offsetFile->read<double>();
		}
	Vrui::Vector offVec(offsets);
	coordTransform=new Vrui::AffineCoordinateTransform(Vrui::ATransform::translate(-offVec));
	Vrui::getCoordinateManager()->setCoordinateTransform(coordTransform);
	
	/* Apply the transformation to any additional scene graphs: */
	if(sceneGraphRoot!=0)
		{
		sceneGraphRoot->translation.setValue(-offVec);
		sceneGraphRoot->update();
		}
	
	/* Create the sun lightsource: */
	sun=Vrui::getLightsourceManager()->createLightsource(false);
	bool newEnableSun=enableSun;
	enableSun=false;
	sun->disable();
	setEnableSun(newEnableSun);
	
	/* Create the GUI: */
	mainMenu=createMainMenu();
	Vrui::setMainMenu(mainMenu);
	if(numOctrees>1)
		octreeDialog=createOctreeDialog();
	renderDialog=createRenderDialog();
	interactionDialog=createInteractionDialog();
	
	/* Initialize the navigation transformation: */
	centerDisplayCallback(0);
	
	/* Register the custom tool classes with the Vrui tool manager: */
	LidarToolFactory* lidarToolFactory=new LidarToolFactory(*Vrui::getToolManager());
	Vrui::getToolManager()->addClass(lidarToolFactory,LidarToolFactory::factoryDestructor);
	ProfileToolFactory* profileToolFactory=new ProfileToolFactory(*Vrui::getToolManager());
	Vrui::getToolManager()->addClass(profileToolFactory,ProfileToolFactory::factoryDestructor);
	
	/* Initialize the scene graph: */
	createSceneGraph();
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
	
	/* Destroy the scene graph: */
	destroySceneGraph();
	
	/* Delete the GUI: */
	delete mainMenu;
	delete octreeDialog;
	delete renderDialog;
	delete interactionDialog;
	
	/* Delete the viewer headlight states: */
	delete[] viewerHeadlightStates;
	
	/* Delete all octrees: */
	for(int i=0;i<numOctrees;++i)
		delete octrees[i];
	delete[] octrees;
	delete[] showOctrees;
	}

void LidarViewer::initContext(GLContextData& contextData) const
	{
	/* Create a new context entry: */
	DataItem* dataItem=new DataItem(contextData);
	contextData.addDataItem(this,dataItem);
	
	/* Create the influence sphere display list: */
	glNewList(dataItem->influenceSphereDisplayListId,GL_COMPILE);
	glDisable(GL_CULL_FACE);
	glEnable(GL_BLEND);
	glDepthMask(GL_FALSE);
	glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
	glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
	glColor(brushColor);
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
	const int numPlaneColors=7;
	static const GLColor<GLfloat,3> planeColors[numPlaneColors]=
		{
		GLColor<GLfloat,3>(0.0f,0.0f,0.5f),
		GLColor<GLfloat,3>(0.0f,0.0f,1.0f),
		GLColor<GLfloat,3>(0.0f,1.0f,1.0f),
		GLColor<GLfloat,3>(1.0f,1.0f,1.0f),
		GLColor<GLfloat,3>(1.0f,1.0f,0.0f),
		GLColor<GLfloat,3>(1.0f,0.0f,0.0f),
		GLColor<GLfloat,3>(0.5f,0.0f,0.0f)
		};
	GLColor<GLfloat,3>* planeColorMap=new GLColor<GLfloat,3>[1024];
	for(int i=0;i<1024;++i)
		{
		int ci0=(i*(numPlaneColors-1))/1023;
		int ci1=ci0<numPlaneColors-1?ci0+1:numPlaneColors-1;
		float cd=float(i*(numPlaneColors-1)-ci0*1023)/1023.0f;
		for(int j=0;j<3;++j)
			planeColorMap[i][j]=planeColors[ci0][j]*(1.0f-cd)+planeColors[ci1][j]*cd;
		}
	
	#if 0
	
	/* Add notches to the color map: */
	for(int i=0;i<=20;++i)
		{
		planeColorMap[(i*1023+10)/20]=GLColor<GLfloat,3>(0.0f,0.0f,0.0f);
		}
	
	#endif
	
	/* Create the color map texture image: */
	glBindTexture(GL_TEXTURE_1D,dataItem->planeColorMapTextureId);
	glTexParameteri(GL_TEXTURE_1D,GL_TEXTURE_MIN_FILTER,GL_NEAREST);
	glTexParameteri(GL_TEXTURE_1D,GL_TEXTURE_MAG_FILTER,GL_NEAREST);
	glTexParameteri(GL_TEXTURE_1D,GL_TEXTURE_BASE_LEVEL,0);
	glTexParameteri(GL_TEXTURE_1D,GL_TEXTURE_MAX_LEVEL,0);
	glTexParameteri(GL_TEXTURE_1D,GL_TEXTURE_WRAP_S,GL_CLAMP_TO_EDGE);
	glTexImage1D(GL_TEXTURE_1D,0,GL_RGB,1024,0,GL_RGB,GL_FLOAT,planeColorMap);
	glBindTexture(GL_TEXTURE_1D,0);
	
	delete[] planeColorMap;
	}

void LidarViewer::toolCreationCallback(Vrui::ToolManager::ToolCreationCallbackData* cbData)
	{
	/* Let the base class at it first: */
	Vrui::Application::toolCreationCallback(cbData);
	
	/* Check if the new tool is a locator tool: */
	Vrui::LocatorTool* ltool=dynamic_cast<Vrui::LocatorTool*>(cbData->tool);
	if(ltool!=0)
		{
		/* Create new locator: */
		Locator* newLocator=new SelectorLocator(ltool,this,cbData->cfg);
		
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
	
	/* Check if the new tool is a surface navigation tool: */
	Vrui::SurfaceNavigationTool* surfaceNavigationTool=dynamic_cast<Vrui::SurfaceNavigationTool*>(cbData->tool);
	if(surfaceNavigationTool!=0)
		{
		/* Set the new tool's alignment function: */
		surfaceNavigationTool->setAlignFunction(Misc::createFunctionCall(this,&LidarViewer::alignSurfaceFrame));
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
	Point displayCenter=Point(Vrui::getInverseNavigationTransformation().transform(Vrui::getDisplayCenter()));
	Scalar displaySize=Scalar(Vrui::getInverseNavigationTransformation().getScaling()*Vrui::getDisplaySize());
	for(int i=0;i<numOctrees;++i)
		if(showOctrees[i])
			{
			octrees[i]->startRenderPass();
			octrees[i]->setFocusAndContext(displayCenter,displaySize*Scalar(0.5),fncWeight);
			octrees[i]->setBaseSurfelSize(splatSize,float(Vrui::getNavigationTransformation().getScaling()));
			}
	
	/* Prepare for interaction with all interactors: */
	for(LocatorList::iterator lIt=locators.begin();lIt!=locators.end();++lIt)
		{
		SelectorLocator* sl=dynamic_cast<SelectorLocator*>(*lIt);
		if(sl!=0&&sl->ready)
			{
			for(int i=0;i<numOctrees;++i)
				if(showOctrees[i])
					octrees[i]->interact(LidarOctree::Interactor(sl->modelCenter,sl->modelRadius));
			}
		}
	}

void LidarViewer::display(GLContextData& contextData) const
	{
	/* Retrieve context entry: */
	DataItem* dataItem=contextData.retrieveDataItem<DataItem>(this);
	
	/* Set up basic OpenGL state: */
	glPushAttrib(GL_ENABLE_BIT|GL_LIGHTING_BIT|GL_LINE_BIT|GL_POINT_BIT|GL_TEXTURE_BIT);
	if(pointBasedLighting&&octrees[0]->hasNormalVectors())
		{
		if(useTexturePlane)
			glMaterial(GLMaterialEnums::FRONT_AND_BACK,GLMaterial(GLMaterial::Color(1.0f,1.0f,1.0f),surfaceMaterial.specular,surfaceMaterial.shininess));
		else if(usePointColors)
			{
			glMaterial(GLMaterialEnums::FRONT_AND_BACK,surfaceMaterial);
			glEnable(GL_COLOR_MATERIAL);
			glColorMaterial(GL_FRONT_AND_BACK,GL_AMBIENT_AND_DIFFUSE);
			}
		else
			glMaterial(GLMaterialEnums::FRONT_AND_BACK,surfaceMaterial);
		
		/* Enable the point-based lighting shader: */
		dataItem->pbls.setUsePlaneDistance(useTexturePlane);
		dataItem->pbls.setUsePointColors(usePointColors);
		dataItem->pbls.setUseSplatting(useSplatting);
		dataItem->pbls.enable();
		
		#if 0
		
		/* Set the shader's eye coordinate surfel size: */
		dataItem->pbls.setSurfelSize(float(Vrui::getNavigationTransformation().getScaling()*splatSize));
		
		#endif
		
		if(useTexturePlane)
			{
			/* Set up distance plane texturing: */
			#ifdef LIDARVIEWER_VISUALIZE_WATER
			GPlane offsetTexturePlane=texturePlane;
			offsetTexturePlane.setOffset(offsetTexturePlane.getOffset()+texturePlaneOffset);
			dataItem->pbls.setDistancePlane(0,offsetTexturePlane,texturePlaneScale);
			#else
			dataItem->pbls.setDistancePlane(0,texturePlane,texturePlaneScale);
			#endif
			
			/* Bind the distance map texture: */
			glBindTexture(GL_TEXTURE_1D,dataItem->planeColorMapTextureId);
			}
		}
	else
		{
		/* Turn off lighting: */
		glDisable(GL_LIGHTING);
		
		if(useTexturePlane)
			{
			/* Set up automatic texture coordinate generation: */
			glTexGeni(GL_S,GL_TEXTURE_GEN_MODE,GL_OBJECT_LINEAR);
			GLdouble planeCoeff[4];
			for(int i=0;i<3;++i)
				planeCoeff[i]=texturePlane.getNormal()[i]/texturePlaneScale;
			#ifdef LIDARVIEWER_VISUALIZE_WATER
			planeCoeff[3]=0.5-(texturePlane.getOffset()+texturePlaneOffset)/texturePlaneScale;
			#else
			planeCoeff[3]=0.5-texturePlane.getOffset()/texturePlaneScale;
			#endif
			glTexGendv(GL_S,GL_OBJECT_PLANE,planeCoeff);
			glEnable(GL_TEXTURE_GEN_S);
			
			/* Enable 1D texture mapping: */
			glEnable(GL_TEXTURE_1D);
			glDisable(GL_TEXTURE_2D);
			glDisable(GL_TEXTURE_3D);
			glBindTexture(GL_TEXTURE_1D,dataItem->planeColorMapTextureId);
			glTexEnvi(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,GL_REPLACE);
			}
		}
	
	glPointSize(pointSize);
	
	if(planeDistanceExaggeration!=1.0)
		{
		/* Enable plane distance exaggeration: */
		glMatrixMode(GL_MODELVIEW);
		glPushMatrix();
		GPlane::Vector fTrans=texturePlane.getNormal()*(texturePlane.getOffset()/Geometry::sqr(texturePlane.getNormal()));
		Geometry::Rotation<GPlane::Scalar,3> fRot=Geometry::Rotation<GPlane::Scalar,3>::rotateFromTo(GPlane::Vector(0,0,1),texturePlane.getNormal());
		glTranslate(fTrans);
		glRotate(fRot);
		glScaled(1.0,1.0,planeDistanceExaggeration);
		glRotate(Geometry::invert(fRot));
		glTranslate(-fTrans);
		}
	
	/* Render the LiDAR point tree: */
	LidarOctree::Frustum frustum;
	frustum.setFromGL();
	for(int i=0;i<numOctrees;++i)
		if(showOctrees[i])
			octrees[i]->glRenderAction(frustum,dataItem->pbls,contextData);
	
	if(planeDistanceExaggeration!=1.0)
		glPopMatrix();
	
	/* Reset OpenGL state: */
	if(pointBasedLighting&&octrees[0]->hasNormalVectors())
		{
		/* Disable the point-based lighting shader: */
		dataItem->pbls.disable();
		
		if(useTexturePlane)
			glBindTexture(GL_TEXTURE_1D,0);
		}
	else
		{
		if(useTexturePlane)
			{
			/* Disable 1D texture mapping: */
			glBindTexture(GL_TEXTURE_1D,0);
			glDisable(GL_TEXTURE_1D);
			
			/* Disable automatic texture coordinate generation: */
			glDisable(GL_TEXTURE_GEN_S);
			}
		}
	glPopAttrib();
	
	glPushAttrib(GL_ENABLE_BIT|GL_LIGHTING_BIT|GL_TEXTURE_BIT);
	
	/* Render LiDAR Viewer's own scene graph: */
	renderSceneGraph(contextData);
	
	/* Render any additional scene graphs: */
	if(sceneGraphRoot!=0)
		Vrui::renderSceneGraph(sceneGraphRoot.getPointer(),true,contextData);
	
	glPopAttrib();
	
	/* Render all locators: */
	for(LocatorList::const_iterator lIt=locators.begin();lIt!=locators.end();++lIt)
		(*lIt)->glRenderAction(contextData);
	
	/* Render all extracted primitives: */
	for(PrimitiveList::const_iterator pIt=primitives.begin();pIt!=primitives.end();++pIt)
		(*pIt)->glRenderAction(contextData);
	}

void LidarViewer::alignSurfaceFrame(Vrui::SurfaceNavigationTool::AlignmentData& alignmentData)
	{
	/* Get the frame's base point: */
	Vrui::Point base=alignmentData.surfaceFrame.getOrigin();
	
	Scalar surfaceZ=base[2];
	if(Vrui::isMaster())
		{
		/* Drop a sphere onto the LiDAR point cloud: */
		base[2]+=alignmentData.probeSize+alignmentData.maxClimb;
		FallingSphereProcessor fsp(base,alignmentData.probeSize);
		for(int i=0;i<numOctrees;++i)
			octrees[i]->processPointsInBox(fsp.getBox(),fsp);
		
		/* Get the surface elevation: */
		if(fsp.getMinZ()!=Math::Constants<Scalar>::min)
			surfaceZ=fsp.getMinZ()-alignmentData.probeSize;
		
		if(Vrui::getMainPipe()!=0)
			{
			/* Send the surface elevation to the slaves: */
			Vrui::getMainPipe()->write<Scalar>(surfaceZ);
			}
		}
	else
		{
		/* Read the surface elevation from the master: */
		Vrui::getMainPipe()->read<Scalar>(surfaceZ);
		}
	
	/* Move the frame's base point: */
	base[2]=surfaceZ;
	
	/* Align the frame with the terrain's x and y directions: */
	alignmentData.surfaceFrame=Vrui::NavTransform(base-Vrui::Point::origin,Vrui::Rotation::identity,alignmentData.surfaceFrame.getScaling());
	}

void LidarViewer::centerDisplayCallback(Misc::CallbackData* cbData)
	{
	/* Initialize the navigation transformation: */
	Vrui::setNavigationTransformation(octrees[0]->getDomainCenter(),octrees[0]->getDomainRadius(),Vrui::Vector(0,0,1));
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
	
	/* Update the other selector mode radio box: */
	switch(defaultSelectorMode)
		{
		case SelectorLocator::ADD:
			mainMenuSelectorModes->setSelectedToggle(0);
			interactionDialogSelectorModes->setSelectedToggle(0);
			break;
		
		case SelectorLocator::SUBTRACT:
			mainMenuSelectorModes->setSelectedToggle(1);
			interactionDialogSelectorModes->setSelectedToggle(1);
			break;
		}
	
	if(overrideTools)
		{
		/* Apply current tool settings to all existing tools: */
		for(LocatorList::iterator lIt=locators.begin();lIt!=locators.end();++lIt)
			(*lIt)->updateSettings();
		}
	}

void LidarViewer::classifySelectionCallback(Misc::CallbackData* cbData)
	{
	for(int i=0;i<numOctrees;++i)
		{
		/* Create a point classification functor: */
		// PointClassifier pc(octrees[i],neighborhoodSize);
		RidgeFinder pc(octrees[i],neighborhoodSize);
		
		/* Classify all selected points: */
		octrees[i]->colorSelectedPoints(pc);
		}
	}

void LidarViewer::saveSelectionCallback(Misc::CallbackData* cbData)
	{
	/* Create the default file name and file name filter: */
	const char* fileName;
	const char* filter;
	if(octrees[0]->hasNormalVectors())
		{
		fileName="SelectedPoints.xyzuvwrgb";
		filter=".xyzuvwrgb";
		}
	else
		{
		fileName="SelectedPoints.xyzrgb";
		filter=".xyzrgb";
		}
	
	try
		{
		/* Open the data directory if it doesn't already exist: */
		if(dataDirectory==0)
			dataDirectory=Vrui::openDirectory(".");
		
		/* Create a uniquely-named selection file in the data directory: */
		std::string selectionFileName=dataDirectory->createNumberedFileName(fileName,4);
		
		/* Create a file selection dialog to select an selection file name: */
		Misc::SelfDestructPointer<GLMotif::FileSelectionDialog> saveSelectionDialog(new GLMotif::FileSelectionDialog(Vrui::getWidgetManager(),"Save Selection...",dataDirectory,selectionFileName.c_str(),filter));
		saveSelectionDialog->getOKCallbacks().add(this,&LidarViewer::saveSelectionOKCallback);
		saveSelectionDialog->deleteOnCancel();
		
		/* Show the file selection dialog: */
		Vrui::popupPrimaryWidget(saveSelectionDialog.releaseTarget());
		}
	catch(std::runtime_error err)
		{
		/* Show an error message: */
		Vrui::showErrorMessage("Save Selection...",Misc::printStdErrMsg("Could not save selection due to exception %s",err.what()));
		}
	}

void LidarViewer::saveSelectionOKCallback(GLMotif::FileSelectionDialog::OKCallbackData* cbData)
	{
	bool success=false;
	
	if(Vrui::isMaster())
		{
		try
			{
			/* Create a selection saver functor: */
			LidarSelectionSaver lss(cbData->getSelectedPath().c_str(),offsets);
			if(octrees[0]->hasNormalVectors())
				{
				/* Save all selected points: */
				octrees[0]->processSelectedPointsWithNormals(lss);
				}
			else
				{
				/* Save all selected points: */
				octrees[0]->processSelectedPoints(lss);
				}
			
			if(Vrui::getMainPipe()!=0)
				{
				/* Send success flag to the slaves: */
				Vrui::getMainPipe()->write<unsigned char>(0);
				}
			
			success=true;
			}
		catch(std::runtime_error err)
			{
			if(Vrui::getMainPipe()!=0)
				{
				/* Send error flag and error message to the slaves: */
				Vrui::getMainPipe()->write<unsigned char>(1);
				Misc::writeCString(err.what(),*Vrui::getMainPipe());
				}
			
			/* Show an error message: */
			Vrui::showErrorMessage("Save Selection...",Misc::printStdErrMsg("Could not write selection to file %s due to exception %s",cbData->getSelectedPath().c_str(),err.what()));
			}
		}
	else
		{
		/* Read success code from the master: */
		if(!(success=Vrui::getMainPipe()->read<unsigned char>()==0))
			{
			/* Read the error message and show an error dialog: */
			char* what=Misc::readCString(*Vrui::getMainPipe());
			Vrui::showErrorMessage("Save Selection...",Misc::printStdErrMsg("Could not write selection to file %s due to exception %s",cbData->getSelectedPath().c_str(),what));
			delete[] what;
			}
		}
	
	if(success)
		{
		/* Remember the current directory for next time: */
		dataDirectory=cbData->selectedDirectory;
		}
	
	/* Destroy the file selection dialog: */
	cbData->fileSelectionDialog->close();
	}

void LidarViewer::clearSelectionCallback(Misc::CallbackData* cbData)
	{
	/* Clear the point selection: */
	for(int i=0;i<numOctrees;++i)
		octrees[i]->clearSelection();
	}

void LidarViewer::extractPlaneCallback(Misc::CallbackData* cbData)
	{
	try
		{
		PlanePrimitive* primitive;
		if(Vrui::isMaster())
			primitive=new PlanePrimitive(octrees[0],Primitive::Vector(offsets),extractorPipe);
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

void LidarViewer::extractBruntonCallback(Misc::CallbackData* cbData)
	{
	try
		{
		BruntonPrimitive* primitive;
		if(Vrui::isMaster())
			primitive=new BruntonPrimitive(octrees[0],Primitive::Vector(offsets),extractorPipe);
		else
			primitive=new BruntonPrimitive(extractorPipe);
		
		/* Store the primitive: */
		lastPickedPrimitive=addPrimitive(primitive);
		
		/* Update the texture generation plane: */
		texturePlane=primitive->getPlane();
		texturePlane.normalize();
		}
	catch(std::runtime_error err)
		{
		if(Vrui::isMaster())
			std::cerr<<"Did not extract brunton due to exception "<<err.what()<<std::endl;
		}
	}

void LidarViewer::extractLineCallback(Misc::CallbackData* cbData)
	{
	try
		{
		LinePrimitive* primitive;
		if(Vrui::isMaster())
			primitive=new LinePrimitive(octrees[0],Primitive::Vector(offsets),extractorPipe);
		else
			primitive=new LinePrimitive(extractorPipe);
		
		/* Store the primitive: */
		lastPickedPrimitive=addPrimitive(primitive);
		}
	catch(std::runtime_error err)
		{
		if(Vrui::isMaster())
			std::cerr<<"Did not extract line due to exception "<<err.what()<<std::endl;
		}
	}

void LidarViewer::extractSphereCallback(Misc::CallbackData* cbData)
	{
	try
		{
		SpherePrimitive* primitive;
		if(Vrui::isMaster())
			primitive=new SpherePrimitive(octrees[0],Primitive::Vector(offsets),extractorPipe);
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
			primitive=new CylinderPrimitive(octrees[0],Primitive::Vector(offsets),extractorPipe);
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
				primitive=new LinePrimitive(planes[0],planes[1],Primitive::Vector(offsets),extractorPipe);
			else
				primitive=new LinePrimitive(extractorPipe);
			}
		else if(planes.size()==3&&lines.empty()&&points.empty())
			{
			/* Create a point by intersecting three planes: */
			if(Vrui::isMaster())
				primitive=new PointPrimitive(planes[0],planes[1],planes[2],Primitive::Vector(offsets),extractorPipe);
			else
				primitive=new PointPrimitive(extractorPipe);
			}
		else if(planes.size()==1&&lines.size()==1&&points.empty())
			{
			/* Create a point by intersecting a plane and a line: */
			if(Vrui::isMaster())
				primitive=new PointPrimitive(planes[0],lines[0],Primitive::Vector(offsets),extractorPipe);
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
	try
		{
		/* Open the data directory if it doesn't already exist: */
		if(dataDirectory==0)
			dataDirectory=Vrui::openDirectory(".");
		
		/* Create a file selection dialog to select a primitive file: */
		Misc::SelfDestructPointer<GLMotif::FileSelectionDialog> loadPrimitivesDialog(new GLMotif::FileSelectionDialog(Vrui::getWidgetManager(),"Load Primitives...",dataDirectory,".dat"));
		loadPrimitivesDialog->getOKCallbacks().add(this,&LidarViewer::loadPrimitivesOKCallback);
		loadPrimitivesDialog->deleteOnCancel();
		
		/* Show the file selection dialog: */
		Vrui::popupPrimaryWidget(loadPrimitivesDialog.releaseTarget());
		}
	catch(std::runtime_error err)
		{
		/* Show an error message: */
		Vrui::showErrorMessage("Load Primitives...",Misc::printStdErrMsg("Could not load primitives due to exception %s",err.what()));
		}
	}
 
void LidarViewer::loadPrimitivesOKCallback(GLMotif::FileSelectionDialog::OKCallbackData* cbData)
	{
	try
		{
		/* Open the primitive file: */
		IO::FilePtr primitiveFile(cbData->selectedDirectory->openFile(cbData->selectedFileName));
		primitiveFile->setEndianness(Misc::LittleEndian);
		
		/* Read the file header: */
		char header[40];
		primitiveFile->read<char>(header,sizeof(header));
		if(strcmp(header,"LidarViewer primitive file v1.2       \n")!=0)
			Misc::throwStdErr("File is not a valid version 1.2 primitive file");
		
		/* Read all primitives in the file: */
		while(!primitiveFile->eof())
			{
			/* Read the primitive type: */
			int primitiveType=primitiveFile->read<int>();
			
			/* Create a primitive of the appropriate type: */
			Primitive* newPrimitive=0;
			switch(primitiveType)
				{
				case 0:
					{
					PlanePrimitive* newPlane=new PlanePrimitive(*primitiveFile,-Primitive::Vector(offsets));
					newPrimitive=newPlane;
					
					/* Update the texture generation plane: */
					texturePlane=newPlane->getPlane();
					texturePlane.normalize();
					break;
					}
				
				case 1:
					newPrimitive=new SpherePrimitive(*primitiveFile,-Primitive::Vector(offsets));
					break;
				
				case 2:
					newPrimitive=new CylinderPrimitive(*primitiveFile,-Primitive::Vector(offsets));
					break;
				
				case 3:
					newPrimitive=new LinePrimitive(*primitiveFile,-Primitive::Vector(offsets));
					break;
				
				case 4:
					newPrimitive=new PointPrimitive(*primitiveFile,-Primitive::Vector(offsets));
					break;
				
				default:
					Misc::throwStdErr("Unknown primitive type %d",primitiveType);
				}
			
			/* Store the primitive: */
			lastPickedPrimitive=addPrimitive(newPrimitive);
			}
		
		/* Remember the current directory for next time: */
		dataDirectory=cbData->selectedDirectory;
		}
	catch(std::runtime_error err)
		{
		/* Show an error message: */
		Vrui::showErrorMessage("Load Primitives...",Misc::printStdErrMsg("Could not load primitives from file %s due to exception %s",cbData->getSelectedPath().c_str(),err.what()));
		}
	
	/* Close the file selection dialog: */
	cbData->fileSelectionDialog->close();
	}

void LidarViewer::savePrimitivesCallback(Misc::CallbackData* cbData)
	{
	try
		{
		/* Open the data directory if it doesn't already exist: */
		if(dataDirectory==0)
			dataDirectory=Vrui::openDirectory(".");

		/* Create a uniquely-named primitive file name in the data directory: */
		std::string primitiveFileName=dataDirectory->createNumberedFileName("SavedPrimitives.dat",4);

		/* Create a file selection dialog to select a primitive file: */
		Misc::SelfDestructPointer<GLMotif::FileSelectionDialog> savePrimitivesDialog(new GLMotif::FileSelectionDialog(Vrui::getWidgetManager(),"Save Primitives...",dataDirectory,primitiveFileName.c_str(),".dat"));
		savePrimitivesDialog->getOKCallbacks().add(this,&LidarViewer::savePrimitivesOKCallback);
		savePrimitivesDialog->deleteOnCancel();
		
		/* Show the file selection dialog: */
		Vrui::popupPrimaryWidget(savePrimitivesDialog.releaseTarget());
		}
	catch(std::runtime_error err)
		{
		/* Show an error message: */
		Vrui::showErrorMessage("Save Primitives...",Misc::printStdErrMsg("Could not save primitives due to exception %s",err.what()));
		}
	}

void LidarViewer::savePrimitivesOKCallback(GLMotif::FileSelectionDialog::OKCallbackData* cbData)
	{
	try
		{
		/* Open the primitive file: */
		IO::FilePtr primitiveFile(cbData->selectedDirectory->openFile(cbData->selectedFileName,IO::File::WriteOnly));
		primitiveFile->setEndianness(Misc::LittleEndian);
		
		/* Write the file header: */
		char header[40];
		snprintf(header,sizeof(header),"LidarViewer primitive file v1.2       \n");
		primitiveFile->write<char>(header,sizeof(header));
		
		/* Write all primitives to the file: */
		for(std::vector<Primitive*>::const_iterator pIt=primitives.begin();pIt!=primitives.end();++pIt)
			{
			/* Determine and write the primitive type: */
			if(dynamic_cast<const PlanePrimitive*>(*pIt)!=0)
				primitiveFile->write<int>(0);
			else if(dynamic_cast<const SpherePrimitive*>(*pIt)!=0)
				primitiveFile->write<int>(1);
			else if(dynamic_cast<const CylinderPrimitive*>(*pIt)!=0)
				primitiveFile->write<int>(2);
			else if(dynamic_cast<const LinePrimitive*>(*pIt)!=0)
				primitiveFile->write<int>(3);
			else if(dynamic_cast<const PointPrimitive*>(*pIt)!=0)
				primitiveFile->write<int>(4);
			
			/* Write the primitive: */
			(*pIt)->write(*primitiveFile,Primitive::Vector(offsets));
			}
		}
	catch(std::runtime_error err)
		{
		/* Show an error message: */
		Vrui::showErrorMessage("Save Primitives...",Misc::printStdErrMsg("Could not write primitives to file %s due to exception %s",cbData->getSelectedPath().c_str(),err.what()));
		}
	
	/* Close the file selection dialog: */
	cbData->fileSelectionDialog->close();
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

void LidarViewer::showOctreeDialogCallback(Misc::CallbackData* cbData)
	{
	if(!Vrui::getWidgetManager()->isManaged(octreeDialog))
		{
		/* Open the octree selection dialog: */
		Vrui::popupPrimaryWidget(octreeDialog);
		}
	}

void LidarViewer::octreeSelectionCallback(GLMotif::ToggleButton::ValueChangedCallbackData* cbData,const int& octreeIndex)
	{
	/* Update the octree visibility flag: */
	showOctrees[octreeIndex]=cbData->set;
	}

void LidarViewer::showRenderDialogCallback(Misc::CallbackData* cbData)
	{
	if(!Vrui::getWidgetManager()->isManaged(renderDialog))
		{
		/* Select the first page (LOD): */
		static_cast<GLMotif::Pager*>(renderDialog->getFirstChild())->setCurrentChildIndex(0);
		
		/* Open the render dialog: */
		Vrui::popupPrimaryWidget(renderDialog);
		}
	}

void LidarViewer::renderQualitySliderCallback(GLMotif::TextFieldSlider::ValueChangedCallbackData* cbData)
	{
	/* Get the new render quality: */
	renderQuality=Scalar(cbData->value);
	
	/* Set all octrees' render quality: */
	for(int i=0;i<numOctrees;++i)
		octrees[i]->setRenderQuality(renderQuality);
	}

void LidarViewer::fncWeightSliderCallback(GLMotif::TextFieldSlider::ValueChangedCallbackData* cbData)
	{
	/* Get the new focus+context weight: */
	fncWeight=Scalar(cbData->value);
	}

void LidarViewer::pointSizeSliderCallback(GLMotif::TextFieldSlider::ValueChangedCallbackData* cbData)
	{
	/* Get the new point size: */
	pointSize=cbData->value;
	}

void LidarViewer::enableLightingCallback(GLMotif::ToggleButton::ValueChangedCallbackData* cbData)
	{
	pointBasedLighting=cbData->set;
	}

void LidarViewer::usePointColorsCallback(GLMotif::ToggleButton::ValueChangedCallbackData* cbData)
	{
	usePointColors=cbData->set;
	}

void LidarViewer::useSplattingCallback(GLMotif::ToggleButton::ValueChangedCallbackData* cbData)
	{
	useSplatting=cbData->set;
	}

void LidarViewer::splatSizeSliderCallback(GLMotif::TextFieldSlider::ValueChangedCallbackData* cbData)
	{
	/* Update the splat size: */
	splatSize=cbData->value;
	}

void LidarViewer::enableSunCallback(GLMotif::ToggleButton::ValueChangedCallbackData* cbData)
	{
	/* Enable the sun light source: */
	setEnableSun(cbData->set);
	}

void LidarViewer::sunAzimuthSliderCallback(GLMotif::TextFieldSlider::ValueChangedCallbackData* cbData)
	{
	/* Update the sun azimuth angle: */
	sunAzimuth=Vrui::Scalar(cbData->value);
	
	/* Update the sun light source: */
	updateSun();
	}

void LidarViewer::sunElevationSliderCallback(GLMotif::TextFieldSlider::ValueChangedCallbackData* cbData)
	{
	/* Update the sun elevation angle: */
	sunElevation=Vrui::Scalar(cbData->value);
	
	/* Update the sun light source: */
	updateSun();
	}

void LidarViewer::enableTexturePlaneCallback(GLMotif::ToggleButton::ValueChangedCallbackData* cbData)
	{
	useTexturePlane=cbData->set;
	}

void LidarViewer::texturePlaneScaleSliderCallback(GLMotif::TextFieldSlider::ValueChangedCallbackData* cbData)
	{
	/* Get the new texture plane scaling factor: */
	texturePlaneScale=cbData->value;
	}

#ifdef LIDARVIEWER_VISUALIZE_WATER

void LidarViewer::texturePlaneOffsetSliderCallback(GLMotif::TextFieldSlider::ValueChangedCallbackData* cbData)
	{
	/* Get the new texture plane offset: */
	texturePlaneOffset=cbData->value;
	}

#endif

void LidarViewer::distanceExaggerationSliderCallback(GLMotif::TextFieldSlider::ValueChangedCallbackData* cbData)
	{
	/* Get new plane distance exaggeration factor: */
	planeDistanceExaggeration=cbData->value;
	
	/* Update the affine coordinate transformer to reflect the new exaggeration value: */
	Vrui::Vector tn=Vrui::Vector(texturePlane.getNormal());
	Vrui::Vector fTrans=tn*(Vrui::Scalar(texturePlane.getOffset())/tn.sqr());
	Vrui::Rotation fRot=Vrui::Rotation::rotateFromTo(Vrui::Vector(0,0,1),tn);
	Vrui::ATransform newTransform=Vrui::ATransform::translate(fTrans);
	newTransform*=Vrui::ATransform::rotate(fRot);
	newTransform*=Vrui::ATransform::scale(Vrui::ATransform::Scale(1.0,1.0,planeDistanceExaggeration));
	newTransform*=Vrui::ATransform::rotate(Geometry::invert(fRot));
	newTransform*=Vrui::ATransform::translate(-(fTrans+Vrui::Vector(offsets)));
	coordTransform->setTransform(newTransform);
	}

void LidarViewer::materialCallback(GLMotif::MaterialEditor::ValueChangedCallbackData* cbData)
	{
	/* Set the surface material properties: */
	surfaceMaterial=cbData->newMaterial;
	}

void LidarViewer::showInteractionDialogCallback(Misc::CallbackData* cbData)
	{
	if(!Vrui::getWidgetManager()->isManaged(interactionDialog))
		{
		/* Open the interaction dialog: */
		Vrui::popupPrimaryWidget(interactionDialog);
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

void LidarViewer::brushSizeSliderCallback(GLMotif::TextFieldSlider::ValueChangedCallbackData* cbData)
	{
	/* Get the new brush size: */
	brushSize=cbData->value;
	
	if(overrideTools)
		{
		/* Apply current tool settings to all existing tools: */
		for(LocatorList::iterator lIt=locators.begin();lIt!=locators.end();++lIt)
			(*lIt)->updateSettings();
		}
	}

void LidarViewer::neighborhoodSizeSliderCallback(GLMotif::TextFieldSlider::ValueChangedCallbackData* cbData)
	{
	/* Get the new neighborhood size: */
	neighborhoodSize=Scalar(cbData->value);
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

