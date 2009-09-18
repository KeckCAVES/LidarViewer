/***********************************************************************
BruntonPrimitive - Class for planes extracted from point clouds, with
additional direct visualization of strike and dip angles.
Copyright (c) 2009 Oliver Kreylos

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

#include "BruntonPrimitive.h"

#include <stdio.h>
#include <Math/Math.h>
#include <Geometry/Point.h>
#include <Geometry/AffineCombiner.h>
#include <SceneGraph/TransformNode.h>
#include <SceneGraph/BillboardNode.h>
#include <SceneGraph/ColorNode.h>
#include <SceneGraph/CoordinateNode.h>
#include <SceneGraph/IndexedLineSetNode.h>
#include <SceneGraph/FontStyleNode.h>
#include <SceneGraph/TextNode.h>
#include <SceneGraph/ShapeNode.h>

#include "SceneGraph.h"

/*********************************
Methods of class BruntonPrimitive:
*********************************/

void BruntonPrimitive::buildBrunton(void)
	{
	/* Create the root node: */
	SceneGraph::TransformNode* rootT=new SceneGraph::TransformNode;
	root=rootT;
	getSceneGraphRoot().children.appendValue(root);
	getSceneGraphRoot().update();
	
	/* Calculate the plane primitive's centroid: */
	Point::AffineCombiner cc;
	for(int i=0;i<4;++i)
		cc.addPoint(getPoint(i));
	Point centroid=cc.getPoint();
	Scalar bruntonScale=Math::mid(Geometry::dist(getPoint(2),getPoint(0)),Geometry::dist(getPoint(3),getPoint(1)));
	
	/* Calculate the plane primitive's dip angle and strike vector: */
	Vector normal=getPlane().getNormal();
	if(normal[2]<Scalar(0))
		normal=-normal;
	normal.normalize();
	Vector strike=normal;
	Scalar dipAngle=Math::acos(strike[2]/Geometry::mag(strike));
	strike[2]=Scalar(0);
	strike.normalize();
	Scalar strikeAngle=Math::atan2(-strike[0],strike[1]);
	
	/* Set the root node's transformation: */
	rootT->translation.setValue(centroid-Point::origin);
	
	/* Create the dip and strike indicator: */
	SceneGraph::ShapeNode* s1=new SceneGraph::ShapeNode;
	rootT->children.appendValue(s1);
	{
	SceneGraph::IndexedLineSetNode* ils=new SceneGraph::IndexedLineSetNode;
	s1->geometry.setValue(ils);
	
	SceneGraph::ColorNode* color=new SceneGraph::ColorNode;
	ils->color.setValue(color);
	color->color.appendValue(SceneGraph::Color(0.0f,0.5f,1.0f));
	color->color.appendValue(SceneGraph::Color(0.0f,1.0f,0.5f));
	color->update();
	
	SceneGraph::CoordinateNode* coord=new SceneGraph::CoordinateNode;
	ils->coord.setValue(coord);
	coord->point.appendValue(Point::origin);
	coord->point.appendValue(Point::origin+normal*bruntonScale);
	coord->point.appendValue(Point::origin+strike*bruntonScale);
	coord->update();
	
	ils->coordIndex.appendValue(0);
	ils->coordIndex.appendValue(1);
	ils->coordIndex.appendValue(-1);
	ils->coordIndex.appendValue(0);
	ils->coordIndex.appendValue(2);
	
	ils->colorPerVertex.setValue(false);
	ils->lineWidth.setValue(3.0f);
	ils->update();
	}
	s1->update();
	
	SceneGraph::ShapeNode* s2=new SceneGraph::ShapeNode;
	rootT->children.appendValue(s2);
	{
	SceneGraph::IndexedLineSetNode* ils=new SceneGraph::IndexedLineSetNode;
	s2->geometry.setValue(ils);
	
	SceneGraph::ColorNode* color=new SceneGraph::ColorNode;
	ils->color.setValue(color);
	color->color.appendValue(SceneGraph::Color(0.0f,0.5f,1.0f));
	color->color.appendValue(SceneGraph::Color(0.0f,0.5f,1.0f));
	color->color.appendValue(SceneGraph::Color(0.0f,1.0f,0.5f));
	color->color.appendValue(SceneGraph::Color(0.0f,1.0f,0.5f));
	color->update();
	
	SceneGraph::CoordinateNode* coord=new SceneGraph::CoordinateNode;
	ils->coord.setValue(coord);
	
	coord->point.appendValue(Point::origin);
	coord->point.appendValue(Point(0,0,bruntonScale));
	coord->point.appendValue(Point(0,bruntonScale,0));
	ils->coordIndex.appendValue(0);
	ils->coordIndex.appendValue(1);
	ils->coordIndex.appendValue(-1);
	for(Scalar a=Scalar(0);a<dipAngle;a+=Math::rad(Scalar(10)))
		{
		ils->coordIndex.appendValue(coord->point.getNumValues());
		coord->point.appendValue(Point::origin+(Vector(0,0,1)*Math::cos(a)+strike*Math::sin(a))*(bruntonScale*Scalar(0.9)));
		}
	ils->coordIndex.appendValue(coord->point.getNumValues());
	coord->point.appendValue(Point::origin+(Vector(0,0,1)*Math::cos(dipAngle)+strike*Math::sin(dipAngle))*(bruntonScale*Scalar(0.9)));
	ils->coordIndex.appendValue(-1);
	
	ils->coordIndex.appendValue(0);
	ils->coordIndex.appendValue(2);
	ils->coordIndex.appendValue(-1);
	Scalar aInc=Math::rad(Scalar(10));
	if(strikeAngle<Scalar(0))
		aInc=-aInc;
	for(Scalar a=Scalar(0);Math::abs(a)<Math::abs(strikeAngle);a+=aInc)
		{
		ils->coordIndex.appendValue(coord->point.getNumValues());
		coord->point.appendValue(Point::origin+Vector(-Math::sin(a),Math::cos(a),0)*(bruntonScale*Scalar(0.9)));
		}
	ils->coordIndex.appendValue(coord->point.getNumValues());
	coord->point.appendValue(Point::origin+Vector(-Math::sin(strikeAngle),Math::cos(strikeAngle),0)*(bruntonScale*Scalar(0.9)));
	
	coord->update();
	
	ils->colorPerVertex.setValue(false);
	ils->lineWidth.setValue(1.0f);
	ils->update();
	}
	s2->update();
	
	SceneGraph::TransformNode* t2=new SceneGraph::TransformNode;
	rootT->children.appendValue(t2);
	t2->translation.setValue((Vector(0,0,1)*Math::cos(Math::div2(dipAngle))+strike*Math::sin(Math::div2(dipAngle)))*bruntonScale);
	{
	SceneGraph::BillboardNode* bb=new SceneGraph::BillboardNode;
	t2->children.appendValue(bb);
	bb->axisOfRotation.setValue(Vector::zero);
	{
	SceneGraph::ShapeNode* s=new SceneGraph::ShapeNode;
	bb->children.appendValue(s);
	
	SceneGraph::TextNode* text=new SceneGraph::TextNode;
	s->geometry.setValue(text);
	
	SceneGraph::FontStyleNode* fs=new SceneGraph::FontStyleNode;
	text->fontStyle.setValue(fs);
	fs->size.setValue(bruntonScale*Scalar(0.25));
	fs->justify.appendValue("MIDDLE");
	fs->justify.appendValue("MIDDLE");
	
	fs->update();
	
	char buffer[40];
	snprintf(buffer,sizeof(buffer),"%.2f",Math::deg(dipAngle));
	text->string.appendValue(buffer);
	text->update();
	
	s->update();
	}
	bb->update();
	}
	t2->update();
	
	SceneGraph::TransformNode* t3=new SceneGraph::TransformNode;
	rootT->children.appendValue(t3);
	t3->translation.setValue(Vector(-Math::sin(Math::div2(strikeAngle)),Math::cos(Math::div2(strikeAngle)),0)*bruntonScale);
	{
	SceneGraph::BillboardNode* bb=new SceneGraph::BillboardNode;
	t3->children.appendValue(bb);
	bb->axisOfRotation.setValue(Vector::zero);
	{
	SceneGraph::ShapeNode* s=new SceneGraph::ShapeNode;
	bb->children.appendValue(s);
	
	SceneGraph::TextNode* text=new SceneGraph::TextNode;
	s->geometry.setValue(text);
	
	SceneGraph::FontStyleNode* fs=new SceneGraph::FontStyleNode;
	text->fontStyle.setValue(fs);
	fs->size.setValue(bruntonScale*Scalar(0.25));
	fs->justify.appendValue("MIDDLE");
	fs->justify.appendValue("MIDDLE");
	
	fs->update();
	
	char buffer[40];
	Scalar sa=-Math::deg(strikeAngle);
	if(sa<Scalar(0))
		sa+=360.0;
	snprintf(buffer,sizeof(buffer),"%.2f",sa);
	text->string.appendValue(buffer);
	text->update();
	
	s->update();
	}
	bb->update();
	}
	t3->update();
	
	rootT->update();
	}

BruntonPrimitive::BruntonPrimitive(const LidarOctree* octree,Comm::MulticastPipe* pipe)
	:PlanePrimitive(octree,pipe)
	{
	buildBrunton();
	}

BruntonPrimitive::BruntonPrimitive(Comm::MulticastPipe* pipe)
	:PlanePrimitive(pipe)
	{
	buildBrunton();
	}

BruntonPrimitive::BruntonPrimitive(Misc::File& file,const Primitive::Vector& translation)
	:PlanePrimitive(file,translation)
	{
	buildBrunton();
	}

BruntonPrimitive::~BruntonPrimitive(void)
	{
	/* Remove the brunton root node from the scene graph: */
	getSceneGraphRoot().children.removeValue(root);
	getSceneGraphRoot().update();
	}
