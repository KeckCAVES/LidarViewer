/***********************************************************************
PointBasedLightingShader - Class to maintain a GLSL point-based lighting
shader that tracks the current OpenGL lighting state.
Copyright (c) 2008-2012 Oliver Kreylos

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

#ifndef POINTBASEDLIGHTINGSHADER_INCLUDED
#define POINTBASEDLIGHTINGSHADER_INCLUDED

#include <Geometry/Plane.h>
#include <GL/gl.h>
#include <GL/Extensions/GLARBShaderObjects.h>

/* Forward declarations: */
class GLContextData;

class PointBasedLightingShader
	{
	/* Embedded classes: */
	public:
	typedef Geometry::Plane<double,3> Plane; // Type for texture-mapping planes
	
	/* Elements: */
	private:
	GLContextData& contextData; // The OpenGL context with which this shader is associated
	unsigned int lightStateVersion; // Version of light tracker's state reflected in the current shader program
	unsigned int clipPlaneStateVersion; // Version of clip plane tracker's state reflected in the current shader program
	unsigned int shaderSettingsVersion; // Version of other shader settings reflected in the current shader program
	unsigned int settingsVersion; // Version of other shader settings
	bool usePlaneDistance; // Flag whether the point renderer uses distance from a user-defined plane to texture-map points
	bool usePointColors; // Flag whether the point renderer uses point colors as ambient and diffuse color
	GLhandleARB vertexShader,fragmentShader; // Handle for the vertex and fragment shaders
	GLhandleARB programObject; // Handle for the linked program object
	int planeDistancePlaneLocation; // Location of distance plane equation uniform variable
	int planeDistanceMapLocation; // Location of distance plane texture map uniform variable
	
	/* Private methods: */
	void compileShader(void); // Recompiles the point-based lighting shader based on the current states of all OpenGL light sources and clipping planes
	
	/* Constructors and destructors: */
	public:
	PointBasedLightingShader(GLContextData& sContextData);
	~PointBasedLightingShader(void);
	
	/* Methods: */
	void setUsePlaneDistance(bool newUsePlaneDistance); // Sets the plane distance texturing flag
	void setUsePointColors(bool newUsePointColors); // Sets the point coloring flag
	void enable(void); // Enables point-based lighting in the current OpenGL context
	void setDistancePlane(int textureUnit,const Plane& distancePlane,double distancePlaneScale) const; // Sets the plane distance texturing plane equation and texture unti index
	void disable(void); // Disables point-based lighting in the current OpenGL context
	};

#endif
