/***********************************************************************
PointBasedLightingShader - Class to maintain a GLSL point-based lighting
shader that tracks the current OpenGL lighting state.
Copyright (c) 2008 Oliver Kreylos

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
#include <stdio.h>
#include <iostream>
#include <Misc/ThrowStdErr.h>
#include <GL/gl.h>
#include <GL/GLExtensionManager.h>
#include <GL/Extensions/GLARBShaderObjects.h>
#include <GL/Extensions/GLARBVertexShader.h>
#include <GL/Extensions/GLARBFragmentShader.h>

#include "PointBasedLightingShader.h"

/*************************************************
Static elements of class PointBasedLightingShader:
*************************************************/

const char* PointBasedLightingShader::applyLightTemplate=
	"\
	vec4 applyLight<lightIndex>(const vec4 vertexEc,const vec3 normalEc,const vec4 ambient,const vec4 diffuse,const vec4 specular)\n\
		{\n\
		vec4 color;\n\
		vec3 lightDirEc;\n\
		float nl;\n\
		\n\
		/* Calculate per-source ambient light term: */\n\
		color=gl_LightSource[<lightIndex>].ambient*ambient;\n\
		\n\
		/* Compute the light direction (works both for directional and point lights): */\n\
		lightDirEc=gl_LightSource[<lightIndex>].position.xyz*vertexEc.w-vertexEc.xyz*gl_LightSource[<lightIndex>].position.w;\n\
		lightDirEc=normalize(lightDirEc);\n\
		\n\
		/* Compute the diffuse lighting angle: */\n\
		nl=dot(normalEc,lightDirEc);\n\
		if(nl>0.0)\n\
			{\n\
			vec3 eyeDirEc;\n\
			float nhv;\n\
			\n\
			/* Calculate per-source diffuse light term: */\n\
			color+=(gl_LightSource[<lightIndex>].diffuse*diffuse)*nl;\n\
			\n\
			/* Compute the eye direction: */\n\
			eyeDirEc=normalize(-vertexEc.xyz);\n\
			\n\
			/* Compute the specular lighting angle: */\n\
			nhv=max(dot(normalEc,normalize(eyeDirEc+lightDirEc)),0.0);\n\
			\n\
			/* Calculate per-source specular lighting term: */\n\
			color+=(gl_LightSource[<lightIndex>].specular*specular)*pow(nhv,gl_FrontMaterial.shininess);\n\
			}\n\
		\n\
		/* Return the result color: */\n\
		return color;\n\
		}\n\
	\n";

const char* PointBasedLightingShader::applyAttenuatedLightTemplate=
	"\
	vec4 applyLight<lightIndex>(const vec4 vertexEc,const vec3 normalEc,const vec4 ambient,const vec4 diffuse,const vec4 specular)\n\
		{\n\
		vec4 color;\n\
		vec3 lightDirEc;\n\
		float lightDist,nl,att;\n\
		\n\
		/* Calculate per-source ambient light term: */\n\
		color=gl_LightSource[<lightIndex>].ambient*ambient;\n\
		\n\
		/* Compute the light direction (works both for directional and point lights): */\n\
		lightDirEc=gl_LightSource[<lightIndex>].position.xyz*vertexEc.w-vertexEc.xyz*gl_LightSource[<lightIndex>].position.w;\n\
		lightDist=length(lightDirEc);\n\
		lightDirEc=normalize(lightDirEc);\n\
		\n\
		/* Compute the diffuse lighting angle: */\n\
		nl=dot(normalEc,lightDirEc);\n\
		if(nl>0.0)\n\
			{\n\
			vec3 eyeDirEc;\n\
			float nhv;\n\
			\n\
			/* Calculate per-source diffuse light term: */\n\
			color+=(gl_LightSource[<lightIndex>].diffuse*diffuse)*nl;\n\
			\n\
			/* Compute the eye direction: */\n\
			eyeDirEc=normalize(-vertexEc.xyz);\n\
			\n\
			/* Compute the specular lighting angle: */\n\
			nhv=max(dot(normalEc,normalize(eyeDirEc+lightDirEc)),0.0);\n\
			\n\
			/* Calculate per-source specular lighting term: */\n\
			color+=(gl_LightSource[<lightIndex>].specular*specular)*pow(nhv,gl_FrontMaterial.shininess);\n\
			}\n\
		\n\
		/* Attenuate the per-source light terms: */\n\
		att=(gl_LightSource[<lightIndex>].quadraticAttenuation*lightDist+gl_LightSource[<lightIndex>].linearAttenuation)*lightDist+gl_LightSource[<lightIndex>].constantAttenuation;\n\
		color*=1.0/att;\n\
		\n\
		/* Return the result color: */\n\
		return color;\n\
		}\n\
	\n";

const char* PointBasedLightingShader::applySpotLightTemplate=
	"\
	vec4 applyLight<lightIndex>(const vec4 vertexEc,const vec3 normalEc,const vec4 ambient,const vec4 diffuse,const vec4 specular)\n\
		{\n\
		vec4 color;\n\
		vec3 lightDirEc;\n\
		float sl,nl,att;\n\
		\n\
		/* Calculate per-source ambient light term: */\n\
		color=gl_LightSource[<lightIndex>].ambient*ambient;\n\
		\n\
		/* Compute the light direction (works both for directional and point lights): */\n\
		lightDirEc=gl_LightSource[<lightIndex>].position.xyz*vertexEc.w-vertexEc.xyz*gl_LightSource[<lightIndex>].position.w;\n\
		lightDirEc=normalize(lightDirEc);\n\
		\n\
		/* Calculate the spot light angle: */\n\
		sl=-dot(lightDirEc,normalize(gl_LightSource[<lightIndex>].spotDirection));\n\
		\n\
		/* Check if the point is inside the spot light's cone: */\n\
		if(sl>=gl_LightSource[<lightIndex>].spotCosCutoff)\n\
			{\n\
			/* Compute the diffuse lighting angle: */\n\
			nl=dot(normalEc,lightDirEc);\n\
			if(nl>0.0)\n\
				{\n\
				vec3 eyeDirEc;\n\
				float nhv;\n\
				\n\
				/* Calculate per-source diffuse light term: */\n\
				color+=(gl_LightSource[<lightIndex>].diffuse*diffuse)*nl;\n\
				\n\
				/* Compute the eye direction: */\n\
				eyeDirEc=normalize(-vertexEc.xyz);\n\
				\n\
				/* Compute the specular lighting angle: */\n\
				nhv=max(dot(normalEc,normalize(eyeDirEc+lightDirEc)),0.0);\n\
				\n\
				/* Calculate per-source specular lighting term: */\n\
				color+=(gl_LightSource[<lightIndex>].specular*specular)*pow(nhv,gl_FrontMaterial.shininess);\n\
				}\n\
			\n\
			/* Calculate the spot light attenuation: */\n\
			att=pow(sl,gl_LightSource[<lightIndex>].spotExponent);\n\
			color*=att;\n\
			}\n\
		\n\
		/* Return the result color: */\n\
		return color;\n\
		}\n\
	\n";

const char* PointBasedLightingShader::applyAttenuatedSpotLightTemplate=
	"\
	vec4 applyLight<lightIndex>(const vec4 vertexEc,const vec3 normalEc,const vec4 ambient,const vec4 diffuse,const vec4 specular)\n\
		{\n\
		vec4 color;\n\
		vec3 lightDirEc;\n\
		float sl,nl,att;\n\
		\n\
		/* Calculate per-source ambient light term: */\n\
		color=gl_LightSource[<lightIndex>].ambient*ambient;\n\
		\n\
		/* Compute the light direction (works both for directional and point lights): */\n\
		lightDirEc=gl_LightSource[<lightIndex>].position.xyz*vertexEc.w-vertexEc.xyz*gl_LightSource[<lightIndex>].position.w;\n\
		lightDirEc=normalize(lightDirEc);\n\
		\n\
		/* Calculate the spot light angle: */\n\
		sl=-dot(lightDirEc,normalize(gl_LightSource[<lightIndex>].spotDirection))\n\
		\n\
		/* Check if the point is inside the spot light's cone: */\n\
		if(sl>=gl_LightSource[<lightIndex>].spotCosCutoff)\n\
			{\n\
			/* Compute the diffuse lighting angle: */\n\
			nl=dot(normalEc,lightDirEc);\n\
			if(nl>0.0)\n\
				{\n\
				vec3 eyeDirEc;\n\
				float nhv;\n\
				\n\
				/* Calculate per-source diffuse light term: */\n\
				color+=(gl_LightSource[<lightIndex>].diffuse*diffuse)*nl;\n\
				\n\
				/* Compute the eye direction: */\n\
				eyeDirEc=normalize(-vertexEc.xyz);\n\
				\n\
				/* Compute the specular lighting angle: */\n\
				nhv=max(dot(normalEc,normalize(eyeDirEc+lightDirEc)),0.0);\n\
				\n\
				/* Calculate per-source specular lighting term: */\n\
				color+=(gl_LightSource[<lightIndex>].specular*specular)*pow(nhv,gl_FrontMaterial.shininess);\n\
				}\n\
			\n\
			/* Calculate the spot light attenuation: */\n\
			att=pow(sl,gl_LightSource[<lightIndex>].spotExponent);\n\
			color*=att;\n\
			}\n\
		\n\
		/* Attenuate the per-source light terms: */\n\
		att=(gl_LightSource[<lightIndex>].quadraticAttenuation*lightDist+gl_LightSource[<lightIndex>].linearAttenuation)*lightDist+gl_LightSource[<lightIndex>].constantAttenuation;\n\
		color*=1.0/att;\n\
		\n\
		/* Return the result color: */\n\
		return color;\n\
		}\n\
	\n";

/*****************************************
Methods of class PointBasedLightingShader:
*****************************************/

std::string PointBasedLightingShader::createApplyLightFunction(const char* functionTemplate,int lightIndex) const
	{
	std::string result;
	
	/* Create the light index string: */
	char index[10];
	snprintf(index,sizeof(index),"%d",lightIndex);
	
	/* Replace all occurrences of <lightIndex> in the template string with the light index string: */
	const char* match="<lightIndex>";
	const char* matchStart;
	int matchLen=0;
	for(const char* tPtr=functionTemplate;*tPtr!='\0';++tPtr)
		{
		if(matchLen==0)
			{
			if(*tPtr=='<')
				{
				matchStart=tPtr;
				matchLen=1;
				}
			else
				result.push_back(*tPtr);
			}
		else if(matchLen<12)
			{
			if(*tPtr==match[matchLen])
				{
				++matchLen;
				if(matchLen==strlen(match))
					{
					result.append(index);
					matchLen=0;
					}
				}
			else
				{
				for(const char* cPtr=matchStart;cPtr!=tPtr;++cPtr)
					result.push_back(*cPtr);
				matchLen=0;
				--tPtr;
				}
			}
		}
	
	return result;
	}

PointBasedLightingShader::PointBasedLightingShader(void)
	:colorMaterial(false),
	 lightStates(0),
	 vertexShader(0),fragmentShader(0),
	 programObject(0)
	{
	/* Determine the maximum number of light sources supported by the local OpenGL: */
	glGetIntegerv(GL_MAX_LIGHTS,&maxNumLights);
	
	/* Initialize the light state array: */
	lightStates=new LightState[maxNumLights];
	updateLightingState();
	
	/* Check for the required OpenGL extensions: */
	if(!GLARBShaderObjects::isSupported())
		Misc::throwStdErr("GLShader::GLShader: GL_ARB_shader_objects not supported");
	if(!GLARBVertexShader::isSupported())
		Misc::throwStdErr("GLShader::GLShader: GL_ARB_vertex_shader not supported");
	if(!GLARBFragmentShader::isSupported())
		Misc::throwStdErr("GLShader::GLShader: GL_ARB_fragment_shader not supported");
	
	/* Initialize the required extensions: */
	GLARBShaderObjects::initExtension();
	GLARBVertexShader::initExtension();
	GLARBFragmentShader::initExtension();
	
	/* Create the vertex and fragment shaders: */
	vertexShader=glCreateShaderObjectARB(GL_VERTEX_SHADER_ARB);
	fragmentShader=glCreateShaderObjectARB(GL_FRAGMENT_SHADER_ARB);
	
	/* Create the program object: */
	programObject=glCreateProgramObjectARB();
	glAttachObjectARB(programObject,vertexShader);
	glAttachObjectARB(programObject,fragmentShader);
	}

PointBasedLightingShader::~PointBasedLightingShader(void)
	{
	glDeleteObjectARB(programObject);
	glDeleteObjectARB(vertexShader);
	glDeleteObjectARB(fragmentShader);
	delete[] lightStates;
	}

bool PointBasedLightingShader::updateLightingState(void)
	{
	bool mustRecompile=false;
	
	/* Check the color material flag: */
	bool newColorMaterial=glIsEnabled(GL_COLOR_MATERIAL);
	if(newColorMaterial!=colorMaterial)
		mustRecompile=true;
	colorMaterial=newColorMaterial;
	
	for(int lightIndex=0;lightIndex<maxNumLights;++lightIndex)
		{
		GLenum light=GL_LIGHT0+lightIndex;
		
		/* Get the light's enabled flag: */
		bool enabled=glIsEnabled(light);
		if(enabled!=lightStates[lightIndex].enabled)
			mustRecompile=true;
		lightStates[lightIndex].enabled=enabled;
		
		if(enabled)
			{
			/* Determine the light's attenuation and spot light state: */
			bool attenuated=false;
			bool spotLight=false;
			
			/* Get the light's position: */
			GLfloat pos[4];
			glGetLightfv(light,GL_POSITION,pos);
			if(pos[3]!=0.0f)
				{
				/* Get the light's attenuation coefficients: */
				GLfloat att[3];
				glGetLightfv(light,GL_CONSTANT_ATTENUATION,&att[0]);
				glGetLightfv(light,GL_LINEAR_ATTENUATION,&att[1]);
				glGetLightfv(light,GL_QUADRATIC_ATTENUATION,&att[2]);
				
				/* Determine whether the light is attenuated: */
				if(att[0]!=1.0f||att[1]!=0.0f||att[2]!=0.0f)
					attenuated=true;
				
				/* Get the light's spot light cutoff angle: */
				GLfloat spotLightCutoff;
				glGetLightfv(light,GL_SPOT_CUTOFF,&spotLightCutoff);
				spotLight=spotLightCutoff<=90.0f;
				}
			
			if(attenuated!=lightStates[lightIndex].attenuated)
				mustRecompile=true;
			lightStates[lightIndex].attenuated=attenuated;
			
			if(spotLight!=lightStates[lightIndex].spotLight)
				mustRecompile=true;
			lightStates[lightIndex].spotLight=spotLight;
			}
		}
	
	return mustRecompile;
	}

void PointBasedLightingShader::compileShader(void)
	{
	std::string vertexShaderFunctions;
	std::string vertexShaderMain;
	
	/* Create the main vertex shader starting boilerplate: */
	vertexShaderMain+=
		"\
		void main()\n\
			{\n\
			vec4 vertexEc;\n\
			vec3 normalEc;\n\
			vec4 ambient,diffuse,specular;\n\
			vec4 color;\n\
			\n\
			/* Compute the vertex position in eye coordinates: */\n\
			vertexEc=gl_ModelViewMatrix*gl_Vertex;\n\
			\n\
			/* Compute the normal vector in eye coordinates: */\n\
			normalEc=normalize(gl_NormalMatrix*gl_Normal);\n\
			\n\
			/* Let the normal vector always point towards the eye: */\n\
			if(dot(normalEc,vertexEc.xyz)>0.0)\n\
				normalEc=-normalEc;\n\
			\n";
	
	/* Get the material components: */
	if(colorMaterial)
		{
		vertexShaderMain+=
			"\
			/* Get the material properties from the current color: */\n\
			ambient=gl_Color;\n\
			diffuse=gl_Color;\n\
			specular=gl_FrontMaterial.specular;\n\
			\n";
		}
	else
		{
		vertexShaderMain+=
			"\
			/* Get the material properties from the material state: */\n\
			ambient=gl_FrontMaterial.ambient;\n\
			diffuse=gl_FrontMaterial.diffuse;\n\
			specular=gl_FrontMaterial.specular;\n\
			\n";
		}
	
	/* Continue the main vertex shader: */
	vertexShaderMain+=
		"\
		/* Calculate global ambient light term: */\n\
		color=gl_LightModel.ambient*ambient;\n\
		\n\
		/* Apply all enabled light sources: */\n";
	
	/* Create light application functions for all enabled light sources: */
	for(int lightIndex=0;lightIndex<maxNumLights;++lightIndex)
		if(lightStates[lightIndex].enabled)
			{
			/* Create the appropriate light application function: */
			if(lightStates[lightIndex].spotLight)
				{
				if(lightStates[lightIndex].attenuated)
					vertexShaderFunctions+=createApplyLightFunction(applyAttenuatedSpotLightTemplate,lightIndex);
				else
					vertexShaderFunctions+=createApplyLightFunction(applySpotLightTemplate,lightIndex);
				}
			else
				{
				if(lightStates[lightIndex].attenuated)
					vertexShaderFunctions+=createApplyLightFunction(applyAttenuatedLightTemplate,lightIndex);
				else
					vertexShaderFunctions+=createApplyLightFunction(applyLightTemplate,lightIndex);
				}
			
			/* Call the light application function from the shader's main function: */
			char call[256];
			snprintf(call,sizeof(call),"\t\t\tcolor+=applyLight%d(vertexEc,normalEc,ambient,diffuse,specular);\n",lightIndex);
			vertexShaderMain+=call;
			}
	
	/* Finish the main function: */
	vertexShaderMain+=
		"\
			\n\
			/* Compute final vertex color: */\n\
			gl_FrontColor=color;\n\
			\n\
			/* Use standard vertex position: */\n\
			gl_Position=ftransform();\n\
			}\n";
	
	/* Compile the vertex shader: */
	std::string vertexShaderSource=vertexShaderFunctions+vertexShaderMain;
	glCompileShaderFromString(vertexShader,vertexShaderSource.c_str());
	
	/* Compile the standard fragment shader: */
	const char* fragmentShaderSource=
		"\
		void main()\n\
			{\n\
			gl_FragColor=gl_Color;\n\
			}\n";
	glCompileShaderFromString(fragmentShader,fragmentShaderSource);
	
	/* Link the program object: */
	glLinkProgramARB(programObject);
	
	/* Check if the program linked successfully: */
	GLint linkStatus;
	glGetObjectParameterivARB(programObject,GL_OBJECT_LINK_STATUS_ARB,&linkStatus);
	if(!linkStatus)
		{
		/* Get some more detailed information: */
		GLcharARB linkLogBuffer[2048];
		GLsizei linkLogSize;
		glGetInfoLogARB(programObject,sizeof(linkLogBuffer),&linkLogSize,linkLogBuffer);
		
		/* Signal an error: */
		Misc::throwStdErr("Error \"%s\" while linking shader program",linkLogBuffer);
		}
	}

void PointBasedLightingShader::enable(void)
	{
	/* Update light states and recompile the shader if necessary: */
	if(updateLightingState())
		compileShader();
	
	/* Enable the shader: */
	glUseProgramObjectARB(programObject);
	}

void PointBasedLightingShader::disable(void)
	{
	/* Disable the shader: */
	glUseProgramObjectARB(0);
	}
