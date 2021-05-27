#include "Elements.h"
#include "glew.h"
#include <GL/gl.h>
#define _USE_MATH_DEFINES //for using M_PI
#include <math.h>
#include <vector>
float* CrossProduct(float *a, float *b)
{
	float Product[3];

	//Cross product formula 
	Product[0] = (a[1] * b[2]) - (a[2] * b[1]);
	Product[1] = (a[2] * b[0]) - (a[0] * b[2]);
	Product[2] = (a[0] * b[1]) - (a[1] * b[0]);

	return Product;
}
void genCone(float height, float radius, unsigned slices) {
	float angle = 2 * M_PI / slices;
	//for the circular sector
	glBegin(GL_TRIANGLES);
	for (size_t i = 0; i < slices; i++) {
		glNormal3d(0, -1, 0);
		glVertex3f(0.0f, 0.0f, 0.0f);
		glVertex3f(radius * cos(angle * i), 0.0f, -radius * sin(angle * i));
		glVertex3f(radius * cos(angle * (i+1)), 0.0f, -radius * sin(angle * (i+1)));

	}
	glEnd();
	//for the surface
	glBegin(GL_TRIANGLES);
	for (size_t i = 0; i < slices; i++) {
		float a[3], b[3];
		a[0] = radius * cos(angle * i);
		a[1] = -height;
		a[2] = -radius * sin(angle * i);
		b[0] = radius * cos(angle * (i + 1)) - radius * cos(angle * i);
		b[1] = 0.0;
		b[2] = -radius * sin(angle * (i + 1)) + radius * sin(angle * i);
		float* c = CrossProduct(a, b);
		glNormal3d(c[0], c[1], c[2]);
		glVertex3f(0.0f, height, 0.0f);
		glVertex3f(radius * cos(angle * i), 0.0f, -radius * sin(angle * i));
		glVertex3f(radius * cos(angle * (i + 1)), 0.0f, -radius * sin(angle * (i + 1)));
	}
	glEnd();
}
void genCylinder(float height, float radius, size_t slices) {
	float angle = 2 * M_PI / slices;
	//for the circular sector 1
	glBegin(GL_TRIANGLES);
	for (size_t i = 0; i < slices; i++) {
		glNormal3d(0, -1, 0);
		glVertex3f(0.0f, 0.0f, 0.0f);
		glVertex3f(radius * cos(angle * i), 0.0f, -radius * sin(angle * i));
		glVertex3f(radius * cos(angle * (i + 1)), 0.0f, -radius * sin(angle * (i + 1)));
	}
	glEnd();
	//for the circular sector 2
	glBegin(GL_TRIANGLES);
	for (size_t i = 0; i < slices; i++) {
		glNormal3d(0, 1, 0);
		glVertex3f(0.0f, height, 0.0f);
		glVertex3f(radius * cos(angle * i), height, -radius * sin(angle * i));
		glVertex3f(radius * cos(angle * (i + 1)), height, -radius * sin(angle * (i + 1)));
	}
	glEnd();
	//for the surface
	for (size_t i = 0; i < slices; i++) {	
		glBegin(GL_POLYGON);
		float a[3], b[3];
		a[0] = 0;
		a[1] = -height;
		a[2] = 0;
		b[0] = radius * cos(angle * (i + 1)) - radius * cos(angle * i);
		b[1] = 0.0;
		b[2] = -radius * sin(angle * (i + 1)) + radius * sin(angle * i);
		float* c = CrossProduct(a, b);
		glNormal3d(c[0], c[1], c[2]);
		glVertex3f(radius * cos(angle * i), height, -radius * sin(angle * i));
		glVertex3f(radius * cos(angle * i), 0.0f, -radius * sin(angle * i));
		glVertex3f(radius * cos(angle * (i + 1)), 0.0f, -radius * sin(angle * (i + 1)));
		glVertex3f(radius * cos(angle * (i + 1)), height, -radius * sin(angle * (i + 1)));
		glEnd();
	}

}
void genSphere(size_t stackCount, size_t sectorCount, float radius) {
	//theta [-M_PI, MPI] fai [0, M_PI]
	//vertex position (_radius * cos(theta1) * sin(fai), -_radius * cos(fai), -_radius * sin(theta1) * sin(fai));
	double theta_interval = 2 * M_PI / sectorCount;
	double fai_interval = M_PI / stackCount;
	for (size_t i = 0; i < sectorCount; i++)
	{
		glBegin(GL_QUAD_STRIP);
		for (size_t j = 0; j <= stackCount; j++)
		{
			double theta1 = -M_PI + i * theta_interval;
			double theta2 = -M_PI + (i + 1) * theta_interval;
			double fai = 0 + j * fai_interval;
			double u1 = i * 1.0 / sectorCount;
			double v1 = j * 1.0 / stackCount;
			double u2 = (i + 1) * 1.0 / sectorCount;
			double v2 = j * 1.0 / stackCount;

			glTexCoord2f(u1, v1);
			glNormal3d(cos(theta1) * sin(fai), -cos(fai), -sin(theta1) * sin(fai));
			glVertex3d(radius * cos(theta1) * sin(fai), -radius * cos(fai), -radius * sin(theta1) * sin(fai));

			glTexCoord2f(u2, v2);
			glNormal3d(cos(theta2) * sin(fai), -cos(fai), -sin(theta2) * sin(fai));
			glVertex3d(radius * cos(theta2) * sin(fai), -radius * cos(fai), -radius * sin(theta2) * sin(fai));
		}
		glEnd();
	}
}
void genCube(float w, float l, float h) {
	//first face
	glBegin(GL_POLYGON);
	glNormal3d(0, 0, 1);
	glVertex3d(-l/2, h/2, w/2);
	glVertex3d(-l/2, -h/2, w/2);
	glVertex3d(l/2, -h/2, w/2);
	glVertex3d(l/2, h/2, w/2);
	glEnd();
	//second face
	glBegin(GL_POLYGON);
	glNormal3d(1, 0, 0);
	glVertex3d(l / 2, h / 2, w / 2);
	glVertex3d(l / 2, -h / 2, w / 2);
	glVertex3d(l/2, -h/2, -w/2);
	glVertex3d(l/2, h/2, -w/2);
	glEnd();
	//
	glBegin(GL_POLYGON);
	glNormal3d(0, 0, -1);
	glVertex3d(-l/2, -h/2, -w/2);
	glVertex3d(l/2, -h/2, -w/2);
	glVertex3d(l/2,h/2, -w/2);
	glVertex3d(-l/2, h/2, -w/2);
	glEnd();
	//
	glBegin(GL_POLYGON);
	glNormal3d(-1, 0, 0);
	glVertex3d(-l / 2, h / 2, -w / 2);
	glVertex3d(-l/2, -h/2, -w/2);
	glVertex3d(-l/2, -h/2, w/2);
	glVertex3d(-l/2, h/2, w/2);
	glEnd();
	//
	glBegin(GL_POLYGON);
	glNormal3d(0, 1, 0);
	glVertex3d(-l/2, h/2, -w/2);
	glVertex3d(-l/2, h/2, w/2);
	glVertex3d(l/2, h/2, w/2);
	glVertex3d(l/2, h/2, -w/2);
	glEnd();
	//
	glBegin(GL_POLYGON);
	glNormal3d(0, -1, 0);
	glVertex3d(-l/2, -h/2, w/2);
	glVertex3d(-l/2, -h/2, -w/2);
	glVertex3d(l/2, -h/2, -w/2);
	glVertex3d(l/2, -h/2, w/2);
	glEnd();
}
float gau(float x, float sigma) {
	return (1.0 / (sigma*sqrt(2 * M_PI))) * exp((-1.0/2) * (pow(x,2)/pow(sigma,2)));
}
void gauss1d(size_t cell_num, float r, float g, float b, float sigma, float scale) {
	float step = 1.0 / cell_num;
	for (size_t i = 0; i < cell_num; i++) {
		for (size_t j = 0; j < cell_num; j++) {
			glBegin(GL_POLYGON);
			float x = 1.0*i / (cell_num - 1) - 0.5;
			float g_c = gau(x * scale, sigma);
			glColor3f(g_c*r, g_c*g, g_c*b);
			glVertex3d(j*step - 0.5, i*step - 0.5, 0);
			glVertex3d((j+1)*step - 0.5, i*step - 0.5, 0);
			x = 1.0*(i+1) / (cell_num - 1) - 0.5;
			g_c = gau(x * scale, sigma);
			glColor3f(g_c*r, g_c*g, g_c*b);
			glVertex3d((j + 1)*step - 0.5, (i+1)*step - 0.5, 0);
			glVertex3d(j*step - 0.5, (i+1)*step - 0.5, 0);
			glEnd();
		}
	}
}
void exampleArrows1() {
	genCone(1, 1, 200);
	glPushMatrix();
	glTranslatef(0, -1, 0);
	genCylinder(1, 0.5, 200);
	glPopMatrix();
}
void exampleArrows2() {
	genCone(1, 1, 200);
	glPushMatrix();
	glTranslatef(0, -1, 0);
	glScalef(1.0, 2.3, 1.0);
	genSphere(200, 200, 0.5);
	glPopMatrix();
}
void exampleArrows3() {
	genCone(1, 1, 200);
	glPushMatrix();
	glTranslatef(0, -0.5, 0);
	genCube(0.5, 0.5, 1);
	glPopMatrix();
}