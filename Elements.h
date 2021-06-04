#ifndef ELEMENTS_H_
#define ELEMENTS_H_
/**
* generate cone with height,  raidus and slices
*/
void genCone(float height, float radius, size_t slices);
/**
* generate cylinder with height,  raidus and slices
*/
void genCylinder(float height, float radius, size_t slices);
/**
* generate sphere with stack count, sect count and radius
*/
void genSphere(size_t stackCount, size_t sectorCount, float radius);
/*
* generate cube with width, length, and heigth
*/
void genCube(float w, float l, float h);
/*
* 1-d gaussian
*/
void gauss1d(size_t cell_num,float r, float g, float b, float sigma, float scale);
/*
* generate square
*/
void genSquare(size_t cell_num);
/*
* example arrows
*/
void exampleArrows1();
void exampleArrows2();
void exampleArrows3();
#endif // !ELEMENTS_H_
