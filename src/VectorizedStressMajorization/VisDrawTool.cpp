#include "VectorizedStressMajorization.h"



void drawGrid(float gap) {

	glColor4f(0.9, 0.9, 0.9, 1.0);

	glBegin(GL_LINES);
	glVertex2f(0, -1);
	glVertex2f(0, 1);
	glVertex2f(-1, 0);
	glVertex2f(1, 0);
	for (int i = 0; i < 1.0 / gap; i++) {

		//vertical
		glVertex2f(0 - gap*i, -1);
		glVertex2f(0 - gap*i, 1);
		glVertex2f(0 + gap*i, -1);
		glVertex2f(0 + gap*i, 1);

		//horizontal
		glVertex2f(-1, 0 - gap*i);
		glVertex2f(1, 0 - gap*i);
		glVertex2f(-1, 0 + gap*i);
		glVertex2f(1, 0 + gap*i);
	}

	glColor4f(0.0, 0.9, 0.9, 1.0);
	float a = tanf(M_PI / 8);
	glVertex2f(1, a);
	glVertex2f(-1, -a);

	glVertex2f(a, 1);
	glVertex2f(-a, -1);

	glVertex2f(-a, 1);
	glVertex2f(a, -1);

	glVertex2f(-1, a);
	glVertex2f(1, -a);
	glEnd();
}
void drawTextArea(char* str) {
	int n = strlen(str);
	//设置要在屏幕上显示字符的起始位置  
	glRasterPos2i(-1, -1);
	//逐个显示字符串中的每个字符  
	for (int i = 0; i < n; i++)
		glutBitmapCharacter(GLUT_BITMAP_8_BY_13, *(str + i));
}