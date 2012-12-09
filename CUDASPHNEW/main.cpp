#include <Windows.h>
#include <gl\glew.h>
#include <gl\glut.h>
#include <gl\freeglut_ext.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <cuda_gl_interop.h>
#include "particle_simulation.h"
#include "particle_renderer.h"
#include <stdexcept>
#include <iostream>

void display();
void keyboard(unsigned char key, int x, int y);
void init();


ParticleSimulation* gSim;
ParticleRenderer* gRenderer;


int main(int argc, char* argv[]) 
{
    //cudaGLSetGLDevice(0);    
    glutInit(&argc, argv);
    glutInitContextVersion(3, 3);
    glutInitContextFlags(GLUT_FORWARD_COMPATIBLE);
    glutInitContextProfile(GLUT_CORE_PROFILE);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
	//int width = glutGet(GLUT_SCREEN_WIDTH);
	//int height = glutGet(GLUT_SCREEN_HEIGHT);
	glutInitWindowSize(1280, 800);
    glutCreateWindow("SPH CUDA");
	glewExperimental = TRUE;
	glewInit();
    init();
	glutDisplayFunc(display);
	glutKeyboardFunc(keyboard);
	glutMainLoop();

    return 0;
}

void init() 
{
    try {
        gSim = ParticleSimulation::example01();
        gSim->init();
        gSim->bind();

        gRenderer = new ParticleRenderer(*gSim, 1280, 800);
        gRenderer->setCamera(0.0f, 0.3f, 2.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f);
	    gRenderer->setPerspective(60.0f, 1280.0f/800.0f, 0.1f, 100.0f);

        //gSim->saveInfoTable("tralala.txt");
        //gSim->saveParticleInfo("particle_info.txt");

        printf("+++ finished +++\n");

        //delete gRenderer;
        //delete gSim;

    } catch (std::runtime_error e) {
        std::cout << e.what();
    }
}

void display() 
{
    gRenderer->render();
    gSim->advance();
 //   system("pause");
}

void keyboard(unsigned char key, int x, int y)
{

}