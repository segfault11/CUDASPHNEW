#include <Windows.h>
#include <gl\glew.h>
#include <gl\glut.h>
#include <gl\freeglut_ext.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <cuda_gl_interop.h>
#include "particle_simulation.h"
#include "particle_renderer01.h"
#include "particle_renderer02.h"
#include "triangle_mesh.h"
#include "obstacle_renderer.h"
#include "obstacle_grid.h"
#include <stdexcept>
#include <iostream>
#include <exception>
#include "portable_pixmap.h"
#include "sparse_voxel_map.h"
#include "boundary_map.h"
#include "twoscalestate.h"
#include "twoscalestatesub.h"
#include <iomanip>

using namespace std;

void display();
void keyboard(unsigned char key, int x, int y);
void Init();

ParticleSimulation* gSim;
TwoScaleStateSub* gRenderer;
TriangleMesh* gObstacle;
ObstacleRenderer* gObstacleRenderer;
ObstacleGrid* gObstacleGrid;
BoundaryMap* gBoundaryMap;
bool gPause;

int main (int argc, char* argv[]) 
{
    cudaGLSetGLDevice(0);    
    glutInit(&argc, argv);
    glutInitContextVersion(3, 3);
    glutInitContextFlags(GLUT_FORWARD_COMPATIBLE);
    glutInitContextProfile(GLUT_CORE_PROFILE);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
	//int width = glutGet(GLUT_SCREEN_WIDTH);
	//int height = glutGet(GLUT_SCREEN_HEIGHT);
    glutInitWindowSize(1280, 800);
    glutCreateWindow("SPH CUDA");
    //glutFullScreen();
	glewExperimental = TRUE;
	glewInit();
    Init();
	glutDisplayFunc(display);
	glutKeyboardFunc(keyboard);
	glutMainLoop();

    return 0;
}

void Init() 
{
    //try {
        gPause = false;
        
        gSim = ParticleSimulation::Example01();
        gSim->Init();
        gSim->Bind();
        
        gRenderer = new TwoScaleStateSub(*gSim, 1280, 800);
        gRenderer->setCamera(0.0f, 0.0f, 1.3f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f);
	    gRenderer->setPerspective(60.0f, 1280.0f/800.0f, 0.1f, 100.0f);
          /*
        gObstacle = new TriangleMesh("icosphere.obj");
        gObstacle->fit(Rectangle3f(Vector3f(-0.75f, -0.75f, -0.75f), 
            Vector3f(0.75f, 0.75f, 0.75f)));

        BoundaryMapConfiguration config;
        config.compactSupport = 0.017f;
        config.dx = config.compactSupport/3.0f;
        config.restDistance = config.compactSupport*2.0f;

        gBoundaryMap = new BoundaryMap(config);


         gBoundaryMap->AddCanvas(*gObstacle);
         std::cout << "saving boundary map" << std::endl;
         gBoundaryMap->Save("icosphere.txt");
         std::cout << "finished mitm rotz" << std::endl;
         system("pause");
  */

  /*      gObstacleRenderer = new ObstacleRenderer();
        gObstacleRenderer->setCamera(1.0f, 2.0f, 4.6f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f);
	    gObstacleRenderer->setPerspective(60.0f, 1280.0f/800.0f, 0.1f, 100.0f);
        gObstacleRenderer->setObstacle(*gObstacle);
*/

     //   //gSim->SaveInfoTable("tralala.txt");
     //   //gSim->SaveParticleInfo("particle_info.txt");

     //   printf("+++ finished +++\n");

     //   //delete gRenderer;
     //   //delete gSim;

    //} 
    //catch (std::runtime_error e)
    //{
    //    std::cout << e.what();
    //}
}

void display () 
{
    //gObstacleRenderer->draw();

    gRenderer->render();
    
    if (!gPause) 
    {
        gSim->Advance();
    }
    
    //   system("pause");
}

void keyboard (unsigned char key, int x, int y)
{
    switch (key) {
    case 'p':
        if (gPause) {
            gPause = false;
        } else {
            gPause = true;
        }
        break;
    case 'r':
        gSim->SetNPartThresh(1.0f);
        break;
    case 'f':
        gSim->SetNPartThresh(-1.0f);
        break;
    case 't':
        gSim->IncreaseCmDistanceThresh();
        break;
    case 'g':
        gSim->DecreaseCmDistanceThresh();
        break;
    case 27:
        delete gSim;
		exit(0);
		break;
	default:
		return;
	}
}