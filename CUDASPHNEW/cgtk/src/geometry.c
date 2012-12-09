#include "../include/geometry.h"
#include "../include/error.h"
#ifdef _WIN32
    typedef int bool;
    #define false 0
    #define true 1
#elif
    #include <stdbool.h>
#endif
#include <memory.h>
#include <stdio.h>
#include <stdlib.h>

bool read_obj_file(const char* filename, CGTKObjFile* objFile);

CGTKObjFile* cgtkObjFileAlloc(const char* filename) {
	CGTKObjFile* objFile;

	if (filename == NULL) {
	    return NULL;
	}

	objFile = (CGTKObjFile*) malloc(sizeof(CGTKObjFile));

	if (objFile == NULL) {
	    return NULL;
	}

	memset(objFile, 0, sizeof(CGTKObjFile));

	/* read file and save to object, clean up and return if the file is not a
    ** valid .obj file.
	*/
	if (!read_obj_file(filename, objFile)) {
		CGTK_DUMP_ERROR("Could not read obj file.");
	    free(objFile);
		return NULL;
	}

	return objFile;
}

void cgtkObjFileFree(CGTKObjFile** objFile) {
	if (*objFile == NULL) {
	    return;
	}
	
	free(*objFile);

	*objFile = NULL;
}

/*
** Aux. enum for parsing the object file. It represents the states of a state
** machine, that checks the obj file. 
*/
typedef enum {
	START = 0,
    END,
	PARSE_VERTEX,
    PARSE_FACE,
	IDLE,
	ERROR
} State;

/*
** Aux. function for parsing the object file. It represents the transitions of 
** a state machine, that checks the obj file. 
*/
State delta(State s, char c) {
    switch (s) {

        /* transitions when current state is: START */
        case START:
            switch (c) {
                case 'v':
                    return PARSE_VERTEX;
                case 'f':
                    return PARSE_FACE;
                case '\n':
                    return START;
                case EOF:
                    return END;
                default:
                    return IDLE;
            }

        /* transitions when current state is: PARSE_VERTEX */
        case PARSE_VERTEX:
            switch (c) {
                default:
                    return IDLE;
            }


        /* transitions when current state is: PARSE_VERTEX */
        case PARSE_FACE:
            switch (c) {
                default:
                    return IDLE;
            }

        /* transitions when current state is: IDLE */
        case IDLE:
            switch (c) {
                case '\n':
                    return START;
                case EOF:
                    return END;
                default:
                    return IDLE;
            }
        
    }   

}

/*
** Reads contents of the obj file and stores them in the struct.
*/
bool read_obj_file(const char* filename, CGTKObjFile* objFile) {
	FILE* pFile = NULL;
	char c;
	State state = START;
	float v[3];
	unsigned int vIdx = 0, fIdx = 0;
		
	/* open file */
	pFile = fopen(filename, "r");
	
	if (pFile == NULL) {
		fclose(pFile);
		return false;
	}
	
	/* parse file: 1st run to get number of vertices, faces etc. */	
	while (state != END) {
		c = fgetc(pFile);
        
		switch (state) {
			case START:
                state = delta(state, c);
				break;
			case IDLE:
                state = delta(state, c);
				break;
			case PARSE_VERTEX:
				objFile->nVertices++;
			    state = delta(state, c);
				break;		
            case PARSE_FACE:
                objFile->nFaces++;
                state = delta(state, c);
                break;
            case ERROR:
				CGTK_DUMP_ERROR("Error while reading file.");
                return false;
		}
		
	}

	
	/* allocate memory for geometry data */
	objFile->vertices = (float*) malloc(3*objFile->nVertices*sizeof(float));
	
	if (objFile->vertices == NULL) {
		fclose(pFile);
		return false;
	}
	
	objFile->indices = (unsigned int*) malloc(3*objFile->nFaces*
		sizeof(unsigned int));
		
	if (objFile->indices == NULL) {
		free(objFile->vertices);
		fclose(pFile);
		return false;
	}
	
	/* parse file: 2nd run to get fill geometry data */	
	rewind(pFile);
	state = START;
	while (state != END) {
		c = fgetc(pFile);
        
		switch (state) {
			case START:
                state = delta(state, c);
				break;
			case IDLE:
                state = delta(state, c);
				break;
			case PARSE_VERTEX:
				fscanf(pFile, " %f %f %f", &objFile->vertices[vIdx],
					&objFile->vertices[vIdx + 1],
					&objFile->vertices[vIdx + 2]);
				vIdx += 3;
			    state = delta(state, c);
				break;		
            case PARSE_FACE:
				fscanf(pFile, " %d %d %d", &objFile->indices[fIdx],
					&objFile->indices[fIdx + 1],
					&objFile->indices[fIdx + 2]);

				/* adjust indices to c array style, i.e. substract 1 */	
				objFile->indices[fIdx]--;
				objFile->indices[fIdx + 1]--;
				objFile->indices[fIdx + 2]--;
				
				fIdx += 3;
				state = delta(state, c);
                break;
            case ERROR:
				CGTK_DUMP_ERROR("Error while reading file.");
                return false;
		}
		
	}
	fclose(pFile);
	
	
	
	return true;	
}


void cgtkObjFileDump(CGTKObjFile* objFile) {
	unsigned int i;

	if (objFile == NULL) {
	    return;
	}

	printf("Dumping vertices ...\nTotal # of vertices: %d\n", 
		objFile->nVertices);

	for (i = 0; i < objFile->nVertices; i++) {
	    printf("[%.3f %.3f %.3f]\n", objFile->vertices[3*i], 
			objFile->vertices[3*i + 1], objFile->vertices[3*i + 2]);
	}

	printf("Dumping face indices ...\nTotal # of faces: %d\n", 
		objFile->nFaces);

	for (i = 0; i < objFile->nFaces; i++) {
	    printf("[%d %d %d]\n", objFile->indices[3*i], 
			objFile->indices[3*i + 1], objFile->indices[3*i + 2]);
	}
}

