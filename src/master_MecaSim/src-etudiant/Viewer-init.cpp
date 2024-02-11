/** \file Viewer-init.cpp
 * \brief Fonctions d initialisation des Mesh des objets a afficher.
 */

#include <cassert>
#include <cmath>
#include <cstdio>
#include <iostream>

#include "draw.h"
#include "Viewer.h"

#include "wavefront.h"

using namespace std;


Vector terrainNormal(const Image& im, const int i, const int j)
{
    int ip = i - 1;
    int in = i + 1;
    int jp = j - 1;
    int jn = j + 1;

    Vector a(ip,im(ip,j).r,j);
    Vector b(in,im(in,j).r,j);
    Vector c(i,im(i,jp).r,jp);
    Vector d(i,im(i,jn).r,jn);

    Vector ab = normalize(b-a);
    Vector cd = normalize(d-c);

    Vector n =  -cross(ab,cd);

    return n;
}


void Viewer::init_terrain(Mesh& m_terrain ,const Image& im)
{
	m_terrain = Mesh(GL_TRIANGLE_STRIP);

    for (int i = 1; i < im.width() - 2; i++)
    {
        for (int j = 1; j < im.height() - 1; j++)
        {
			
            //terrainNormal fonction pour faire un affichage correct

            
            m_terrain.normal(terrainNormal(im, i + 1, j));                         // Calcul de la normal en (i+1, j).
            //m_terrain.texcoord(float(i + 1) / im.width(), float(j) / im.height()); // Coordonnées de textures par rapport à l'image.
            m_terrain.vertex(Point(i + 1, 25.f * im(i + 1, j).r, j));       
            m_terrain.texcoord(float(i + 1) / im.width(), float(j) / im.height());   // Sommet à partir du point (i+1, j) de mon image.

            
            m_terrain.normal(terrainNormal(im, i, j));                         // Calcul de la normal en (i, j).
            //m_terrain.texcoord(float(i) / im.width(), float(j) / im.height()); // Coordonnées de textures par rapport à l'image.
            m_terrain.vertex(Point(i, 25.f * im(i, j).r, j));                  // Sommet à partir du point (i, j) de mon image.
            m_terrain.texcoord(float(i) / im.width(), float(j) / im.height());
            
            
            
		}
		m_terrain.restart_strip();
	}


}



/*
 * Creation du maillage des axes pour les visualiser facilement.
 */
void Viewer::init_axe()
{
    m_axe = Mesh(GL_LINES);
    
    // Axe x en rouge
    m_axe.color( Color(1, 0, 0));
    m_axe.vertex( 0,  0, 0);
    m_axe.vertex( 1,  0, 0);
    
    // Axe y en vert
    m_axe.color( Color(0, 1, 0));
    m_axe.vertex( 0,  0, 0);
    m_axe.vertex( 0,  1, 0);
    
    // Axe z en bleu
    m_axe.color( Color( 0, 0, 1));
    m_axe.vertex( 0,  0, 0);
    m_axe.vertex( 0,  0, 1);
    
}


/*
 * Creation du maillage d'une grille.
 */
void Viewer::init_grid()
{
    m_grid = Mesh(GL_LINES);
    
    // Couleur de la grille
    m_grid.color( Color(1, 1, 1));
    
    for(int i=-5;i<=5;++i)
        for(int j=-5;j<=5;++j)
        {
            m_grid.vertex( -5, 0, j);
            m_grid.vertex( 5, 0,  j);
            
            m_grid.vertex( i, 0, -5);
            m_grid.vertex( i, 0, 5);
            
        }
}


/*
 * Creation du maillage du plan de collision.
 */
void Viewer::init_plan(float x, float y, float z)
{
    m_plan = Mesh(GL_LINES);
    
    // Couleur de la grille
    m_plan.color( Color(1, 0, 0));
    
    for(int i=-5;i<=5;++i)
        for(int j=-5;j<=5;++j)
        {
            m_plan.vertex( -5+x, y, j+z);
            m_plan.vertex( 5+x, y,  j+z);
            
            m_plan.vertex( i+x, y, -5+z);
            m_plan.vertex( i+x, y, 5+z);
            
        }
}


/*
 * Creation du maillage d un cube.
 */
void Viewer::init_cube()
{
    // Definition des 8 sommets du cube
    //                          0           1           2       3           4           5       6           7
    // 8 sommets - chaque sommet a coodonnees en x, y, z
    static float pt[8][3] = { {-1,-1,-1}, {1,-1,-1}, {1,-1,1}, {-1,-1,1}, {-1,1,-1}, {1,1,-1}, {1,1,1}, {-1,1,1} };
   
    // 6 faces - chacune des faces a 4 sommets
    static int f[6][4] = {    {0,1,2,3}, {5,4,7,6}, {2,1,5,6}, {0,3,7,4}, {3,2,6,7}, {1,0,4,5} };
    
    // 6 normales (une normale par face) - chaque normale a coodonnees en x, y, z
    static float n[6][3] = { {0,-1,0}, {0,1,0}, {1,0,0}, {-1,0,0}, {0,0,1}, {0,0,-1} };

    
    // Maillage de type GL_TRIANGLE_STRIP
    m_cube = Mesh(GL_TRIANGLE_STRIP);
    
    // Couleur du cube
    m_cube.color( Color(1, 1, 1) );
    
    // Boucle sur les 6 faces du cube
    for (int i=0;i<6;i++)
    {
        // Definition de la normale a la face i
        m_cube.normal(  n[i][0], n[i][1], n[i][2] );
        
        // Definition du sommet 0 de la face i : coordonnees texture + coordonnees geometrique
        m_cube.texcoord( 0,0 );
        m_cube.vertex( pt[ f[i][0] ][0], pt[ f[i][0] ][1], pt[ f[i][0] ][2] );
        
        // Definition du sommet 1 de la face i  : coordonnees texture + coordonnees geometrique
        m_cube.texcoord( 1,0);
        m_cube.vertex( pt[ f[i][1] ][0], pt[ f[i][1] ][1], pt[ f[i][1] ][2] );
        
        // Definition du sommet 3 de la face i  : coordonnees texture + coordonnees geometrique
        m_cube.texcoord(0,1);
        m_cube.vertex(pt[ f[i][3] ][0], pt[ f[i][3] ][1], pt[ f[i][3] ][2] );
        
        // Definition du sommet 2 de la face i : coordonnees texture + coordonnees geometrique
        m_cube.texcoord(1,1);
        m_cube.vertex( pt[ f[i][2] ][0], pt[ f[i][2] ][1], pt[ f[i][2] ][2] );
        
        m_cube.restart_strip();
        
    }
}


void Viewer::init_cubeMap()
{
    static float pt[8][3] = { {-20, -20, -20}, {20, -20, -20}, {20, -20, 20}, {-20, -20, 20}, {-20, 20, -20}, {20, 20, -20}, {20, 20, 20}, {-20, 20, 20} };
    static int f[6][4] = { {3, 2, 1, 0}, {6, 7, 4, 5}, {6, 5, 1, 2}, {4, 7, 3, 0}, {7, 6, 2, 3}, {5, 4, 0, 1} };

    static float text[6][4][2] =
    {
        // Face {3,2,1,0}.
        {{1 / 4.f, 0}, {1 / 2.f, 0}, {1 / 2.f, 1 / 3.f}, {1 / 4.f, 1 / 3.f}},
        // Face {6,7,4,5}.
        {{1 / 2.f, 1}, {1 / 4.f, 1}, {1 / 4.f, 2 / 3.f}, {1 / 2.f, 2 / 3.f}},
        // Face {6,5,1,2}.
        {{3 / 4.f, 2 / 3.f}, {1 / 2.f, 2 / 3.f}, {1 / 2.f, 1 / 3.f}, {3 / 4.f, 1 / 3.f}},
        // Face {4,7,3,0}.
        {{1 / 4.f, 2 / 3.f}, {0, 2 / 3.f}, {0, 1 / 3.f}, {1 / 4.f, 1 / 3.f}},
        // Face {7,6,2,3}.
        {{1, 2 / 3.f}, {3 / 4.f, 2 / 3.f}, {3 / 4.f, 1 / 3.f}, {1, 1 / 3.f}},
        // Face {5,4,0,1}.
        {{1 / 2.f, 2 / 3.f}, {1 / 4.f, 2 / 3.f}, {1 / 4.f, 1 / 3.f}, {1 / 2.f, 1 / 3.f}},
    };

    static float n[6][3] = { {0, 1, 0}, {0, -1, 0}, {-1, 0, 0}, {1, 0, 0}, {0, 0, -1}, {0, 0, 1} };

    m_cubeMap = Mesh(GL_TRIANGLE_STRIP);

    for (int i = 0; i < 6; i++)
    {
        m_cubeMap.normal(n[i][0], n[i][1], n[i][2]);

        m_cubeMap.texcoord(text[i][0][0], text[i][0][1]);
        m_cubeMap.vertex(pt[f[i][0]][0], pt[f[i][0]][1], pt[f[i][0]][2]);

        m_cubeMap.texcoord(text[i][3][0], text[i][3][1]);
        m_cubeMap.vertex(pt[f[i][3]][0], pt[f[i][3]][1], pt[f[i][3]][2]);

        m_cubeMap.texcoord(text[i][1][0], text[i][1][1]);
        m_cubeMap.vertex(pt[f[i][1]][0], pt[f[i][1]][1], pt[f[i][1]][2]);

        m_cubeMap.texcoord(text[i][2][0], text[i][2][1]);
        m_cubeMap.vertex(pt[f[i][2]][0], pt[f[i][2]][1], pt[f[i][2]][2]);

        m_cubeMap.restart_strip();
    }
}







/*
 * Creation du maillage d'une sphere - centre (0, 0, 0) - rayon = r.
 */
void Viewer::init_sphere()
{
    const int divBeta = 26;
    const int divAlpha = divBeta/2;
    int i,j;
    float beta, alpha, alpha2;
    float rayon = 0.05;
    
    m_sphere = Mesh(GL_TRIANGLE_STRIP);

    //m_sphere.color( Color(0, 0, 1) );
        
    for(i=0;i<divAlpha;++i)
    {
        alpha = -0.5f*M_PI + float(i)*M_PI/divAlpha;
        alpha2 = -0.5f*M_PI + float(i+1)*M_PI/divAlpha;
        
        for(j=0;j<divBeta;++j)
        {
            beta = float(j)*2.f*M_PI/(divBeta-1);
            
            m_sphere.texcoord(beta/(2.0f*M_PI), 0.5f+alpha2/M_PI);
            
            m_sphere.normal( - Vector(rayon * cos(alpha2)*cos(beta),
                                    rayon * sin(alpha2),
                                    rayon * cos(alpha2)*sin(beta)) );
            
            m_sphere.vertex( Point(rayon * cos(alpha2)*cos(beta),
                                   rayon*sin(alpha2),
                                   rayon * cos(alpha2)*sin(beta))   );
            
            
            m_sphere.texcoord(beta/(2.0f*M_PI), 0.5f+alpha/M_PI);
            
            m_sphere.normal( - Vector(rayon * cos(alpha)*cos(beta),
                                    rayon * sin(alpha),
                                    rayon * cos(alpha)*sin(beta)) );
            
            m_sphere.vertex( Point(rayon * cos(alpha)*cos(beta),
                                   rayon * sin(alpha),
                                   rayon * cos(alpha)*sin(beta)) );
        }

        
        m_sphere.restart_strip();
    }
}

void Viewer::init_vega()
{
    m_vega = read_mesh("data//Vega/vega.obj");
}





/*
 * Fonction d initialisation.
 */
int Viewer::init()
{
    // Etats par defaut openGL
    glClearColor(0.5f, 0.6f, 0.8f, 1);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glClearDepthf(1);
    glDepthFunc(GL_LESS);
    glEnable(GL_DEPTH_TEST);
    glFrontFace(GL_CW);
    glCullFace(GL_BACK);
    
    if (mb_cullface)
        glEnable(GL_CULL_FACE);
    else
        glDisable(GL_CULL_FACE);        // good for debug
    glEnable(GL_TEXTURE_2D);
    
    // Camera
    m_camera.lookat( Point(0,0,0), 30);
    
    //Orbiter() : m_center(), m_position(), m_rotation(), m_size(5.f), m_radius(5.f) {}
    // m_camera.lookat( Point(0,0,0), 150 );
    
    // Lumiere
    //gl.light( Point(10, -200, 3),Blue( ));
    gl.light(Point(0, -10, 0), White());    // Deuxième source lumineuse
    //gl.light(Point(-50, 50, 0), Blue());


    

    m_terrainAlti = read_image("data/terrain/terrain.png");
    m_terrain_texture = read_texture(0, "data/terrain/terrain_texture.png");
    m_sphere_texture = read_texture(0, "data/monde.jpg");
    earth_texture = read_texture(0, "data/monde.jpg");
    m_cubeMap_texture = read_texture(0, "data/cubemap/space4.png");
    planet_texture = read_texture(0, "data/planet/planet.png"); 
    mytunisia = read_texture(0, "data/planet/Tunisia.png");


    
    // Chargement des textures utilisees dans la scene
    // Texture pour le cube
   //  m_cube_texture = read_texture(0, "../data/debug2x2red.png");
    m_cube_texture = read_texture(0, smart_path("data/debug2x2red.png"));

    // Texture pour le tissu
    //m_tissu_texture = read_texture(0, "data/papillon.png");
    //m_tissu_texture = read_texture(0, "data/textures/tissu1.png");
   // m_tissu_texture = read_texture(0, "data/textures/tissu2.jpg");
    m_tissu_texture = read_texture(0, smart_path("data/textures/dollars.jpg"));


    // Appel des procedures d initialisation des objets de la scene
    // Pour les objets non simules
    // Exemple : init_votreObjet();
    init_axe();
    init_grid();
    init_cube();
    init_sphere();
    init_terrain(m_terrain,m_terrainAlti);
    init_cubeMap();
    init_vega();
    //glGetError();
    
    // Creation du plan (x, y, z) - plan utilise pour les ObjetSimule::Collision(x, y, z);
    // Rq : pas vraiment le plan, mais < x, < y, < z
    init_plan(0, 0, 0);
    
    // Initialisation du Tps
    Tps = 0;
    
    // Point interaction
    MousePos = Vector(0, 0, 0);
    
    return 0;
    
}


