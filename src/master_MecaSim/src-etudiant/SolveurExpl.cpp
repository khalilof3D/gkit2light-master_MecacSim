/*
 * SolveurExpl.cpp : Application schemas semi-implicite sur les objets.
 * Copyright (C) 2016 Florence Zara, LIRIS
 *               florence.zara@liris.univ-lyon1.fr
 *               http://liris.cnrs.fr/florence.zara/
 *
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 2.1 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

/** \file Calculs.cpp
 Fonctions de calculs communes aux objets simules.
 \brief Fonctions de calculs communes aux objets simules.
 */

#include <stdio.h>
#include <math.h>
#include <vector>
#include <iostream>

#include "vec.h"
#include "ObjetSimule.h"
#include "Viewer.h"
#include "SolveurExpl.h"


using namespace std;



/**
 * Calcul de l acceleration des particules
 * avec ajout de la gravite aux forces des particules
 * et ajout de la force due au vent sur une des particules du maillage
 * et reinitialisation des forces.
 */
void SolveurExpl::CalculAccel_ForceGravite(Vector g,
    int nb_som,
    std::vector<Vector>& A,
    std::vector<Vector>& Force,
    std::vector<float>& M)
{
 


 			//// Cas SPH
            // On a calcule dans Force[i] : fij / rho_i
            // Il ne reste qu a ajoute le vecteur g
            // a_i = fij / rho_i + g

        //float side = sqrt(nb_som);
        for (auto i=0; i<nb_som; i++){


            // div et mod 
            // int div = i / side;
            //     
            


                
                if(M[i]> 0){
                    A[i] = Force[i] + M[i] * g;
                    Force[i] = Vector(0,0,0);
                }
        }

            

            
    
}//void


/*! Calcul des vitesses et positions : 
 *  Formule d Euler semi-implicite :
 *  x'(t+dt) = x'(t) + dt x"(t)
 *  x(t+dt) = x(t) + dt x'(t+dt)
 */
void SolveurExpl::Solve(float visco,
    int nb_som,
    int Tps,
    std::vector<Vector>& A,
    std::vector<Vector>& V,
    std::vector<Vector>& P)

    
{
    
    
    float c1 = 0.8f;
    float c2 = 0.06f;
    Vector center = Vector(0,0,0);

    float max = 0.0f;
    
    
    for (auto i = 0; i < nb_som; i++) {
        center = center + P[i] ;
    }
    center = center/nb_som;

    for (auto i = 0; i < nb_som; i++) {
		float distance = length2(center - P[i]);
        if (distance > max) {
			max = distance;
		}
	}


    
    
    for (auto i=0; i<nb_som; i++){
        

        //force de déchirure
        
        float distance2center = -600*exp(-5*length2(center - P[i]));

        //float distance2center = 0.5*(max - length2(center - P[i]));



        //Vector normals = ObjetSimuleMSS->_vertex_normals[i];



        //portance et résistance de l'air

        // Résistance de l'air, sens inverse de la vitesse
        Vector air = -c1 * V[i] * length(V[i]);

        // Portance, sens vertical
        Vector portance = (c2 * length2(V[i]) * distance2center) * Vector(0, 1, 0);

        A[i] = A[i] + air + portance;




// Bonnes valeurs pour le cas horizontal

float c1 = 0.4f;
float c2 = 0.2f;
Vector n = Vector(0, 1, 0) * length2(center - P[i]) / nb_som;

//A[i] = A[i] + c2 * length(V[i]) * length(V[i]) * Vector(0, 1, 0) - 500 * n - c1 * V[i] * length(V[i]);
        

        
        V[i] = visco * (V[i] + _delta_t*A[i]);
        P[i] = P[i] + _delta_t*V[i];
    }



       
}//void

