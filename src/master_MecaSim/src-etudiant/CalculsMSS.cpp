/*
 * CalculsMSS.cpp :
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

/** \file CalculsMSS.cpp
Programme calculant pour chaque particule i d un MSS son etat au pas de temps suivant 
 (methode d 'Euler semi-implicite) : principales fonctions de calculs.
\brief Fonctions de calculs de la methode semi-implicite sur un systeme masses-ressorts.
*/ 

#include <stdio.h>
#include <math.h>
#include <vector>
#include <iostream>

#include "vec.h"
#include "ObjetSimule.h"
#include "ObjetSimuleMSS.h"
#include "Viewer.h"

#include "mesh.h"

using namespace std;





/**
* Calcul des forces appliquees sur les particules du systeme masses-ressorts.
 */
void ObjetSimuleMSS::CalculForceSpring()
{
	// Calcul des forces élastiques appliquées aux particules du système masse-ressort
	for (int i = 0; i < _SystemeMasseRessort->GetNbParticule(); i++)
	{
		Particule* particuleCourante = _SystemeMasseRessort->GetPartList()[i];
		Vector forceTotaleSurParticule = Vector(0.0, 0.0, 0.0); // Initialisation de la force totale sur la particule courante

		for (int j = 0; j < particuleCourante->GetNbVoisins(); j++)
		{
			Ressort* ressortConnecte = particuleCourante->GetRessortList()[j];
			Particule* particuleA = ressortConnecte->GetParticuleA();
			Particule* particuleB = ressortConnecte->GetParticuleB();

			// Si la particule courante est la particule B du ressort, on inverse pour standardiser le calcul
			if (particuleB == _SystemeMasseRessort->GetParticule(i))
			{
				swap(particuleA, particuleB);
			}

			float raideurRessort = ressortConnecte->GetSpring()->_Raideur;
			float longueurActuelle = distance(Point(P[particuleA->_Id]), Point(P[particuleB->_Id]));
			float longueurRepos = ressortConnecte->GetSpring()->_L0;
			Vector directionRessort = normalize(P[particuleB->_Id] - P[particuleA->_Id]);

			// Calcul de la force élastique
			Vector forceElastique = raideurRessort * (longueurActuelle - longueurRepos) * directionRessort;

			// Calcul de l'amortissement
			Vector vecteurAB = particuleB->GetPosition() - particuleA->GetPosition();
			Vector vitesseRelative = vecteurAB / longueurActuelle;
			float coefficientAmortissement = dot(-ressortConnecte->GetAmortissement(), vitesseRelative);
			Vector forceAmortissement = vecteurAB / length(vecteurAB) * coefficientAmortissement;

			// Addition des forces élastique et d'amortissement
			forceTotaleSurParticule = forceTotaleSurParticule + forceElastique + forceAmortissement;

			// Gestion de la déchirure du ressort si la longueur dépasse un seuil spécifique
			if (longueurActuelle > 1) // Le seuil de déchirure est arbitrairement fixé à 1 pour cet exemple
			{
				particuleCourante->GetRessortList().erase(particuleCourante->GetRessortList().begin() + j);
				// Note : Faire attention à la modification de la liste pendant l'itération
			}
		}

		// Mise à jour de la force totale exercée sur la particule courante
		Force[particuleCourante->_Id] = forceTotaleSurParticule;
	}
			
}//void
 //pt fixe 0 70

/**
 * Gestion des collisions avec le sol.
 */
void ObjetSimuleMSS::Collision()
{
	/// Arret de la vitesse quand touche le plan
	for (auto i = 0; i < P.size(); i++) {
		
		
		//9.99 si il y a une collision déchirante avec la sphere pour des raisons esthetiques
		if (P[i].y < -12.99) {
			P[i].y = -12.99;
			V[i] = V[i] / 10;
		}
	}
	//CollisionMesh();
}

//Viewer m_viewer;

/**
 * Gestion des collisions avec le mesh. Non utilisée.
 */
//void ObjetSimuleMSS::CollisionMesh()
	
//    float thr = 2.0;




	//for (int j = 0; j < P.size(); j++)
	
		//for (int i = 0; i < terrain.triangle_count(); i++)
		//for (int i = 0; i < 10; i++)
		
			//TriangleData t = terrain.triangle(i);
			//TriangleData t = TriangleData();

			
			//t = Vector(5, 4, 3);
			//Vector ta = Vector(t.a.x, t.a.y, t.a.z);

			//if (length(P[j] - ta) > thr)
				//continue;

			//Vector tb = Vector(t.b.x, t.b.y, t.b.z);
			//Vector tc = Vector(t.c.x, t.c.y, t.c.z);
			//Vector n = normalize(cross(tb - ta, tc - ta));

			//float dist = dot(n, P[j] - ta);

			// Si le point est sous le plan du triangle (et proche du triangle)
			//if (dist > 0.01)
			
				// Calculer la projection du point sur le plan du triangle
				//Vector proj = P[j] - n * dist;

				// Vérifier si le point est à l'intérieur du triangle
				//Vector v0 = ta - proj;
				//Vector v1 = tb - proj;
				//Vector v2 = tc - proj;

				//float dot00 = dot(v0, v0);
				//float dot01 = dot(v0, v1);
				//float dot02 = dot(v0, v2);
				//float dot11 = dot(v1, v1);
				//float dot12 = dot(v1, v2);

				//float invDenom = 1 / (dot00 * dot11 - dot01 * dot01);
				//float u = (dot11 * dot02 - dot01 * dot12) * invDenom;
				//float v = (dot00 * dot12 - dot01 * dot02) * invDenom;

				// Si le point est à l'intérieur du triangle
				//if ((u >= 0) && (v >= 0) && (u + v < 1))
				
//					P[j] = proj;




/**
 * Gestion des collisions avec la sphere.
 */
void ObjetSimuleMSS::CollisionSphere(Point p, double r) 
{
	


	float px = p.x;
	float py = p.y;
	float pz = p.z;
	float r2 = r * r;
	
	
	for (int i = 0; i < P.size(); i++) {
		
		
		float x = P[i].x - px;
		float y = P[i].y - py;
		float z = P[i].z - pz;
		float d2 = x * x + y * y + z * z;
		if (d2 < r2) {
			float d = sqrt(d2);
			float l = r / d;
			P[i].x = px + x * l;
			P[i].y = py + y * l;
			P[i].z = pz + z * l;

			//La sphère effectue des transformations une fois atteinte par le tissu et le déchire.
			P[i] = (Vector)(RotationY(1))((Point)P[i]);

			V[i] = V[i] /5 ;


			
		}


		
	}

}// void

