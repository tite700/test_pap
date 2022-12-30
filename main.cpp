#include "AMere.h"
#include "EDPcomp.h"
#include "EDPred.h"
#include "Graph.h"

#include <stdio.h>
#include <ostream>
#include <iostream>

int main(){

    /*----------------*/
    /*  TEST ON CranckPut */
    /*-----------------*/
    CranckPut cranck1=CranckPut(1.0,300.0,1000,1000,1);
    cranck1.Methode(100.0,0.1,0.1);
    for(int i=0;i<40;i++){std::cout<<cranck1(1,i)<<std::endl;;}

    /*----------------*/
    /*  TEST ON CranckCall */
    /*-----------------*/

    CranckCall cranck2=CranckCall(1.0,300.0,1000,1000,1);
    cranck2.Methode(100.0,0.1,0.1);
    for(int i=0;i<40;i++){std::cout<<cranck2(1,i)<<std::endl;;}

    /*----------------*/
    /*  TEST ON DiffImpPut */
    /*-----------------*/

    DiffImpPut diffimp1=DiffImpPut(1.0,300.0,1000,1000,1);
    diffimp1.Methode(100.0,0.1,0.1);
    for(int i=0;i<40;i++){std::cout<<diffimp1(1,i)<<std::endl;;}

    /*----------*/
    /* ESSAI D'AFFICHAGE */
    /*----------*/

    // Création de l'objet graphique
    Graph graph1(640, 480, "cranck1");
    Graph graph2(640, 480, "cranck2");
    Graph graph3(640, 480, "diffimp1");

    // Ajout des valeurs au graphique
    for (int i=0;i<40;i++){
        graph1.addPoint(i, cranck1(1,i));
        graph2.addPoint(i, cranck2(1,i));
        graph3.addPoint(i, diffimp1(1,i));

    }

    // Affichage du graphique
    graph1.draw();
    graph2.draw();
    graph3.draw();

    // Attente de fermerture de la fenêtre
    bool running = true;
    while (running) {
        SDL_Event event;
        while (SDL_PollEvent(&event)) {
            if (event.type == SDL_QUIT) {
                running = false;
            }
        }
    }

    return 0;
}
