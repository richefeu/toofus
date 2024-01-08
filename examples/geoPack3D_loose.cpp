#include <iostream>

#include "../geoPack3D.hpp"

int main(int argc, char const *argv[]) {

  GeoPack3D GP(0.3,  // rmin
               1.0,  // rmax
               5000, // k (number of attemps)
               0.0,  // xmin
               10.0, // xmax
               0.0,  // ymin
               10.0, // ymax
               0.0,  // zmin
               10.0, // zmax
               0.0,  // gapTol, a (positive) tolerance on the space between the spheres
               10000 // max number of placed spheres
  );

  // pour avoir un nouveau tirage à chaque fois
  GP.seedTime();

  // grace à ce truc, le z local sera limité à une certaine valeur
  // et on s'attend ici que la densité globale évolue comme z (je sais que ce n'est obligatoire)

  GP.limit_localNumberNeighbors(0.1, // dist max for z computation (mettre une petite valeur comme Rmin/10)
                                6    // maximum coordination
  );

  GP.limit_localSolidFraction(2.0, // radius into whitch the local solid fraction is computed
                              0.2  // max local solid fraction
  );

  GP.distMin = 0.0; // on peut mettre une autre valeur...

  // c'est là qu'on pack en périodique (xmin, xmax...)
  GP.execPeriodic();
  std::cout << "Cell volume: " << (GP.xmax - GP.xmin) * (GP.ymax - GP.ymin) * (GP.zmax - GP.zmin) << '\n';
  std::cout << "Overall solid fraction (by assuming there is no overlap): " << GP.getSolidFraction() << '\n';

  // on sauvegarde la boite (c'est util pour la visualisation)
  GP.saveBox("box.txt");
  GP.save("packing.txt");
  // utiliser seeSphere pour visualiser
  // touche i pour afficher les particules ghost 
  // touche l pour afficher la cellule périodique
  // touche k pour afficher les connections entre particules

  return 0;
}