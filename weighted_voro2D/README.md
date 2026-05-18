# laguerre.hpp

Bibliothèque *single-header* C++17 pour le calcul de **diagrammes de Laguerre** (diagrammes de Voronoï pondérés) en 2D.

Aucune dépendance externe. Il suffit de copier `laguerre.hpp` dans votre projet.

---

## Table des matières

1. [Concept](#concept)
2. [Démarrage rapide](#démarrage-rapide)
3. [LaguerreBuilder](#laguerrebuilder)
4. [DCEL](#dcel)
5. [EdgeInfo et for_each_edge](#edgeinfo-et-for_each_edge)
6. [Export SVG](#export-svg)
7. [Performances](#performances)
8. [Internals](#internals)

---

## Concept

Un **diagramme de Laguerre** généralise le diagramme de Voronoï en associant à chaque site un *poids* (rayon `r`). La frontière entre deux sites `i` et `j` n'est plus la médiatrice de leurs centres, mais la **bissectrice de puissance** :

```
2(xj - xi)·x + 2(yj - yi)·y = xj² + yj² - rj² - xi² - yi² + ri²
```

Quand tous les rayons sont égaux, on retrouve le diagramme de Voronoï classique.

---

## Démarrage rapide

```cpp
#include "laguerre.hpp"

int main() {
    lv::LaguerreBuilder<double> b;

    b.add_site(0.25, 0.25, 0.10);
    b.add_site(0.75, 0.25, 0.05);
    b.add_site(0.50, 0.75, 0.15);

    b.set_bbox(0.0, 0.0, 1.0, 1.0); // optionnel
    b.compute_fast2();               // calcul des cellules

    // Export SVG direct
    lv::SVGOptions opt;
    lv::write_svg(b, "out.svg", opt);

    return 0;
}
```

Compilation :

```bash
g++ -std=c++17 -O2 main.cpp -o main
```

---

## LaguerreBuilder

```cpp
lv::LaguerreBuilder<double> b;
```

Le type `T` peut être `float` ou `double`.

### Ajouter des sites

```cpp
b.add_site(x, y, r);

// ou directement via le vecteur public :
b.sites.push_back({x, y, r});
```

Chaque site est un `struct Site { T x, y, r; }`.

### Bounding box

Par défaut, la bbox est calculée automatiquement pour englober tous les disques (`fit_bbox()`). Vous pouvez la fixer manuellement :

```cpp
b.set_bbox(xmin, ymin, xmax, ymax);
```

Les cellules sont clippées sur cette bbox.

### Calcul des cellules

Trois variantes disponibles selon vos besoins :

| Méthode | Complexité théorique | Usage recommandé |
|---|---|---|
| `compute()` | O(n²) | n < 500, ou référence |
| `compute_fast()` | O(n² log n) pire cas, O(n log n) pratique | n < 2 000 |
| `compute_fast2()` | O(n log n) pratique (grille spatiale) | n ≥ 2 000 |

```cpp
b.compute();        // référence, simple
b.compute_fast();   // tri par distance + pruning
b.compute_fast2();  // grille spatiale + pruning par anneau
```

> **Note** : les trois méthodes produisent des résultats identiques. Seule la performance diffère.

### Résultat

Après `compute*()`, les cellules sont disponibles dans :

```cpp
b.cells; // std::vector<std::vector<Vec2<T>>>
         // cells[i] = polygone convexe de la cellule i
```

L'ordre de `cells` correspond à l'ordre de `sites`.

---

## DCEL

La **DCEL** (*Doubly Connected Edge List*) est une structure topologique qui représente le diagramme comme un graphe planaire : sommets, demi-arêtes et faces sont tous reliés entre eux.

### Construction

```cpp
lv::DCEL<double> dcel;
lv::build_true_laguerre_dcel<double>(b.cells, dcel);
```

### Structure interne

```
dcel.vertices   — liste des sommets uniques (Vec2<T>)
dcel.halfedges  — liste des demi-arêtes (stockées par paires)
dcel.faces      — liste des faces (une par cellule)
```

**Vertex**

```cpp
dcel.vertices[i].p  // Vec2<T> : position du sommet
```

**HalfEdge**

```cpp
dcel.halfedges[i].origin  // indice du sommet de départ
dcel.halfedges[i].twin    // indice du demi-arête opposé
dcel.halfedges[i].face    // indice de la face à gauche (-1 si bord)
dcel.halfedges[i].next    // prochain demi-arête dans la face
```

Les demi-arêtes sont stockées par paires : `halfedges[2k]` et `halfedges[2k+1]` sont toujours jumeaux.

**Face**

```cpp
dcel.faces[i].site  // indice du site générateur
dcel.faces[i].edge  // un demi-arête quelconque de cette face
```

### Parcourir le contour d'une face

```cpp
int start = dcel.faces[i].edge;
int e     = start;
do {
    const auto &v = dcel.vertices[dcel.halfedges[e].origin].p;
    // utiliser v.x, v.y
    e = dcel.halfedges[e].next;
} while (e != start);
```

### Arêtes de bord

Un demi-arête est **sur le bord de la bbox** si son jumeau n'appartient à aucune face :

```cpp
if (dcel.halfedges[e].face < 0) {
    // bord extérieur
}
```

---

## EdgeInfo et for_each_edge

`for_each_edge` itère sur chaque arête géométrique **une seule fois** (pas de doublon demi-arête / jumeau) et fournit toutes les informations utiles via un `EdgeInfo`.

```cpp
lv::for_each_edge<T>(dcel, [&](const lv::EdgeInfo<T> &info) {

    info.A          // Vec2<T> : premier sommet
    info.B          // Vec2<T> : second sommet
    info.mid        // Vec2<T> : milieu de l'arête
    info.length     // T       : longueur

    info.he         // int     : indice du demi-arête canonique
    info.site0      // int     : site à gauche de he  (-1 si bord)
    info.site1      // int     : site à gauche du twin (-1 si bord)
    info.border     // bool    : true si arête de bord

    // utilitaire : retourne (min(site0,site1), max(site0,site1))
    auto [smin, smax] = info.sorted_sites();
});
```

`T` est le type utilisé pour le `DCEL` (`float` ou `double`). En pratique avec `DCEL<double>` :

```cpp
double total = 0.0;
lv::for_each_edge<double>(dcel, [&](const lv::EdgeInfo<double> &info) {
    if (!info.border) total += info.length;
});
```

---

## Export SVG

### `write_svg` — rendu rapide des cellules

```cpp
lv::SVGOptions opt;
opt.width       = 800;   // largeur en pixels
opt.height      = 800;   // hauteur en pixels
opt.margin      = 30;    // marge en pixels

opt.draw_cells  = true;  // polygones colorés
opt.draw_sites  = true;  // points des sites
opt.draw_labels = true;  // indices des sites
opt.draw_radii  = true;  // cercles de poids (tiretés)
opt.draw_limits = true;  // bounding box (tiretée)

lv::write_svg(b, "out.svg", opt);
```

### `write_svg_dcel_debug` — rendu avec topologie DCEL

Affiche les arêtes issues de la DCEL avec leur classification (bord / intérieure) et le label des sites adjacents.

```cpp
lv::write_svg_dcel_debug(dcel, b, "debug.svg", opt);
```

Les arêtes intérieures apparaissent en bleu, les arêtes de bord en rose.

---

## Performances

Mesures indicatives sur distribution uniforme aléatoire, compilé avec `-O2` :

| n sites | `compute()` | `compute_fast()` | `compute_fast2()` |
|---:|---:|---:|---:|
| 500 | ~29 ms | ~15 ms | ~2 ms |
| 2 000 | ~259 ms | ~157 ms | ~9 ms |
| 5 000 | ~1 493 ms | ~926 ms | ~19 ms |

> Les gains de `compute_fast2()` sont maximaux sur les distributions uniformes. En cas de rayons très hétérogènes (certaines cellules très grandes), la grille spatiale reste correcte mais le pruning est moins efficace.

### Choisir la bonne méthode

```
n < 300   → compute()        (plus simple à déboguer)
n < 2000  → compute_fast()   (bon compromis)
n ≥ 2000  → compute_fast2()  (grille spatiale, meilleur cas pratique)
```

---

## Internals

Détails des choix d'implémentation internes, utiles si vous souhaitez modifier le header.

### DCEL : tables de hachage

`vmap` (vertices) et `emap` (halfedges) utilisent `std::unordered_map` avec des hashers dédiés (`DCELKeyHash`, `PairHash`) pour des insertions et recherches en O(1) amorti, contre O(log n) avec `std::map`.

### for_each_edge : pas de déduplication explicite

Les demi-arêtes étant stockées par paires `(2k, 2k+1)`, l'index pair est toujours le canonique. La boucle itère avec `i += 2` — pas besoin de `std::set` pour éviter les doublons.

### for_each_edge : callable template

La signature utilise un template `Fn` au lieu de `std::function` pour éviter l'overhead d'indirection et l'éventuelle allocation dynamique :

```cpp
template <typename T, typename Fn>
void for_each_edge(const DCEL<T> &dcel, Fn &&fn);
```

N'importe quel callable est accepté : lambda, foncteur, pointeur de fonction.

### Export SVG : SVGContext

Les deux fonctions `write_svg` et `write_svg_dcel_debug` partagent un helper interne `SVGContext<T>` qui encapsule le calcul de la transformation monde → pixels et les sections communes (cellules, sites, rayons, bbox). Modifier l'apparence d'un élément commun ne nécessite qu'un seul changement.
