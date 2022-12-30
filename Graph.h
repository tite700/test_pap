#include <SDL2/SDL.h>
#include <stdio.h>
#include <ostream>
#include <iostream>
#include <vector>
#include <string>


class Graph {
 public:
  // Constructeur qui initialise la fenêtre SDL et la renderer
  Graph(int width, int height, const std::string& title);
  // Destructeur qui libère les ressources SDL
  ~Graph();

  // Méthode pour ajouter une nouvelle paire de valeurs x et y à la courbe
  void addPoint(float x, float y);
  // Méthode pour tracer la courbe à l'écran
  void draw();

 private:
  // Largeur et hauteur de la fenêtre
  int width_, height_;

  
  // Titre de la fenêtre
  std::string title_;
  // Fenêtre SDL et renderer
  SDL_Window* window_;
  SDL_Renderer* renderer_;
  // Séries de valeurs x et y
  std::vector<float> x_values_, y_values_;
};