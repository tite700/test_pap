#include "Graph.h"

// Implémentation du constructeur
Graph::Graph(int width, int height, const std::string& title) :
    width_(width), height_(height), title_(title) {
  // Initialisation de la bibliothèque SDL
  SDL_Init(SDL_INIT_VIDEO);
  // Création de la fenêtre
  window_ = SDL_CreateWindow(title.c_str(), SDL_WINDOWPOS_UNDEFINED,
                             SDL_WINDOWPOS_UNDEFINED, width, height,
                             SDL_WINDOW_SHOWN);
  // Création du renderer
  renderer_ = SDL_CreateRenderer(window_, -1, SDL_RENDERER_ACCELERATED);
}

// Implémentation du destructeur
Graph::~Graph() {
  // Libération des ressources SDL
  SDL_DestroyRenderer(renderer_);
  SDL_DestroyWindow(window_);
  SDL_Quit();
}

// Implémentation de la méthode pour ajouter une nouvelle paire de valeurs
void Graph::addPoint(float x, float y) {
  x_values_.push_back(x);
  y_values_.push_back(y);
}

// Implémentation de la méthode pour tracer la courbe à l'écran
void Graph::draw() {
    // Effacement de l'écran
    SDL_SetRenderDrawColor(renderer_, 255, 255, 255, 255);
    SDL_RenderClear(renderer_);

    // Tracé de la courbe
    SDL_SetRenderDrawColor(renderer_, 0, 0, 0, 255);
    for (size_t i = 0; i < x_values_.size() - 1; ++i) {
        int x1 = static_cast<int>(x_values_[i] * width_);
        int y1 = static_cast<int>(y_values_[i] * height_);
        int x2 = static_cast<int>(x_values_[i+1] * width_);
        int y2 = static_cast<int>(y_values_[i+1] * height_);
        SDL_RenderDrawLine(renderer_, x1, y1, x2, y2);
    }

    // Mise à jour de l'affichage
    SDL_RenderPresent(renderer_);

}


